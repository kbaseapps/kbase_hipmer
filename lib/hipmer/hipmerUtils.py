# -*- coding: utf-8 -*-
import os
import sys
import uuid
import numpy as np
from pprint import pformat

from Bio import SeqIO

from installed_clients.ReadsUtilsClient import ReadsUtils  # @IgnorePep8
from installed_clients.baseclient import ServerError
from installed_clients.ReadsAPIClient import ReadsAPI  # @IgnorePep8
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.kb_quastClient import kb_quast


class hipmerUtils:
    def __init__(self, config, token):
        #BEGIN_CONSTRUCTOR
        self.scratch = os.path.abspath(config['scratch'])
        self.callbackURL = os.environ.get('SDK_CALLBACK_URL')
        self.token = token

    # target is a list for collecting log messages
    def log(self, target, message):
        # we should do something better here...
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    def _validate_inputs(self, params):
        # do some basic checks
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'reads' not in params:
            raise ValueError('reads parameter is required')
        if 'output_contigset_name' not in params:
            raise ValueError('output_contigset_name parameter is required')
        # Check mer_sizes
        if 'mer_sizes' not in params:
            raise ValueError('mer_sizes is a required parameter')
        # Remove white spaces
        mer_sizes = params['mer_sizes'].replace(' ', '').replace('\t', '')
        mer_sizes_int = []
        mer_max = 0
        for mer in mer_sizes.split(','):
            try:
                meri = int(mer)
                mer_sizes_int.append(meri)
                if meri > mer_max:
                    mer_max = meri
            except:
                raise ValueError('mer sizes should be an integer')
            if meri < 10 or meri > 100:
                raise ValueError('mer sizes should be between 10 and 100')

        params['mer_sizes_int'] = mer_sizes_int
        params['mer_max'] = mer_max
        params['la_mer_max'] = mer_max + 2
        return True

    def check_reads(self, refs, console):
        # Hipmer requires some parameters to be set for the reads.
        # Let's check those first before wasting time with downloads.
        rapi = ReadsAPI(self.callbackURL, token=self.token,
                        service_ver='dev')
        err_msg = '%s does not specify a mean insert size\n'
        err_msg += 'Please re-upload the reads and used the advanced parameters'
        err_msg += 'options to specify the insert size parameters.\n'
        err_msg += 'This is required to run HipMer.'
        for ref in refs:
            p = {'workspace_obj_ref': ref}
            info = rapi.get_reads_info_all_formatted(p)
            if info['Type'] != 'Paired End':
                continue
            if info['Insert_Size_Mean'] == 'Not Specified' or \
               info['Insert_Size_Std_Dev'] == 'Not Specified':
                sys.stderr.write(err_msg % (info['Name']))
                return False

        return True

    def fix_pe_fq(self, fn):
        inf = open(fn, 'r')
        line = inf.readline()
        fields = line.split()
        rid = fields[0]
        if rid.startswith('@SRR') and rid.count('.') == 1 \
           and rid.find('/') < 0:
            print("Need to convert SRA fastq")
        else:
            inf.close()
            return
        os.rename(fn, fn + ".orig")
        out = open(fn, 'w')
        inf.seek(0)
        ln = 1
        pair = 1
        for line in inf:
            if ln == 1:
                name = line.split(' ')[0].replace('.', ':').rstrip()
                line = '%s/%s\n' % (name, pair)
                pair += 1
                if pair == 3:
                    pair = 1
            elif ln == 3:
                line = "+\n"
            out.write(line)
            ln += 1
            if ln == 5:
                ln = 1
        inf.close()
        out.close()

    def fixup_reads(self, params):
        for r in params['reads']:
            files_obj = params['readsfiles'][r['ref']]['files']
            if 'rev' not in files_obj or files_obj['rev'] is None:
                self.fix_pe_fq(files_obj['fwd'])

    def get_reads_RU(self, refs, console):
        readcli = ReadsUtils(self.callbackURL, token=self.token)

        typeerr = ('Supported types: KBaseFile.SingleEndLibrary ' +
                   'KBaseFile.PairedEndLibrary ' +
                   'KBaseAssembly.SingleEndLibrary ' +
                   'KBaseAssembly.PairedEndLibrary')
        try:
            reads = readcli.download_reads({'read_libraries': refs,
                                            'interleaved': 'true',
                                            'gzipped': None
                                            })['files']
        except ServerError as se:
            self.log(console, 'logging stacktrace from dynamic client error')
            self.log(console, se.data)
            if typeerr in se.message:
                prefix = se.message.split('.')[0]
                raise ValueError(
                    prefix + '. Only the types ' +
                    'KBaseAssembly.PairedEndLibrary ' +
                    'and KBaseFile.PairedEndLibrary are supported')
            else:
                raise

        self.log(console, 'Got reads data from converter:\n' + pformat(reads))
        return reads

    def generate_config(self, params):
        """
        Generate the HipMer config
        """
        self.config_file = '%s/%s' % (self.scratch, 'hipmer.config')
        with open(self.config_file, 'w') as f:
            # Describe the libraries ( one line per library )
            # lib_seq [ wildcard ][ prefix ][ insAvg ][ insSdev ][ avgReadLen ]
            #         [ hasInnieArtifact ][ isRevComped ][ useForContigging ]
            #         [ onoSetId ][ useForGapClosing ][ 5pWiggleRoom ]
            #         [3pWiggleRoom] [FilesPerPair] [ useForSplinting ]
            #
            # TODO: make these params
            #
            # FilesPerPair = 2
            fmt = 'lib_seq %s %s %d %d   %d %d %d   %d %d %d  %d %d %d %d\n'
            wd = '/kb/module/work/tmp'
            total_bases = 0
            for r in params['reads']:
                # TODO: check read type and set count
                files_obj = params['readsfiles'][r['ref']]['files']
                filelist = [files_obj['fwd'].replace(wd, '.')]
                if 'rev' in files_obj and files_obj['rev'] is not None:
                    rfile = files_obj['rev'].replace(wd, '.')
                    filelist.append(rfile)
                reads_obj = params['readsfiles'][r['ref']]
                print(files_obj)
                if files_obj['otype'] == 'single':
                    r['ins_avg'] = 0
                    r['ins_dev'] = 0
                    count = 0
                else:
                    r['ins_avg'] = int(reads_obj['insert_size_mean'])
                    r['ins_dev'] = int(reads_obj['insert_size_std_dev'])
                    count = len(filelist)
                r['avg_read_len'] = int(reads_obj['read_length_mean'])
                r['is_rev_comped'] = 0
                if 'read_orientation_outward' in reads_obj and \
                        reads_obj['read_orientation_outward'] is not None and \
                        reads_obj['read_orientation_outward'][0].lower() == 't':
                    r['is_rev_comped'] = 1

                files = ','.join(filelist)
                # lib_seq small.forward.fq,small.reverse.fq   small  215  10   \
                #    101 0 0      1 1 1  0 0 2 1
                f.write(fmt % (
                    files, r['prefix'], r['ins_avg'], r['ins_dev'],
                    0, r['has_innie_artifact'],
                    r['is_rev_comped'], r['use_for_contigging'],
                    r['ono_set_id'], r['use_for_gap_closing'],
                    r['fp_wiggle_room'], r['tp_wiggle_room'],
                    count, r['use_for_splinting']))
                total_bases += reads_obj['total_bases']
            f.write('\n')
            paramf = {
                'dynamic_min_depth': 'dynamic_min_depth %f\n',
                'mer_sizes': 'mer_sizes %s\n',
                'min_depth_cutoff': 'min_depth_cutoff %d\n',
                'gap_close_rpt_depth_ratio': 'gap_close_rpt_depth_ratio %f\n',
                'assm_scaff_len_cutoff': 'assm_scaff_len_cutoff %d\n',
                'mer_max': 'scaff_mer_size %d\n',
                'la_mer_max': 'la_mer_size_max %d\n'
            }
            if params['type'] == 'diploid':
                paramf['bubble'] = 'bubble_min_depth_cutoff %d\n'
                params['bubble'] = params['bubble_min_depth_cutoff']
                paramf['diploid'] = 'is_diploid %d\n'
                params['diploid'] = 1
            elif params['type'] == 'metagenome':
                for p in ['alpha', 'beta', 'tau', 'error_rate']:
                    paramf[p] = '%s %%f\n' % (p)
                    params[p] = params[p]
                paramf['bubble'] = 'bubble_min_depth_cutoff %d\n'
                params['bubble'] = params['bubble_min_depth_cutoff']
                paramf['metagenome'] = 'is_metagenome %d\n'
                params['metagenome'] = 1
            else:
                params['diploid'] = 0
                params['metagenome'] = 0

            for param in paramf:
                f.write(paramf[param] % (params[param]))
            f.close()

        return total_bases

    def generate_submit(self, tsize, debug=False):
        """
        Generate SLURM submit script
        """
        bpn = 500000000
        nodes = int((tsize + bpn - 1) / bpn)
        # It seems like hipmer fails with odd numbers of nodes
        # So let's add one if it is odd.
        if nodes % 2:
            nodes += 1
        self.submit = '%s/%s' % (self.scratch, 'slurm.submit')
        with open(self.submit, 'w') as f:
            f.write('#!/bin/bash\n')
            if debug:
                f.write('#SBATCH -q debug\n')
                f.write('#SBATCH --time=0:30:00\n')
            else:
                f.write('#SBATCH -q regular\n')
                f.write('#SBATCH --time=02:00:00\n')

            f.write('#SBATCH --nodes=%d\n' % nodes)
            f.write('#SBATCH --ntasks-per-node=32\n')
            f.write('#SBATCH --job-name=HipMer\n')
            f.write('export CORES_PER_NODE=${CORES_PER_NODE:=${SLURM_TASKS_PER_NODE%%\(*}}\n')
            f.write('N=${N:=${SLURM_NTASKS}}\n')
            f.write('HIPMER_INSTALL=${HIPMER_INSTALL:=$(pwd)/hipmer-v0.9.6}\n')
            f.write('INST=${HIPMER_INSTALL:=$1}\n')
            f.write('. $INST/env.sh\n')
            f.write('\n')
            f.write('export RUNDIR=${RUNDIR:=$(pwd)}\n')
            f.write('${INST}/bin/run_hipmer.sh ${RUNDIR}/hipmer.config\n')
            f.close()

    def save_assembly(self, wsname, output_contigs, token, name, console):
        self.log(console, 'Uploading FASTA file to Assembly')
        assemblyUtil = AssemblyUtil(self.callbackURL, token=token,
                                    service_ver='dev')
        assemblyUtil.save_assembly_from_fasta({'file': {'path': output_contigs},
                                               'workspace_name': wsname,
                                               'assembly_name': name
                                               })

    def prepare_run(self, params):
        """
        run HPC enabled hipmer
        Returns: report object
        """
        console = []
        self.log(console, 'Running run_hipmer_hpc with params=')
        self.log(console, pformat(params))

        # Validate parameters.  This will raise an error if there
        # is a problem.
        self._validate_inputs(params)
        ws_name = params['workspace_name']

        # Get the read library
        print("Running pre stage")
        refs = []
        for read in params['reads']:
            read_name = read['read_library_name']
            if '/' in read_name:
                ref = read_name
            else:
                ref = ws_name + '/' + read_name
            refs.append(ref)
            read['ref'] = ref
        if not self.check_reads(refs, console):
            raise ValueError('The reads failed validation\n')

        params['readsfiles'] = self.get_reads_RU(refs, console)
        self.fixup_reads(params)

        # Generate submit script
        ts = self.generate_config(params)
        debug = False
        if 'usedebug' in params and params['usedebug'] > 0:
            debug = True
        self.generate_submit(ts, debug=debug)

    def finish_run(self, params):
        """
        Finish up the run by uploading output and
        creating the report
        """
        console = []
        self.log(console, 'Running post')

        # run hipmer, capture output as it happens
        self.log(console, 'running hipmer:')

        output_contigs = os.path.join(self.scratch, 'results', 'final_assembly.fa')
        output_name = params['output_contigset_name']
        if not os.path.exists(output_contigs):
            print("It looks like HipMER failed for some reason.")
            print("Show errors in log file")
            logfile = ''
            for fn in os.listdir('.'):
                if fn.startswith('slurm-'):
                    logfile = fn
            if logfile != '':
                with open(logfile, 'r') as f:
                    for line in f:
                        if line.lower().find('error') >= 0:
                            print(line)
            raise RuntimeError("Error in HipMER execution")

        wsname = params['workspace_name']
        self.log(console, 'Uploading FASTA file to Assembly')
        assemblyUtil = AssemblyUtil(self.callbackURL, token=self.token)
        save_input = {'file': {'path': output_contigs},
                      'workspace_name': wsname,
                      'assembly_name': output_name
                      }
        output_data_ref = assemblyUtil.save_assembly_from_fasta(save_input)
        # create a Report
        # compute a simple contig length distribution for the report
        lengths = []
        for seq_record in SeqIO.parse(output_contigs, 'fasta'):
            lengths.append(len(seq_record.seq))

        report = ''
        report += 'ContigSet saved to: ' + params['workspace_name'] + '/'
        report += params['output_contigset_name'] + '\n'
        report += 'Assembled into ' + str(len(lengths)) + ' contigs.\n'
        report += 'Avg Length: ' + str(sum(lengths) / float(len(lengths))) + ' bp.\n'

        bins = 10
        counts, edges = np.histogram(lengths, bins)
        report += 'Contig Length Distribution (# of contigs -- min to max basepairs):\n'
        for c in range(bins):
            report += '   \%d\t--\t%d' % (counts[c], edges[c])
            report += ' to %d bp\n' % (edges[c + 1])

        print('Running QUAST')
        kbq = kb_quast(self.callbackURL)
        try:
            quastret = kbq.run_QUAST({'files': [{'path': output_contigs,
                                                 'label': params['output_contigset_name']}]})
        except Exception as e:
            # not really any way to test this, all inputs have been checked
            # earlier and should be ok
            print('Logging exception from running QUAST')
            print((str(e)))
            # TODO delete shock node
            raise

        print('Saving report')
        kbr = KBaseReport(self.callbackURL)
        try:
            report_info = kbr.create_extended_report(
                {'message': report,
                 'objects_created': [{'ref': output_data_ref,
                                      'description': 'Assembled contigs'}],
                 'direct_html_link_index': 0,
                 'html_links': [{'shock_id': quastret['shock_id'],
                                 'name': 'report.html',
                                 'label': 'QUAST report'}
                                ],
                 'report_object_name': 'kb_megahit_report_' + str(uuid.uuid4()),
                 'workspace_name': params['workspace_name']
                 })
        except Exception as e:
            # not really any way to test this, all inputs have been checked earlier and should be
            # ok
            print('Logging exception from creating report object')
            print((str(e)))
            # TODO delete shock node
            raise

        # STEP 6: contruct the output to send back
        output = {'report_name': report_info['name'],
                  'report_ref': report_info['ref']
                  }
        return output

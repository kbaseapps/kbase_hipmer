# -*- coding: utf-8 -*-
import os
import sys
import uuid
import numpy as np
from pprint import pformat
from pprint import pprint
import subprocess

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
        print("PARAMS \n{}".format(params))

        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'reads' not in params:
            raise ValueError('reads parameter is required')
        if 'output_contigset_name' not in params:
            raise ValueError('output_contigset_name parameter is required')
        if 'read_type' not in params['reads'][0]:
            raise ValueError('read_type parameter is required')
        if 'is_rev_comped' not in params['reads'][0]:
            raise ValueError('is_rev_comped parameter is required')
        if params['is_meta'] is not None:
            if not params['is_meta']['min_depth_cutoff']:
                raise ValueError('If this is a metagenome, there needs to be a cutoff value set for min_depth_cutoff')

        # Check mer_sizes
        if 'mer_sizes' not in params:
            raise ValueError('mer_sizes is a required parameter')

        # Parse the kmer string into a list.  This step verifies that there
        # were integers separated by commas. I don't actually use the created list (mer_sizes_int).
        mer_sizes = params['mer_sizes'].replace(' ', '').replace('\t', '')
        mer_sizes_int = []
        for mer in mer_sizes.split(','):
            try:
                meri = int(mer)
                mer_sizes_int.append(meri)
            except:
                raise ValueError('mer sizes should be an integer')

            if meri < 10 or meri > 127:
                raise ValueError('mer sizes should be between 10 and 127.')

        # just in case I need the list someday.
        params['mer_sizes_int'] = mer_sizes_int

        # I check that the user input insert sizes and insert size standard deviation
        # if they are using paired-end reads.
        # However, this is checked in the "check_reads" function.


        return True

    def check_reads(self, refs, console, params):
        # Hipmer requires some parameters to be set for the reads.
        # Let's check those first before wasting time with downloads.
        # TODO: Shane gets read info from read object.  Jeff will check
        #  that there were insert sizes input by the user by checking params
        # rapi = ReadsAPI(self.callbackURL, token=self.token,
        #                 service_ver='dev')
        # err_msg = '%s does not specify a mean insert size\n'
        # err_msg += 'Please re-upload the reads and use the advanced parameters'
        # err_msg += 'options to specify the insert size parameters.\n'
        # err_msg += 'This is required to run HipMer.'
        # for ref in refs:
        #     p = {'workspace_obj_ref': ref}
        #     info = rapi.get_reads_info_all_formatted(p)
        #     if info['Type'] != 'Paired End':
        #         continue
        #     if info['Insert_Size_Mean'] == 'Not Specified' or \
        #        info['Insert_Size_Std_Dev'] == 'Not Specified':
        #         sys.stderr.write(err_msg % (info['Name']))
        #         return False
        # TODO: Jeffs section using params to check for insert sizes if paired-end reads.
        for r in params['reads']:
            if r['read_type'] == 'paired':
                # TODO: try and get read name from read ref
                # # get read names in case we need to write error msg
                # read_ref = r['read_library_name']
                # rapi = ReadsAPI(self.callbackURL, token=self.token, service_ver='dev')
                # p = {'workspace_obj_ref': r}
                # info = rapi.get_reads_info_all_formatted(p)
                # read_name = info['Name']

                if r['ins_avg'] is not None and r['ins_dev'] is not None:
                    return True
                else:
                    err_msg = 'It looks like the user did not specify an average read insert size and/or the insert size stdev '
                    err_msg += 'which is required for paired-end reads.'
                    err_msg += 'The offending reads are: %s'
                    sys.stderr.write(err_msg % (r['read_library_name']))
                    return False

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



    def get_total_bases(self, params):
        """
        The total number of bases for all sequence files will be used to estimate the number of
        nodes required to run hipmer.
        """
        total_bases=0
        for r in params['reads']:
            # TODO: check read type and set count
            reads_obj = params['readsfiles'][r['ref']]
            total_bases += reads_obj['total_bases']
        return total_bases


    def generate_command(self, params, nodes):
        """
        SBATCH will run this command in the "generate_submit" function
        An example command:

        hipmer --threads=$((4 * 32)) \
               -k 21 \
               --interleaved frags.25K.fastq: i180:s10 \
               --interleaved jumps.25K.fastq: rc:i3000:s500 \
               --min-depth 7
        """
        iterator=0
        final_read_args=''
        hipmer_command=''

        for r in params['reads']:
            kmer_str = params['mer_sizes']
            read_ref = r['ref']

            # we are not running the command in the docker container so the path to the fastq
            # needs to be pointing to somewhere outside and not /kb/module/work/tmp
            ori_file_name = params['readsfiles'][read_ref]['files']['fwd']
            file_name = os.path.basename(ori_file_name)

            # build base command
            hipmer_command = "hipmer --threads=$(({} * 32)) -k {}".format(nodes, params['mer_sizes'])

            #
            # Test if we have a metagenome
            #
            min_depth=None
            metagenome_opts=''
            plant_opts=''

            if params['is_meta'] is not None:
                #
                # if metagenome
                #
                min_depth = params['is_meta']['min_depth_cutoff']
                if params['is_meta']['aggressive']:
                    # if metagenome and aggressive algorithm should be used
                    metagenome_opts = "--aggressive --meta --min-depth {} ".format(min_depth)
                else:
                    metagenome_opts = "--meta --min-depth {} ".format(min_depth)

            #
            # Test if we have plant data
            #
            if params['is_plant'] is not None:
                diploid_value = params['is_plant']['diploid']
                is_high_heterozygosity = params['is_plant']['high_heterozygosity']
                bubble_min_depth_cutoff = params['is_plant']['bubble_min_depth_cutoff']

                plant_opts = "--bubble-depth-cutoff {} --diploidy {} ".format(bubble_min_depth_cutoff, diploid_value)

            #
            # format argument for "reads" section
            #
            if r['read_type'] == "paired":
                # if paired end reads
                # test if reverse compliment ("outtie") or "innie" library
                if r['is_rev_comped'] == 'outtie':
                    read_args = " --interleaved {}:rc:i{}:s{}".format(file_name, r['ins_avg'], r['ins_dev'])
                else:
                    read_args = " --interleaved {}:i{}:s{}".format(file_name, r['ins_avg'], r['ins_dev'])
            else:
                if r['is_rev_comped'] == 'outtie':
                    read_args = " --single jumps.{}.fastq rc".format(iterator)
                else:
                    read_args = " --single single.{}.fastq".format(iterator)

            final_read_args += read_args

        hipmer_command += metagenome_opts + plant_opts + final_read_args


        print("HIPMER COMMAND : {}".format(hipmer_command))

        return hipmer_command


    def generate_submit(self, total_bases, params, debug=False):
        """
        Generate SLURM submit script
        """
        bpn = 500000000
        nodes = int((total_bases + bpn - 1) / bpn)


        # It seems like hipmer fails with odd numbers of nodes
        # So let's add one if it is odd.
        if nodes % 2:
            nodes += 1

        hipmer_command = self.generate_command(params,nodes)

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
            f.write('#SBATCH -C haswell\n')
            f.write('#SBATCH --ntasks-per-node=32\n')
            f.write('#SBATCH --job-name=HipMer\n')
            f.write('#SBATCH --license=SCRATCH\n')
            f.write('#SBATCH -o slurm.out\n')
            f.write('set -e\n')
            f.write('source /etc/profile.d/modules.sh\n')
            f.write('export CRAY_CPU_TARGET=haswell\n')
            f.write('printenv\n')
            f.write('module switch bupc-narrow bupc-narrow')
            f.write('export CORES_PER_NODE=${CORES_PER_NODE:=${SLURM_TASKS_PER_NODE%%\(*}}\n')
            f.write('export THREADS=${THREADS:=${SLURM_NTASKS}}\n')
            f.write('echo "Detected CORES_PER_NODE=${CORES_PER_NODE} and THREADS=${THREADS}"\n')
            f.write('echo Executing ' + hipmer_command + '\n')
            f.write('ls\n')
            f.write('echo "at $(date) on $(uname -n)"\n')
            f.write('HIPMER_INSTALL=$(pwd)/v1.0-378*/bin\n')
            f.write('${HIPMER_INSTALL}/' + hipmer_command + '\n')
            f.close()

            return self.submit

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
        self.log(console, 'Running run_kbase_hipmer with params=')
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
        if not self.check_reads(refs, console, params):
            raise ValueError('The reads failed validation\n')

        params['readsfiles'] = self.get_reads_RU(refs, console)
        self.fixup_reads(params)

        # Generate submit script
        (total_bases) = self.get_total_bases(params)


        debug = False
        if 'usedebug' in params and params['usedebug'] > 0:
            debug = True
        submit_file = self.generate_submit(total_bases, params, debug=debug)


    def finish_run(self, params):
        """
        Finish up the run by uploading output and
        creating the report
        """
        console = []
        self.log(console, 'Running post')

        # run hipmer, capture output as it happens
        self.log(console, 'running hipmer:')

        print("RUNDIR {}".format(os.environ['RUNDIR']))
        test = os.environ['RUNDIR']
        testid = test + "/hipmer-run-*/results/final_assembly.fa"
        output_contigs = os.path.join(os.environ['RUNDIR'], "hipmer-run-*/results/final_assembly.fa")
        print("OUTPUT CONTIGS: {}".format(output_contigs))

        if os.path.exists(testid):
            print("OUTPUT CONTIGS ARE HERE {}".format(testid))
        else:
            print("testid not found")

        output_name = params['output_contigset_name']
        slurm_out = os.path.join(self.scratch, 'slurm.out')
        if not os.path.exists(output_contigs):
            self.log(console, "It looks like HipMER failed. Could not find the output contigs.")
            self.log(console, "Show errors in log file")
            with open(slurm_out, 'r') as f:
                for line in f:
                    if line.lower().find('error') >= 0:
                        self.log(console, line)
            raise RuntimeError("Error in HipMER execution")

#        sys.exit()
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

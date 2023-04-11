# -*- coding: utf-8 -*-
import os
import sys
import uuid
import glob
import numpy as np
from pprint import pformat
from pprint import pprint
import subprocess
from Bio import SeqIO

from installed_clients.ReadsUtilsClient import ReadsUtils  # @IgnorePep8
from installed_clients.baseclient import ServerError
from installed_clients.ReadsAPIServiceClient import ReadsAPI  # @IgnorePep8
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.kb_quastClient import kb_quast
from installed_clients.specialClient import special

class hipmerUtils:
    def __init__(self, config, token):
        #BEGIN_CONSTRUCTOR
        self.scratch = os.path.abspath(config['scratch'])
        self.callbackURL = os.environ.get('SDK_CALLBACK_URL')
        print(config['service-wizard'])
        self.rapi = ReadsAPI(config['service-wizard'], token=token)
        self.readcli = ReadsUtils(self.callbackURL, token=token)
        # rapi = ReadsAPI(self.srv_wizard, token=self.token, service_ver='dev')
        self.sr = special(self.callbackURL, token=token)
        self.submit_script = 'slurm2.sl'
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

        # Check mer_sizes
        if 'mer_sizes' not in params:
            raise ValueError('mer_sizes is a required parameter. Recommended vals for metagenome: 21,33,55,77,99')
        # Check scaff_mer_lens
        if 'scaff_mer_lens' not in params:
            raise ValueError('scaff_mer_lens is a required parameter. Recommended vals for metagenome: 99,33')
            
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

        # Parse the scaff kmer string into a list.  This step verifies that there
        # were integers separated by commas. I don't actually use the created list (scaff_mer_sizes_int).
        if params.get('scaff_mer_lens'):
            scaff_mer_lens = params['scaff_mer_lens'].replace(' ', '').replace('\t', '')
            scaff_mer_lens_int = []
            for mer in scaff_mer_lens.split(','):
                try:
                    meri = int(mer)
                    scaff_mer_lens_int.append(meri)
                except:
                    raise ValueError('scaff mer sizes should be an integer')

                if meri <= 10:
                    raise ValueError('scaff mer sizes should be above 10.')

            # just in case I need the list someday.
            params['scaff_mer_lens_int'] = scaff_mer_lens_int

        
        # I check that the user input insert sizes and insert size standard deviation
        # if they are using paired-end reads.
        # However, this is checked in the "check_reads" function.


        return True


    # _validate_input_reads_sizelimit()
    #
    def _validate_input_reads_sizelimit (self, refs, console, params):
        # Make sure user isn't trying to launch a job that's too big
        self.log(console, "GOT TO A")
        if params.get('read_Gbp_limit'):
            self.log(console, "GOT TO B")
            read_Gbp_limit = int(params['read_Gbp_limit']) * 1000000000
            self.log(console, "GOT TO C")

            
            total_bases = 0
            for ref in refs:
                self.log(console, "GOT TO D")
                this_reads_metadata = self.rapi.get_reads_info_all_formatted ({'workspace_obj_ref':ref})
                self.log(console, "GOT TO E")
                these_bases = this_reads_metadata['Total_Number_of_Bases'].replace(',','')
                total_bases += int(these_bases)
                self.log(console, "GOT TO F")

            if total_bases > read_Gbp_limit:
                self.log(console, "GOT TO G")
                err_msg = "ERROR: reads size exceeds limit for running MetaHipMer.  Input reads total bp {} > {}".format(total_bases, read_Gbp_limit)
                err_msg += "\nYou can either reduce the number of libraries in your input or increase the limit in the advanced settings.  If you choose the latter, please get approval from KBase support at http://www.kbase.us/support prior to submitting your job"
                self.log(console, err_msg)
                raise ValueError (err_msg)
            else:
                self.log(console, "GOT TO H")
                success_msg = "reads size does not exceed limit for running MetaHipMer.  Input reads total bp {} <= {}".format(total_bases, read_Gbp_limit)
                self.log(console, success_msg)
        self.log(console, "GOT TO I")
        return True
    
            
    def check_reads(self, refs, console, params):
        # Hipmer requires some parameters to be set for the reads.
        # Let's check those first before wasting time with downloads.
        # We will check whether the reads are paired end or single end from the read object.
        err_msg = '%s must be a paired-end read library.\n'
        err_msg += 'This is required to run HipMer.'
        index = 0

        print("REFS \n{}".format(refs))

        for ref in refs:
            p = {'workspace_obj_ref': ref}

            if info['Type'] == 'Paired End':
                # TODO: someday we should be able to handle single end reads
                params['reads'][index]['read_type'] = 'paired'
                params['reads'][index]['info'] = info
            elif info['Type'] == 'Single End':
                params['reads']['read_type'] = 'single'
                continue
            else:
                sys.stderr.write(err_msg % (info['Name']))
                return False

            index+=1

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
            files_obj = params['readsfiles'][r]['files']
            if 'rev' not in files_obj or files_obj['rev'] is None:
                self.fix_pe_fq(files_obj['fwd'])

    def get_reads_RU(self, refs, console):

        typeerr = ('Supported types: KBaseFile.SingleEndLibrary ' +
                   'KBaseFile.PairedEndLibrary ' +
                   'KBaseAssembly.SingleEndLibrary ' +
                   'KBaseAssembly.PairedEndLibrary')
        try:
            reads = self.readcli.download_reads({'read_libraries': refs,
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
            reads_obj = params['readsfiles'][r]
            total_bases += reads_obj['total_bases']
        return total_bases


    def get_total_gigs(self, params):
        """
        The total number of bases for all sequence files will be used to estimate the number of
        nodes required to run hipmer.
        """
        total_size_gigs=0
        for r in params['reads']:
            # we are not running the command in the docker container so the path to the fastq
            # needs to be pointing to somewhere outside and not /kb/module/work/tmp
            fastq_path = params['readsfiles'][r]['files']['fwd']

            fsize = os.stat(fastq_path).st_size

            size_gigs = int(fsize)/1024/1024/1024
            print("file size: {} {}".format(fastq_path, size_gigs))

            total_size_gigs += size_gigs

        return total_size_gigs


    def generate_command(self, params, nodes):
        """
        SBATCH will run this command in the "generate_submit" function
        """
        iterator=0
        final_read_args=''
        hipmer_command=''

        kmer_str = params['mer_sizes']
        scaff_kmer_str = None
        if params.get('scaff_mer_lens'):
            scaff_kmer_str = params['scaff_mer_lens']
            
        # Test if we have a metagenome
        #

        for read_ref in params['reads']:

            # we are not running the command in the docker container so the path to the fastq
            # needs to be pointing to somewhere outside and not /kb/module/work/tmp
            read_file = params['readsfiles'][read_ref]
            ori_file_name = read_file['files']['fwd']
            file_name = os.path.basename(ori_file_name)

            #
            # format argument for "reads" section
            #
            read_args = " -r {}".format(file_name)
            #if r['read_type'] == "paired":
            #    # if paired end reads
            #    # test if reverse compliment ("outtie") or "innie" library
            #    read_args = "-r {}".format(file_name)
            #else:
            #    print("Ignoring %s" % (file_name))
                    

            final_read_args += read_args

        # build base command
        hipmer_command = "mhm2.py --kmer-lens {} ".format(kmer_str)
        if scaff_kmer_str is not None:
            hipmer_command += "--scaff-kmer-lens {} ".format(scaff_kmer_str) 
        hipmer_command += final_read_args

        return hipmer_command


    def generate_submit(self, total_size_gigs, params, debug=False):
        """
        Generate SLURM submit script
        """
        # the formula for estimating number of nodes required
        # nodes = 40*Gb-reads/80G
        nodes = round( ((total_size_gigs * 40) / 30) + 1 )
        if nodes < 1:
            nodes = 2 # can't be an odd number

        # It seems like hipmer fails with odd numbers of nodes
        # So let's add one if it is odd.
        if nodes % 2:
            nodes += 1

        hipmer_command = self.generate_command(params,nodes)
        print("HIPMER CMD: {}".format(hipmer_command))

        submit = '%s/%s' % (self.scratch, self.submit_script)
        with open(submit, 'w') as f:
            f.write('#!/bin/bash\n')
            if debug:
                f.write('#SBATCH -q debug\n')
                f.write('#SBATCH --time=0:30:00\n')
            else:
                f.write('#SBATCH -q regular\n')
                f.write('#SBATCH --time=02:00:00\n')

            f.write('#SBATCH --nodes=%d\n' % nodes)
            f.write('#SBATCH -C knl,quad,cache\n')
            f.write('#SBATCH --ntasks-per-node=68\n')
            f.write('#SBATCH --job-name=HipMer\n')
            f.write('#SBATCH --license=SCRATCH\n')
            f.write('set -e\n\n')
            f.write('if [ -z "$MODULEPATH" ] ; then\n')
            f.write('    . /etc/profile.d/modules.sh\n')
            f.write('fi\n')
            f.write('module load PrgEnv-cray/6.0.10\n')
            f.write('module load upcxx/2022.3.0\n')
            f.write('HIPMER_INSTALL=$(pwd)/v2*/bin\n')
            f.write('${HIPMER_INSTALL}/' + hipmer_command + '\n')
            f.close()

            return submit

    def fasta_filter_contigs_generator(self, fasta_record_iter, min_contig_length):
        """ generates SeqRecords iterator for writing from a legacy contigset object """
        rows = 0
        rows_added = 0
        for record in fasta_record_iter:
            rows += 1
            if len(record.seq) >= min_contig_length:
                rows_added += 1
                yield record
        #print(f' - filtered out {rows - rows_added} of {rows} contigs that were shorter than {(min_contig_length)} bp.')
        print(' - filtered out '+str(rows - rows_added)+' of '+str(rows)+' contigs that were shorter than ('+str(min_contig_length)+') bp.')

    def filter_contigs_by_length(self, fasta_file_path, min_contig_length):
        """ removes all contigs less than the min_contig_length provided """
        filtered_fasta_file_path = os.path.abspath(fasta_file_path).split('.fa')[0] + "_filtered.fa"

        fasta_record_iter = SeqIO.parse(fasta_file_path, 'fasta')
        SeqIO.write(self.fasta_filter_contigs_generator(fasta_record_iter, min_contig_length),
                    filtered_fasta_file_path, 'fasta')

        return filtered_fasta_file_path

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
        for read_name in params['reads']:
            # read_name = read['read_library_name']
            if '/' in read_name:
                ref = read_name
            else:
                ref = ws_name + '/' + read_name
            refs.append(ref)

        #if not self.check_reads(refs, console, params):
        #    raise ValueError('The reads failed validation\n')

        # make sure not too big to run
        self.log(console, "PRE VALIDATE")
        self._validate_input_reads_sizelimit (refs, console, params)
        self.log(console, "POST VALIDATE")

        # download reads to filesystem
        params['readsfiles'] = self.get_reads_RU(refs, console)
        # self.fixup_reads(params)

        # Generate submit script
        (total_size_gigs) = self.get_total_gigs(params)


        debug = False
        if 'usedebug' in params and params['usedebug'] > 0:
            debug = True
        submit_file = self.generate_submit(total_size_gigs, params, debug=debug)
        return submit_file

    def submit(self):
        p = {'submit_script': self.submit_script}
        res = self.sr.slurm(p)
        print('slurm:'+str(res))

    def finish_run(self, params):
        """
        Finish up the run by uploading output and
        creating the report
        """
        console = []
        self.log(console, 'Running post')

        # run hipmer, capture output as it happens
        self.log(console, 'running hipmer:')

        # grab path of output contigs
        output_contigs = ''
        for root, subdirs, files in os.walk(self.scratch):
            for f in files:
                if f == 'final_assembly.fasta':
                    output_contigs = os.path.join(root,f)
                    print("found OUTPUT CONTIGS {}".format(output_contigs))
                    continue

        output_name = params['output_contigset_name']

        if not os.path.exists(output_contigs):
            self.log(console, "It looks like HipMER failed. Could not find the output contigs.")
            self.log(console, "Check the log for errors")
            raise RuntimeError("Error in HipMER execution")

        wsname = params['workspace_name']

        self.log(console, 'Filtering short length contigs from HipMer assembly')

        assemblyUtil = AssemblyUtil(self.callbackURL, token=self.token)

        assembly_size_filter = params['assembly_size_filter']

        filtered_fasta_file_path = self.filter_contigs_by_length(output_contigs, assembly_size_filter)

        if os.stat(filtered_fasta_file_path).st_size == 0:
            raise ValueError("Error: Using input parameters, you have filtered all contigs from the HipMer \
                             assembly. Decrease the minimum contig size and try again.")
        else:
            output_contigs = filtered_fasta_file_path

        self.log(console, 'Uploading FASTA file to Assembly')

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
            # TODO delete shock node
            raise e

        # STEP 6: contruct the output to send back
        output = {'report_name': report_info['name'],
                  'report_ref': report_info['ref']
                  }
        return output

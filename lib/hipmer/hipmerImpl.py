#BEGIN_HEADER
import os
import sys
import hashlib
import subprocess
import requests
import re
import traceback
import uuid
import numpy as np
from pprint import pformat

from Bio import SeqIO

from biokbase.workspace.client import Workspace as workspaceService
from ReadsUtils.ReadsUtilsClient import ReadsUtils  # @IgnorePep8
from ReadsUtils.baseclient import ServerError
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from KBaseReport.KBaseReportClient import KBaseReport
from KBaseReport.baseclient import ServerError as _RepError
from kb_quast.kb_quastClient import kb_quast
from kb_quast.baseclient import ServerError as QUASTError
#from kb_ea_utils.kb_ea_utilsClient import kb_ea_utils

import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()
#END_HEADER


class hipmer:
    '''
    Module Name:
    hipmer

    Module Description:
    A KBase module: hipmer
    '''

    # WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    VERSION = "0.0.1"
    GIT_URL = "git@github.com:kbaseapps/kbase_hipmer.git"
    GIT_COMMIT_HASH = "cf9790a0f15a0f8a08a6e1a654af969bd072cb5a"

    #BEGIN_CLASS_HEADER
    workspaceURL = None

    # target is a list for collecting log messages
    def log(self, target, message):
        # we should do something better here...
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    def get_pe_library_deinterleaved(self, ws_data, ws_info, forward, reverse):
        pass

    def get_reads_RU(self, ctx, refs, console):
        readcli = ReadsUtils(self.callbackURL, token=ctx['token'],
                             service_ver='dev')

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

    def get_reads(self, ctx, ref, console):
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': ref}])
            data = objects[0]['data']
            info = objects[0]['info']
            type_name = info[2].split('.')[1].split('-')[0]
        except Exception as e:
            raise ValueError(
                'Unable to fetch read library object from workspace: ' +
                str(e))
            # to get the full stack trace: traceback.format_exc()

        # Download the paired end library
        if type_name == 'PairedEndLibrary':
            try:
                if 'lib1' in data:
                    forward_reads = data['lib1']['file']
                elif 'handle_1' in data:
                    forward_reads = data['handle_1']
                if 'lib2' in data:
                    reverse_reads = data['lib2']['file']
                elif 'handle_2' in data:
                    reverse_reads = data['handle_2']
                else:
                    reverse_reads = {}

                fr_file_name = forward_reads['id']
                if 'file_name' in forward_reads:
                    fr_file_name = forward_reads['file_name']

                forward_reads_file_loc = os.path.join(self.scratch,
                                                      fr_file_name)
                forward_reads_file = open(forward_reads_file_loc, 'w', 0)
                self.log(console, 'downloading reads file: ' +
                         str(forward_reads_file_loc))
                headers = {'Authorization': 'OAuth ' + ctx['token']}
                url = forward_reads['url'] + '/node/' + forward_reads['id']
                url += '?download'
                r = requests.get(url, stream=True, headers=headers)
                for chunk in r.iter_content(1024):
                    forward_reads_file.write(chunk)
                forward_reads_file.close()
                self.log(console, 'done')
                # END NOTE

                if 'interleaved' in data and data['interleaved']:
                    self.log(console, 'extracting forward/reverse reads into separate files')
                    if re.search('gz', fr_file_name, re.I):
                        bcmdstring = 'gunzip -c ' + forward_reads_file_loc
                    else:
                        bcmdstring = 'cat ' + forward_reads_file_loc

                    cmdstring = bcmdstring + '| (paste - - - - - - - -  | '
                    cmdstring += 'tee >(cut -f 1-4 | tr "\t" "\n" > '
                    cmdstring += self.scratch
                    cmdstring += '/forward.fastq) | cut -f 5-8 | '
                    cmdstring += 'tr "\t" "\n" > '
                    cmdstring += self.scratch + '/reverse.fastq )'
                    cmdProcess = subprocess.Popen(cmdstring,
                                                  stdout=subprocess.PIPE,
                                                  stderr=subprocess.PIPE,
                                                  shell=True,
                                                  executable='/bin/bash')
                    stdout, stderr = cmdProcess.communicate()
                    message = "cmdstring: " + cmdstring + " stdout: "
                    message += stdout + " stderr: " + stderr
                    self.log(console, message)

                    fr_file_name = 'forward.fastq'
                    forward_reads['file_name'] = fr_file_name
                    rev_file_name = 'reverse.fastq'
                    reverse_reads['file_name'] = rev_file_name
                    reads = [fr_file_name, rev_file_name]
                else:
                    # we need to read in reverse reads file separately
                    rev_file_name = reverse_reads['id']
                    if 'file_name' in reverse_reads:
                        rev_file_name = reverse_reads['file_name']
                    # NOTE: Replace with local method call
                    reverse_reads_file_loc = os.path.join(self.scratch,
                                                          rev_file_name)
                    reverse_reads_file = open(reverse_reads_file_loc, 'w', 0)
                    message = 'downloading reverse reads file: '
                    message += str(reverse_reads_file_loc)
                    self.log(console, message)
                    url = reverse_reads['url'] + '/node/' + reverse_reads['id']
                    url += '?download'
                    r = requests.get(url, stream=True, headers=headers)
                    for chunk in r.iter_content(1024):
                        reverse_reads_file.write(chunk)
                    reverse_reads_file.close()
                    reads = [fr_file_name, rev_file_name]
                    self.log(console, 'done')
                    # END NOTE
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to download paired-end read library files: ' + str(e))
        else:
            raise ValueError('Cannot yet handle library type of: ' + type_name)

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
                r['ins_avg'] = int(reads_obj['insert_size_mean'])
                r['ins_dev'] = int(reads_obj['insert_size_std_dev'])
                r['avg_read_len'] = int(reads_obj['read_length_mean'])
                r['is_rev_comped'] = 0
                if 'read_orientation_outward' in reads_obj and \
                        reads_obj['read_orientation_outward'][0].lower() == 't':
                    r['is_rev_comped'] = 1

                count = len(filelist)
                files = ','.join(filelist)
                # lib_seq small.forward.fq,small.reverse.fq   small  215  10   \
                #    101 0 0      1 1 1  0 0 2 1
                f.write(fmt % (
                    files, r['prefix'], r['ins_avg'], r['ins_dev'],
                    r['avg_read_len'], r['has_innie_artifact'],
                    r['is_rev_comped'], r['use_for_contigging'],
                    r['ono_set_id'], r['use_for_gap_closing'],
                    r['fp_wiggle_room'], r['tp_wiggle_room'],
                    count, r['use_for_splinting']))
                total_bases += reads_obj['total_bases']
            f.write('\n')
            paramf = {
                'is_diploid': 'is_diploid %d\n',
                'dynamic_min_depth': 'dynamic_min_depth %f\n',
                'mer_size': 'mer_size %d\n',
                'min_depth_cutoff': 'min_depth_cutoff %d\n',
                'gap_close_rpt_depth_ratio': 'gap_close_rpt_depth_ratio %f\n',
                'assm_scaff_len_cutoff': 'assm_scaff_len_cutoff %d\n'
            }
            if params['is_diploid'] == 1:
                paramf['bubble_min_depth_cutoff'] = 'bubble_min_depth_cutoff %d'
            for param in paramf:
                f.write(paramf[param] % (params[param]))
            f.close()

        return total_bases

    def generate_submit(self, tsize):
        """
        Generate SLURM submit script
        """
        bpn = 1000000000
        nodes = int((tsize + bpn - 1) / bpn)
        self.submit = '%s/%s' % (self.scratch, 'slurm.submit')
        with open(self.submit, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH --partition=debug\n')
            f.write('#SBATCH --nodes=%d -C haswell\n' % nodes)
            f.write('#SBATCH --ntasks-per-node=32\n')
            f.write('#SBATCH --time=00:30:00\n')
            f.write('#SBATCH --job-name=HipMer\n')
            f.write('export CORES_PER_NODE=${CORES_PER_NODE:=${SLURM_TASKS_PER_NODE%%\(*}}\n')
            f.write('N=${N:=${SLURM_NTASKS}}\n')
            f.write('HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/hipmer-install-cori}\n')
            f.write('INST=${HIPMER_INSTALL:=$1}\n')
            f.write('. $INST/env.sh\n')
            f.write('\n')
            f.write('export RUNDIR=${RUNDIR:=$(pwd)}\n')
            f.write('${INST}/bin/run_hipmer.sh ${RUNDIR}/hipmer.config\n')
            f.close()

    def get_wsid(self, ws_name, token):
        ws = workspaceService(self.workspaceURL, token=token)
        wsinfo = ws.get_workspace_info({'workspace': ws_name})
        return wsinfo[0]

    def save_assembly_old(self, output_contigs, ctx, params, console, ws, ws_name):
        wsid = self.get_wsid()

        # parse the output and save back to KBase

        # Warning: this reads everything into memory!  Will not work if
        # the contigset is very large!
        contigset_data = {
            'id': 'hipmer.contigset',
            'source': 'User assembled contigs from reads in KBase',
            'source_id': 'none',
            'md5': 'md5 of what? concat seq? concat md5s?',
            'contigs': []
        }

        lengths = []
        for seq_record in SeqIO.parse(output_contigs, 'fasta'):
            contig = {
                'id': seq_record.id,
                'name': seq_record.name,
                'description': seq_record.description,
                'length': len(seq_record.seq),
                'sequence': str(seq_record.seq),
                'md5': hashlib.md5(str(seq_record.seq)).hexdigest()
            }
            lengths.append(contig['length'])
            contigset_data['contigs'].append(contig)

        # load the method provenance from the context object
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        input_objects = []
        for read in params['reads']:
            read_name = read['read_library_name']
            if '/' in read_name:
                ref = read_name
            else:
                ref = ws_name + '/' + read_name
            input_objects.append(ref)
        provenance[0]['input_ws_objects'] = input_objects

        # save the contigset output
        wsObj = {
            'id': wsid,  # set the output workspace ID
            'objects': [
                {
                    'type': 'KBaseGenomes.ContigSet',
                    'data': contigset_data,
                    'name': params['output_contigset_name'],
                    'meta': {},
                    'provenance': provenance
                }
            ]
        }
        new_obj_info = ws.save_objects(wsObj)
        if new_obj_info is None:
            self.log(console, "Failed to save object")

    def save_assembly(self, wsname, output_contigs, token, name, console):
        self.log(console, 'Uploading FASTA file to Assembly')
        assemblyUtil = AssemblyUtil(self.callbackURL, token=token,
                                    service_ver='dev')
        assemblyUtil.save_assembly_from_fasta({'file': {'path': output_contigs},
                                               'workspace_name': wsname,
                                               'assembly_name': name
                                               })
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.scratch = os.path.abspath(config['scratch'])
        self.callbackURL = os.environ.get('SDK_CALLBACK_URL')
        print "Callback=%s" % (self.callbackURL)
        #END_CONSTRUCTOR
        pass

    def run_hipmer_hpc(self, ctx, params):
        """
        :param params: instance of type "AssemblyParams" (Run assembler
           workspace_name - the name of the workspace for input/output
           read_library_name - the name of the PE read library (SE library
           support in the future) output_contig_set_name - the name of the
           output contigset extra_params - assembler specific parameters
           min_contig_length - minimum length of contigs to output, default
           200 @optional min_contig_len @optional extra_params) -> structure:
           parameter "workspace_name" of String, parameter
           "read_library_name" of String, parameter "output_contigset_name"
           of String, parameter "min_contig_len" of Long, parameter
           "extra_params" of list of String
        :returns: instance of type "AssemblyOutput" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_hipmer_hpc
        console = []
        self.log(console, 'Running run_hipmer_hpc with params=')
        self.log(console, pformat(params))

        # do some basic checks
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'reads' not in params:
            raise ValueError('reads parameter is required')
        if 'output_contigset_name' not in params:
            raise ValueError('output_contigset_name parameter is required')
        ws_name = params['workspace_name']
        #ws = workspaceService(self.workspaceURL, token=ctx['token'])

        if 'POST' not in os.environ:
            # Get the read library
            print "Running pre stage"
            refs = []
            for read in params['reads']:
                read_name = read['read_library_name']
                if '/' in read_name:
                    ref = read_name
                else:
                    ref = ws_name + '/' + read_name
                refs.append(ref)
                read['ref'] = ref

            params['readsfiles'] = self.get_reads_RU(ctx, refs, console)

            # set the output location
            output_dir = self.scratch
            # Generate submit script
            ts = self.generate_config(params)
            self.generate_submit(ts)
            return

        print "Running POST stage"

        # run hipmer, capture output as it happens
        self.log(console, 'running hipmer:')

        output_dir = self.scratch
        output_contigs = os.path.join(output_dir, 'final_assembly.fa')
        output_name = params['output_contigset_name']
        wsname = params['workspace_name']
        #output_data_ref = self.save_assembly(wsname,
        #                                     output_contigs,
        #                                     ctx['token'],
        #                                     output_name,
        #                                     console)
        self.log(console, 'Uploading FASTA file to Assembly')
        assemblyUtil = AssemblyUtil(self.callbackURL, token=ctx['token'],
                                    service_ver='dev')
        save_input = {'file': {'path': output_contigs},
                      'workspace_name': wsname,
                      'assembly_name': output_name
                      }
        output_data_ref = assemblyUtil.save_assembly_from_fasta(save_input)
        print 'ref: ' + output_data_ref
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
        except QUASTError as qe:
            # not really any way to test this, all inputs have been checked
            # earlier and should be ok
            print('Logging exception from running QUAST')
            print(str(qe))
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
        except _RepError as re:
            # not really any way to test this, all inputs have been checked earlier and should be
            # ok
            print('Logging exception from creating report object')
            print(str(re))
            # TODO delete shock node
            raise

        # STEP 6: contruct the output to send back
        output = {'report_name': report_info['name'],
                  'report_ref': report_info['ref']
                  }
        return [output]
        #END run_hipmer_hpc

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_hipmer_hpc return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION,
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]

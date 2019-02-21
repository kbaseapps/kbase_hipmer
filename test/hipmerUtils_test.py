import unittest
import os
import time
import requests
import shutil

from os import environ
from configparser import ConfigParser
from pprint import pprint
from subprocess import call as Call

from hipmer.hipmerUtils import hipmerUtils


class hipmerTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        print(config_file)
        for nameval in config.items('hipmer'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.scratch = cls.cfg['scratch']
        cls.test_dir = os.path.dirname(os.path.abspath(__file__))
        cls.data_dir = os.path.join(cls.test_dir, '../test/data/')

    def test_validate(self):
        # run hipmer
        params = {
            'mer_sizes': '21,41,127',
            'workspace_name': 'bogus',
            'output_contigset_name': 'hipmer.contigs',
            'is_meta': {
                'aggressive': 1,
                'min_depth_cutoff': 7,
            },
            'usedebug': 1,
            'interleaved': 1,
            'reads': [{
                'read_type': 'paired',
                'ins_avg': 100,
                'ins_dev': 10,
                'is_rev_comped': 0,
                'read_library_name': '24799/3/1'
            }]
        }
        HU = hipmerUtils(self.cfg, self.token)
        result = HU._validate_inputs(params)
        self.assertTrue(result)

        # must be ints
        params['mer_sizes'] = 'asdf,asdfa,asdf'
        with self.assertRaises(ValueError):
            result = HU._validate_inputs(params)

        # can't go over 127 for kmers
        params['mer_sizes'] = '10 31 130'
        with self.assertRaises(ValueError):
            result = HU._validate_inputs(params)

        # can't go over 127 for kmers
        params['mer_sizes'] = '10,31,130'
        with self.assertRaises(ValueError):
            result = HU._validate_inputs(params)

        # ins_avg and ins_dev must not be 'None' when paired-end
        params['reads'][0]['ins_dev'] = None
        result = HU.check_reads('24799/3/1', [], params)
        self.assertFalse(result)

#    def test_config(self):
#       configf = os.path.join(self.scratch, 'hipmer.config')
#       if os.path.exists(configf):
#           os.remove(configf)
#       readobj1 = {
#           'files': {'fwd': 'unmerged.q', 'otype': 'paired'},
#           'insert_size_mean': '1000',
#           'insert_size_std_dev': '100',
#           'read_length_mean': '100',
#           'total_bases': 10000
#       }
#       readobj2 = {
#           'files': {'fwd': 'merged.q', 'otype': 'single'},
#           'insert_size_mean': '1000',
#           'insert_size_std_dev': '100',
#           'read_length_mean': '100',
#           'total_bases': 10000
#       }
#       params = {
#           'mer_sizes': '21,41,127',
#           'workspace_name': 'bogus',
#           'output_contigset_name': 'hipmer.contigs',
#           'is_meta': {
#               'aggressive': 1,
#               'min_depth_cutoff': 7,
#           },
#           'usedebug': 1,
#           'interleaved': 1,
#           'reads': [{
#               'read_type': 'paired',
#               'ins_avg': 100,
#               'ins_dev': 10,
#               'is_rev_comped': 0,
#               'read_library_name': 'bogus'
#           }]
#       }
#
#
#        #reads_obj = params['readsfiles'][r['ref']]
#        HU = kbase_hipmerUtils(self.cfg, self.token)
##        result = HU.generate_config(params)
##        self.assertTrue(result)
#        self.assertTrue(os.path.exists(configf))
#        configs = dict()
#        libs = dict()
#        # Read in config file and parse contents
#        with open(configf) as conf:
#            for line in conf:
#                line = line.rstrip().replace('  ', ' ').replace('  ', ' ')
#                print(line)
#                tl = line.split(' ')
#                if len(tl) == 2:
#                    configs[tl[0]] = tl[1]
#                if len(tl) > 2:
#                    lib = tl[1]
#                    libs[lib] = {'ins': tl[3], 'ct': tl[13]}
#        #lib_seq ./517bf2c4-2e3d-4341-a196-12d0a5f20bbf.inter.fastq small 250
#        #    10   100 0 1   1 1 1  0 0 1 1
#
#        # Confirm metagenome options are set as expected
#        for k in ['alpha', 'beta', 'tau', 'error_rate']:
#            self.assertIn(k, configs)
#        self.assertEqual(configs['is_metagenome'], '1')
#        self.assertEqual(libs['unmerged.q']['ins'], '1000')
#        self.assertEqual(libs['unmerged.q']['ct'], '1')
#        self.assertEqual(libs['merged.q']['ct'], '0')

    def test_fixup(self):
        tfile = self.scratch + '/t.fq'
        params = {
            'reads': [{'ref': '1/2/3'}],
            'readsfiles': {
                '1/2/3': {'files': {'fwd': tfile}}
            }
        }
        shutil.copyfile(self.data_dir + 'sra.fq', tfile)
        HU = hipmerUtils(self.cfg, self.token)
        HU.fixup_reads(params)
        self.assertTrue(os.path.exists(tfile+'.orig'))

        os.remove(tfile + '.orig')
        shutil.copyfile(self.data_dir + 'nosra.fq', tfile)
        HU.fixup_reads(params)
        self.assertFalse(os.path.exists(tfile + '.orig'))

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
from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.AbstractHandleClient import AbstractHandle as HandleService

class hipmerTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = 'bogus'
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('hipmer'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.scratch = cls.cfg['scratch']
        cls.test_dir = os.path.dirname(os.path.abspath(__file__))
        cls.data_dir = os.path.join(cls.test_dir, '../test/data/')

    def test_validate(self):
        # run hipmer
        params = {
            'assm_scaff_len_cutoff': 1000,
            'mer_sizes': "21",
            'workspace_name': 'bogus',
            'output_contigset_name': 'hipmer.contigs',
            'min_depth_cutoff': 7,
            'type': 'uniploid',
            'dynamic_min_depth': 1,
            'gap_close_rpt_depth_ratio': 2,
            'reads': [{
                'use_for_splinting': 1,
                'use_for_gap_closing': 1,
                'has_innie_artifact': 0,
                'use_for_contigging': 1,
                'prefix': 'small',
                'ono_set_id': 1,
                'tp_wiggle_room': 0,
                'fp_wiggle_room': 0,
                'read_library_name': 'bogus'
            }]
        }
        HU = hipmerUtils(self.cfg, self.token)
        result = HU._validate_inputs(params)
        self.assertTrue(result)
        params['mer_sizes'] = 'asdf,asdfa,asdf'
        with self.assertRaises(ValueError):
            result = HU._validate_inputs(params)
        params['mer_sizes'] = '1,11,100'
        with self.assertRaises(ValueError):
            result = HU._validate_inputs(params)

    def test_config(self):
        # run hipmer
        configf = os.path.join(self.scratch, 'hipmer.config')
        if os.path.exists(configf):
            os.remove(configf)
        readobj1 = {
            'files': {'fwd': 'unmerged.q', 'otype': 'paired'},
            'insert_size_mean': '1000',
            'insert_size_std_dev': '100',
            'read_length_mean': '100',
            'total_bases': 10000
        }
        readobj2 = {
            'files': {'fwd': 'merged.q', 'otype': 'single'},
            'insert_size_mean': '1000',
            'insert_size_std_dev': '100',
            'read_length_mean': '100',
            'total_bases': 10000
        }
        params = {
            'assm_scaff_len_cutoff': 1000,
            'mer_sizes': "21",
            'workspace_name': 'bogus',
            'output_contigset_name': 'hipmer.contigs',
            'min_depth_cutoff': 7,
            'type': 'metagenome',
            'alpha': 0.1,
            'beta': 0.2,
            'tau': 2.0,
            'bubble_min_depth_cutoff': 1,
            'error_rate': 0.9,
            'dynamic_min_depth': 1,
            'gap_close_rpt_depth_ratio': 2,
            'mer_max': 250,
            'la_mer_max': 252,
            'mer_sizes_int': [250],
            'readsfiles': {'1/2/3': readobj1, '1/2/4': readobj2},
            'reads': [{
                'use_for_splinting': 1,
                'use_for_gap_closing': 1,
                'has_innie_artifact': 0,
                'use_for_contigging': 1,
                'prefix': 'small',
                'ono_set_id': 1,
                'tp_wiggle_room': 0,
                'fp_wiggle_room': 0,
                'read_library_name': 'bogus',
                'ref': '1/2/3'
            }, {
                'use_for_splinting': 1,
                'use_for_gap_closing': 1,
                'has_innie_artifact': 0,
                'use_for_contigging': 1,
                'prefix': 'small',
                'ono_set_id': 1,
                'tp_wiggle_room': 0,
                'fp_wiggle_room': 0,
                'read_library_name': 'bogus2',
                'ref': '1/2/4'
            }]
        }

        #reads_obj = params['readsfiles'][r['ref']]
        HU = hipmerUtils(self.cfg, self.token)
        result = HU.generate_config(params)
        self.assertTrue(result)
        self.assertTrue(os.path.exists(configf))
        configs = dict()
        libs = dict()
        # Read in config file and parse contents
        with open(configf) as conf:
            for line in conf:
                line = line.rstrip().replace('  ', ' ').replace('  ', ' ')
                print(line)
                tl = line.split(' ')
                if len(tl) == 2:
                    configs[tl[0]] = tl[1]
                if len(tl) > 2:
                    lib = tl[1]
                    libs[lib] = {'ins': tl[3], 'ct': tl[13]}
        #lib_seq ./517bf2c4-2e3d-4341-a196-12d0a5f20bbf.inter.fastq small 250
        #    10   100 0 1   1 1 1  0 0 1 1

        # Confirm metagenome options are set as expected
        for k in ['alpha', 'beta', 'tau', 'error_rate']:
            self.assertIn(k, configs)
        self.assertEqual(configs['is_metagenome'], '1')
        self.assertEqual(libs['unmerged.q']['ins'], '1000')
        self.assertEqual(libs['unmerged.q']['ct'], '1')
        self.assertEqual(libs['merged.q']['ct'], '0')

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

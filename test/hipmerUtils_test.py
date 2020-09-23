import unittest
import os
import time
import requests
import shutil
import json
from copy import deepcopy

from os import environ
from configparser import ConfigParser
from pprint import pprint
from subprocess import call as Call
from unittest.mock import MagicMock

from hipmer.hipmerUtils import hipmerUtils

MOCK_GET_READS = {'files': 
    {
    '52407/2/1': {'base_percentages': {'A': 16.0606,
                                    'C': 34.1326,
                                    'G': 33.952,
                                    'N': 0.0,
                                    'T': 15.8549},
               'files': {'fwd': 'small.inter.fq',
                         'fwd_name': 'small.inter.fq.gz',
                         'otype': 'interleaved',
                         'rev': None,
                         'rev_name': None,
                         'type': 'interleaved'
                         },
               'gc_content': 0.680846,
               'insert_size_mean': None,
               'insert_size_std_dev': None,
               'number_of_duplicates': 11,
               'phred_type': '33',
               'qual_max': 51.0,
               'qual_mean': 43.0606,
               'qual_min': 10.0,
               'qual_stdev': 10.5302,
               'read_count': 2500,
               'read_length_mean': 100.0,
               'read_length_stdev': 0.0,
               'read_orientation_outward': 'false',
               'read_size': None,
               'ref': '52407/2/1',
               'sequencing_tech': 'Illumina',
               'single_genome': 'true',
               'source': None,
               'strain': None,
               'total_bases': 250000
               }
    }
}

MOCK_GET_INFO = {
    'Insert_Size_Std_Dev': 'Not Specified',
    'Phred_Type': '33',
    'Number_of_Duplicate_Reads': '11 (0.44%)',
    'Name': 'small.interlaced_reads',
    'Strain': 'Not Specified',
    'GC_Percentage': '68.08%',
    'Type': 'Paired End',
    'Read_Length_Std_Dev': '0.0',
    'Insert_Size_Mean': 'Not Specified',
    'Quality_Score_Mean_Std_Dev': '43.06 (10.53)',
    'Mean_Read_Length': '100.0',
    'Quality_Score_Min_Max': '10.0/51.0',
    'Total_Number_of_Bases': '250,000',
    'Source': 'Not Specified',
    'Base_Percentages': 'A(16.06%), C(34.13%), G(33.95%), T(15.85%), N(0.0%)',
    'Single_Genome': 'Yes',
    'Platform': 'Illumina',
    'workspace_name': 'KBaseTestData',
    'Number_of_Reads': '2,500',
    'id': 2,
    'Outward_Read_Orientation': 'No'
}
#    'is_meta': {
#        'aggressive': 1
#    },
#    "is_plant": {
#        "diploid": "low"
#    },

PARAMS = {
    'mer_sizes': '21,41,127',
    'workspace_name': 'bogus',
    'output_contigset_name': 'hipmer.contigs',
    'usedebug': 1,
    'reads': [{
        'ins_avg': 100,
        'ins_dev': 10,
        'is_rev_comped': 0,
        'read_library_name': '52407/2/1'
    }]
}


class hipmerTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        print(config_file)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        print(config_file)
        for nameval in config.items('hipmer'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.scratch = cls.cfg['scratch']
        cls.test_dir = os.path.dirname(os.path.abspath(__file__))
        cls.data_dir = os.path.join(cls.test_dir, '../test/data/')

    def _get_cmd(self):
        with open(self.scratch +'/slurm2.sl') as f:
            for line in f:
                if line.startswith('$'):
                    return line

    def test_validate(self):
        params = deepcopy(PARAMS)
        params['is_meta'] = 1
        params['agressive'] = 1
        # run hipmer
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


    def test_reads(self):
        hu = hipmerUtils(self.cfg, self.token)
        params = deepcopy(PARAMS)
        params['is_meta'] = 1
        params['agressive'] = 1
        mock_get_reads = deepcopy(MOCK_GET_READS)
        fn = os.path.join(self.data_dir, 'small.inter.fq')
        mock_get_reads['files']['52407/2/1']['files']['fwd'] = fn

        hu.rapi.get_reads_info_all_formatted = MagicMock(return_value = MOCK_GET_INFO)
        hu.readcli.download_reads = MagicMock(return_value = mock_get_reads)

        # Test with specified settings
        hu.prepare_run(params)
        line = self._get_cmd()
        self.assertIn(':i100:s10', line)

        # Test with no params and nothing in reads objects
        params['reads'][0]['ins_avg'] = None
        params['reads'][0]['ins_dev'] = None
        hu.prepare_run(params)
        line = self._get_cmd()
        self.assertNotIn(':i', line)

        # Test with no params but reads have it
        mock_get_reads['files']['52407/2/1']['insert_size_mean'] = 123
        mock_get_reads['files']['52407/2/1']['insert_size_std_dev'] = 5
        hu.prepare_run(params)
        line = self._get_cmd()
        self.assertIn(':i123:s5', line)

        # Test with no params but reads have it
        params['reads'][0]['ins_avg'] = 234
        params['reads'][0]['ins_dev'] = 6
        hu.prepare_run(params)
        line = self._get_cmd()
        self.assertIn(':i234:s6', line)

        mock_get_reads['files']['52407/2/1']['read_orientation_outward'] = 'true'
        # Test with no params but reads have it
        hu.prepare_run(params)
        line = self._get_cmd()
        self.assertIn(':rc', line)



##    def test_generate_command(self):
##        HU = kbase_hipmerUtils(self.cfg, self.token)
##        result = HU.generate_config(params)
##        self.assertTrue(result)

#    def test_fixup(self):
#        tfile = self.scratch + '/t.fq'
#        params = {
#            'reads': [{'ref': '1/2/3'}],
#            'readsfiles': {
#                '1/2/3': {'files': {'fwd': tfile}}
#            }
#        }
#        shutil.copyfile(self.data_dir + 'sra.fq', tfile)
#        HU = hipmerUtils(self.cfg, self.token)
#        HU.fixup_reads(params)
#        self.assertTrue(os.path.exists(tfile+'.orig'))
#
#        os.remove(tfile + '.orig')
#        shutil.copyfile(self.data_dir + 'nosra.fq', tfile)
#        HU.fixup_reads(params)
#        self.assertFalse(os.path.exists(tfile + '.orig'))
#

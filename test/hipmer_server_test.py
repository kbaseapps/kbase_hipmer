import unittest
import os
import time
import requests

from os import environ
from ConfigParser import ConfigParser
from pprint import pprint
from requests_toolbelt import MultipartEncoder
from subprocess import call as Call

from biokbase.workspace.client import Workspace as workspaceService
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService
from hipmer.hipmerImpl import hipmer
from hipmer.hipmerServer import MethodContext


class hipmerTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'provenance': [
                            {'service': 'hipmer',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('hipmer'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.ws = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = hipmer(cls.cfg)
        cls.shockURL = cls.cfg['shock-url']
        cls.handleURL = cls.cfg['handle-service-url']
        cls.scratch = cls.cfg['scratch']
        print "shock %s" % (cls.shockURL)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.ws.delete_workspace({'workspace': cls.wsName})
            print cls.wsName
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.ws

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_hipmer_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        print ret
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    # call this method to get the WS object info of a Paired End Library (will
    # upload the example data if this is the first time the method is called during tests)
    def getPairedEndLibInfo(self):
        if hasattr(self.__class__, 'pairedEndLibInfo'):
            return self.__class__.pairedEndLibInfo
        # 1) upload files to shock
        token = self.ctx['token']
        forward_shock_file = self.upload_file_to_shock(
            shock_service_url=self.shockURL,
            filePath='data/small.forward.fq',
            token=token)
        reverse_shock_file = self.upload_file_to_shock(
            shock_service_url=self.shockURL,
            filePath='data/small.reverse.fq',
            token=token)
        #pprint(forward_shock_file)
        #pprint(reverse_shock_file)

        # 2) create handle
        hs = HandleService(url=self.handleURL, token=token)
        handf = {
            'id': forward_shock_file['id'],
            'type': 'shock',
            'url': self.shockURL,
            'file_name': forward_shock_file['file']['name'],
            'remote_md5': forward_shock_file['file']['checksum']['md5']
        }
        forward_handle = hs.persist_handle(handf)
        handr = {
            'id': reverse_shock_file['id'],
            'type': 'shock',
            'url': self.shockURL,
            'file_name': reverse_shock_file['file']['name'],
            'remote_md5': reverse_shock_file['file']['checksum']['md5']
        }

        reverse_handle = hs.persist_handle(handr)

        # 3) save to WS
        paired_end_library = {
            'lib1': {
                'file': {
                    'hid': forward_handle,
                    'file_name': forward_shock_file['file']['name'],
                    'id': forward_shock_file['id'],
                    'url': self.shockURL,
                    'type': 'shock',
                    'remote_md5': forward_shock_file['file']['checksum']['md5']
                },
                'encoding': 'UTF8',
                'type': 'fastq',
                'size': forward_shock_file['file']['size']
            },
            'lib2': {
                'file': {
                    'hid': reverse_handle,
                    'file_name': reverse_shock_file['file']['name'],
                    'id': reverse_shock_file['id'],
                    'url': self.shockURL,
                    'type': 'shock',
                    'remote_md5': reverse_shock_file['file']['checksum']['md5']
                },
                'encoding': 'UTF8',
                'type': 'fastq',
                'size': reverse_shock_file['file']['size']

            },
            'interleaved': 0,
            'sequencing_tech': 'artificial reads',
            'read_length_mean': 100,
            'insert_size_mean': 250,
            'insert_size_std_dev': 10,
            'total_bases': 125000,
            'read_orientation_outward': 1
        }
        ws_obj = {
            'workspace': self.getWsName(),
            'objects': [
                {
                    'type': 'KBaseFile.PairedEndLibrary',
                    'data': paired_end_library,
                    'name': 'test.pe.reads',
                    'meta': {},
                    'provenance': [
                        {
                            'service': 'hipmer',
                            'method': 'test_hipmer'
                        }
                    ]
                }]
        }

        new_obj_info = self.ws.save_objects(ws_obj)
        self.__class__.pairedEndLibInfo = new_obj_info[0]
        return new_obj_info[0]

    def createBogus(self, fname):
        outdir = self.scratch + '/results/'
        if os.path.exists(outdir) is False:
            os.makedirs(outdir)
        print 'dest: %s/%s' % (outdir, fname)
        ret = Call(['cp', 'data/output.contig.fa', '%s/%s' % (outdir, fname)])
        print ret

    def test_0validate(self):
        # run hipmer
        params = {
            'assm_scaff_len_cutoff': 1000,
            'mer_sizes': "21",
            'workspace_name': 'bogus',
            'output_contigset_name': 'hipmer.contigs',
            'min_depth_cutoff': 7,
            'is_diploid': None,
            'is_metagenome': None,
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
        result = self.getImpl()._validate_inputs(params)
        self.assertTrue(result)
        params['mer_sizes'] = 'asdf,asdfa,asdf'
        with self.assertRaises(ValueError):
            result = self.getImpl()._validate_inputs(params)
        params['mer_sizes'] = '1,11,100'
        with self.assertRaises(ValueError):
            result = self.getImpl()._validate_inputs(params)
        # Simulate diploid and metagenome both being set
        params['mer_sizes'] = '21'
        params['is_diploid'] = {'bubble_min_depth_cutoff': 1}
        params['is_metagenome'] = {'alpha': 0.1, 'beta': 0.2, 'tau': 2.0}
        with self.assertRaises(ValueError):
            result = self.getImpl()._validate_inputs(params)

    def test_0config(self):
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
            'is_diploid': None,
            'is_metagenome': {'alpha': 0.1, 'beta': 0.2, 'tau': 2.0},
            'dynamic_min_depth': 1,
            'gap_close_rpt_depth_ratio': 2,
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
            },{
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
        result = self.getImpl().generate_config(params)
        self.assertTrue(result)
        self.assertTrue(os.path.exists(configf))
        configs = dict()
        # Read in config file and parse contents
        with open(configf) as conf:
            for lines in conf:
                print lines.rstrip()
                templist = lines.rstrip().split(' ')
                if len(templist) == 2:
                    configs[templist[0]] = templist[1]
        # Confirm metagenome options are set as expected
        self.assertIn('alpha', configs)
        self.assertEquals(configs['is_metagenome'], '1')

    def test_post(self):
        pe_lib_info = self.getPairedEndLibInfo()
        #pprint(pe_lib_info)

        # run hipmer
        params = {
            'assm_scaff_len_cutoff': 1000,
            'mer_sizes': '21',
            'workspace_name': pe_lib_info[7],
            'output_contigset_name': 'hipmer.contigs',
            'min_depth_cutoff': 7,
            'is_diploid': None,
            'is_metagenome': None,
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
                'read_library_name': pe_lib_info[1]
            }]
        }
        self.createBogus('final_assembly.fa')
        os.environ['POST'] = '1'

        result = self.getImpl().run_hipmer_hpc(self.getContext(), params)

        print('RESULT:')
        pprint(result)

        # check the output
        info_list = self.ws.get_object_info([{'ref': pe_lib_info[7] + '/hipmer.contigs'}], 1)
        self.assertEqual(len(info_list), 1)
        contigset_info = info_list[0]
        self.assertEqual(contigset_info[1], 'hipmer.contigs')
        self.assertEqual(contigset_info[2].split('-')[0], 'KBaseGenomeAnnotations.Assembly')
        #self.assertEqual(contigset_info[2].split('-')[0], 'KBaseGenomes.ContigSet')

    def test_run_hipmer(self):

        # figure out where the test data lives
        pe_lib_info = self.getPairedEndLibInfo()
        #pprint(pe_lib_info)
        if 'POST' in os.environ:
            del os.environ['POST']

        # run hipmer
        params = {
            'assm_scaff_len_cutoff': 1000,
            'mer_sizes': '21',
            'workspace_name': pe_lib_info[7],
            'output_contigset_name': 'hipmer.contigs',
            'min_depth_cutoff': 7,
            'is_diploid': None,
            'is_metagenome': None,
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
                'read_library_name': pe_lib_info[1]
            }]
        }

        result = self.getImpl().run_hipmer_hpc(self.getContext(), params)
        self.createBogus('final_assembly.fa')

        print('RESULT:')
        pprint(result)
        self.assertEqual(os.path.exists(self.scratch + '/hipmer.config'), True)
        self.assertEqual(os.path.exists(self.scratch + '/slurm.submit'), True)
        # check the output
        ##info_list = self.ws.get_object_info([{'ref':pe_lib_info[7] + '/output.contigset'}], 1)
        #self.assertEqual(len(info_list),1)
        ##contigset_info = info_list[0]
        #self.assertEqual(contigset_info[1],'output.contigset')
        #self.assertEqual(contigset_info[2].split('-')[0],'KBaseGenomes.ContigSet')

    # Helper script borrowed from the transform service, logger removed
    def upload_file_to_shock(self,
                             shock_service_url=None,
                             filePath=None,
                             ssl_verify=True,
                             token=None):
        """
        Use HTTP multi-part POST to save a file to a SHOCK instance.
        """

        if token is None:
            raise Exception("Authentication token required!")

        #build the header
        header = dict()
        header["Authorization"] = "Oauth {0}".format(token)

        if filePath is None:
            raise Exception("No file given for upload to SHOCK!")

        dataFile = open(os.path.abspath(filePath), 'rb')
        m = MultipartEncoder(fields={'upload': (os.path.split(filePath)[-1], dataFile)})
        header['Content-Type'] = m.content_type

        #logger.info("Sending {0} to {1}".format(filePath,shock_service_url))
        try:
            response = requests.post(shock_service_url + "/node",
                                     headers=header, data=m,
                                     allow_redirects=True,
                                     verify=ssl_verify)
            dataFile.close()
        except:
            dataFile.close()
            raise

        if not response.ok:
            response.raise_for_status()

        result = response.json()

        if result['error']:
            raise Exception(result['error'][0])
        else:
            return result["data"]

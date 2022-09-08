import unittest
import os
import time
import requests
import shutil

from os import environ
from configparser import ConfigParser
from pprint import pprint
from requests_toolbelt import MultipartEncoder
from subprocess import call as Call

from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.AbstractHandleClient import AbstractHandle as HandleService
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
        print("shock %s" % (cls.shockURL))
        cls.testWS = 'KBaseTestData'
        cls.testReads = 'small.interlaced_reads'

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.ws.delete_workspace({'workspace': cls.wsName})
            print(cls.wsName)
            print('Test workspace was deleted')

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_hipmer_" + str(suffix)
        ret = self.ws.create_workspace({'workspace': wsName})
        print(ret)
        self.__class__.wsName = wsName
        return wsName

    def createBogus(self, fname):
        outdir = self.scratch + '/results/'
        if os.path.exists(outdir) is False:
            os.makedirs(outdir)
        print('dest: %s/%s' % (outdir, fname))
        ret = Call(['cp', 'data/output.contig.fa', '%s/%s' % (outdir, fname)])
        print(ret)

    def xtest_post(self):

        # run hipmer
        wsn = self.getWsName()
        params = {
            'mer_sizes': '21,41,127',
            'scaff_mer_sizes': '99,33',
            'workspace_name': wsn,
            'output_contigset_name': 'hipmer.contigs',
            'is_meta': 1,
            'aggressive': 1,
            'usedebug': 0,
            'reads': [{
                'ins_avg': 100,
                'ins_dev': 10,
                'is_rev_comped': 0,
                'read_library_name': 'bogus'
            }]
        }

        self.createBogus('final_assembly.fa')
        os.environ['POST'] = '1'

        result = self.serviceImpl.run_hipmer_hpc(self.ctx, params)

        print('RESULT for POST:')
        pprint(result)

        # check the output
        info_list = self.ws.get_object_info([{'ref': wsn + '/hipmer.contigs'}], 1)
        self.assertEqual(len(info_list), 1)
        contigset_info = info_list[0]
        self.assertEqual(contigset_info[1], 'hipmer.contigs')
        self.assertEqual(contigset_info[2].split('-')[0], 'KBaseGenomeAnnotations.Assembly')

    # TODO: need to mock up special so we can trap the slurm submit
    def xtest_run_hipmer(self):

        # figure out where the test data lives
        wsn = self.getWsName()
        ret = self.ws.copy_object({'from': {'workspace': self.testWS, 'name': self.testReads},
                        'to': {'workspace': wsn, 'name': self.testReads}})
        upa = '{}/{}/{}'.format(ret[6], ret[0], ret[4])

        #pe_lib_info = self.getPairedEndLibInfo()
        pe_lib_info = self.ws.get_object_info([{'ref': upa}], 1)[0]
        pprint(pe_lib_info)
        if 'POST' in os.environ:
            del os.environ['POST']

        # run hipmer
        params = {
            'mer_sizes': '21,41,127',
            'scaff_mer_sizes': '99,33',
            'workspace_name': pe_lib_info[7],
            'output_contigset_name': 'hipmer.contigs',
            'is_meta': 1,
            'aggressive': 1,
            'assembly_size_filter': 1200,
            'usedebug': 1,
            'interleaved': 1,
            'reads': [{
                'ins_avg': 100,
                'ins_dev': 10,
                'is_rev_comped': 0,
                'read_library_name': pe_lib_info[1]
            }]
        }

        result = self.serviceImpl.run_hipmer_hpc(self.ctx, params)
        self.createBogus('final_assembly.fa')

        print('RESULT:')
        pprint(result)
        self.assertEqual(os.path.exists(self.scratch + '/slurm.submit'), True)

#        # check the output
#        ##info_list = self.ws.get_object_info([{'ref':pe_lib_info[7] + '/output.contigset'}], 1)
#        #self.assertEqual(len(info_list),1)
#        ##contigset_info = info_list[0]
#        #self.assertEqual(contigset_info[1],'output.contigset')
#        #self.assertEqual(contigset_info[2].split('-')[0],'KBaseGenomes.ContigSet')

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

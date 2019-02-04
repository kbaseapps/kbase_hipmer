# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import requests

from hipmer.hipmerUtils import hipmerUtils
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

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.4"
    GIT_URL = "git@github.com:kbaseapps/kbase_hipmer.git"
    GIT_COMMIT_HASH = "31a2f650af82b4caa11491f224f2c0c05544660f"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.scratch = os.path.abspath(config['scratch'])
        self.callbackURL = os.environ.get('SDK_CALLBACK_URL')
        self.config = config
        print("Callback=%s" % (self.callbackURL))
        #END_CONSTRUCTOR

    def run_hipmer_hpc(self, ctx, params):
        """
        :param params: instance of unspecified object
        :returns: instance of type "AssemblyOutput" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_hipmer_hpc
        HU = hipmerUtils(self.config, ctx['token'])
        if 'POST' not in os.environ:
            # Pre-Stage
            HU.prepare_run(params)
            return
        else:
            output = HU.finish_run(params)
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

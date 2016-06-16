'''
A KBase client that dynamically configures itself on request to call
local KBase SDK methods.
This client is currently in the proof-of-concept stage.
This client is not concurrency-safe. Do not rely on environment variables
for initializing tokens or callback URLs if running in a threaded or concurrent
programming based environment.
'''

import json as _json
import requests as _requests
import urlparse as _urlparse
import random as _random
import os as _os
import functools

_CT = 'content-type'
_AJ = 'application/json'
_URL_SCHEME = frozenset(['http', 'https'])
_CB_URL = 'SDK_CALLBACK_URL'
_AUTH = 'AUTHORIZATION'
_TOKEN = 'token'
_KB_TOKEN = 'KB_AUTH_TOKEN'

# TODO add async calls (needs callback service support)
# TODO support calls to the NJSW for narrative methods
# TODO add docs from compilation report when available
# TODO see if there's a way better than partial to generate the functions.
# TODO better docs
# TODO don't blow away the call context.


class ServerError(Exception):
    '''
    The call returned an error. Fields:
    name - the name of the error.
    code - the error code.
    message - a human readable error message.
    data - the server side stacktrace.
    '''

    def __init__(self, name, code, message, data=None, error=None):
        self.name = name
        self.code = code
        self.message = message if message else ''
        self.data = data or error or ''
        # data = JSON RPC 2.0, error = 1.1

    def __str__(self):
        return self.name + ': ' + str(self.code) + '. ' + self.message + \
            '\n' + self.data


class _JSONObjectEncoder(_json.JSONEncoder):

    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        if isinstance(obj, frozenset):
            return list(obj)
        return _json.JSONEncoder.default(self, obj)


class KBDynClient(object):
    '''
    The KBase dynamic client.
    Required initialization arguments (positional):
    catalog_url - the url of the KBase catalog service.
    Optional arguments (keywords in positional order):
    call_context - the method call context (TODO better description)
    callback_url - the url of the Narrative Job Runner callback server. If
        not provided, attempts to retrieve from the call_context and the
        environment in that order.
    token - the KBase authorization token for the user. If not provided,
        attempts to retrieve from the call_context and the
        environment in that order. If no token is found, any method calls that
        require authentication will fail.
    timeout - methods will fail if they take longer than this value in seconds.
        Default 1800.
    trust_all_ssl_certificates - set to True to trust self-signed certificates.
        If you don't understand the implications, leave as the default, False.
    '''

    class _Mods(object):
        pass

    def __init__(self, catalog_url, call_context=None, callback_url=None,
                 token=None, timeout=30 * 60,
                 trust_all_ssl_certificates=False):
        self.catalog_url = catalog_url
        self._check_url(self.catalog_url)
        self.callback_url = self._get_callback_url(callback_url, call_context)
        self.timeout = int(timeout)
        self._headers = dict()
        self.trust_all_ssl_certificates = trust_all_ssl_certificates
        # token overrides user_id and password
        if token:
            self._headers[_AUTH] = token
        elif call_context and call_context.get(_TOKEN):
            self._headers[_AUTH] = call_context.get(_TOKEN)
        elif _os.environ.get(_KB_TOKEN):
            self._headers[_AUTH] = _os.environ.get(_KB_TOKEN)
        if self.timeout < 1:
            raise ValueError('Timeout value must be at least 1 second')
        self.mods = self._Mods()

    def _get_callback_url(self, callback_url, call_context):
        if callback_url:
            url = callback_url
        elif call_context and call_context.get(_CB_URL):
            url = call_context.get(_CB_URL)
        else:
            url = _os.environ.get(_CB_URL)
        if not url:
            raise ValueError('Could not get valid callback url from ' +
                             'arguments, call context, or environment')
        self._check_url(url)
        return url

    def _check_url(self, url):
        scheme, _, _, _, _, _ = _urlparse.urlparse(url)
        if scheme not in _URL_SCHEME:
            raise ValueError(url + " isn't a valid http url")

    def _call(self, url, method, params, call_context=None):
        if call_context and type(call_context) is not dict:
            raise ValueError('call_context must be a dict')
        arg_hash = {'method': method,
                    'params': params,
                    'version': '1.1',
                    'id': str(_random.random())[2:]
                    }
        if call_context:
            arg_hash['context'] = call_context

        body = _json.dumps(arg_hash, cls=_JSONObjectEncoder)
        response = _requests.post(url, data=body, headers=self._headers,
                                  timeout=self.timeout,
                                  verify=not self.trust_all_ssl_certificates)
        response.encoding = 'utf-8'
        if response.status_code == 500:
            if _CT in response.headers and response.headers[_CT] == _AJ:
                err = response.json()
                if 'error' in err:
                    raise ServerError(**err['error'])
                else:
                    raise ServerError('Unknown', 0, response.text)
            else:
                raise ServerError('Unknown', 0, response.text)
        if not response.ok:
            response.raise_for_status()
        resp = response.json()
        if 'result' not in resp:
            raise ServerError('Unknown', 0, 'An unknown server error occurred')
        return resp['result']

    def _prep_call(self, method, version, *args, **keywords):
        del keywords
        ret = self._call(self.callback_url, method, args,
                         call_context={'service_ver': version})
        if len(ret) == 1:
            return ret[0]
        return ret

    def load_module(self, module_name, version=None):
        '''
        Load a module from the KBase Catalog. Module functions will be loaded
        into the client namespace as:
        client.mods.[module name].[method name]
        and can be called as:
        result = client.mods.[module name].[method name](param1, param2, ...,
                call_context={..})

        Required parameters (positional):
        module_name - the name of the module to load.
        Optional parameters:
        version - the version of the module, either 'dev', 'beta', or
           'release' (the default).
        '''

        if version not in ['dev', 'beta']:
            version = 'release'
        ret = self._call(self.catalog_url, 'Catalog.get_module_info',
                         [{'module_name': module_name}])[0]
        name = ret['module_name']
        modinfo = ret[version]
        methods = modinfo['local_functions']
        if not methods:
            raise ValueError('No local functions in module ' + name)
        methstore = self._Mods()
        for m in methods:
            meth = functools.partial(self._prep_call, name + '.' + m, version)
            setattr(methstore, m, meth)
        setattr(self.mods, name, methstore)

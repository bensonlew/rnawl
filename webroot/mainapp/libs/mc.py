import memcache
from web.session import Store


class MemcacheStore(Store):
    def __init__(self, config):
        '''
        config = {
        'servers': ['127.0.0.1:11211'],
        'timeout': 1440
        }
        '''
        self.mc = memcache.Client(config['servers'])
        self.timeout = config['timeout']

    def __contains__(self, key):
        return True if self.mc.get(key) else False

    def __getitem__(self, key):
        return self.mc.get(key)

    def __setitem__(self, key, value):
        self.mc.set(key, value, self.timeout)

    def __delitem__(self, key):
        self.mc.delete(key)

    def cleanup(self, timeout):
        '''You need nothing to do. Memcache can handle it.'''
        pass

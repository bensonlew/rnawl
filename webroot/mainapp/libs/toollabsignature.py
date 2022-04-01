# -*- coding: utf-8 -*-
# __author__ = 'HD'
import hashlib
import web
import functools
import random
import time
import json
# from biocluster.config import Config
from mainapp.config.db import Config


class Signature(object):

    def __init__(self):
        data = web.input()
        info = {
            "signature": data.signature if hasattr(data, "signature") else web.ctx.env.get('HTTP_SIGNATURE'),
            "timestamp": data.timestamp if hasattr(data, "timestamp") else web.ctx.env.get('HTTP_TIMESTAMP'),
            "nonce": data.nonce if hasattr(data, "nonce") else web.ctx.env.get('HTTP_NONCE'),
            # "key": Config().rcf.get("toollab", "authkey")
            "key": Config().get_webauth("tool_lab","authkey")
        }

        for data in info.values():
            if not data:
                raise web.badrequest

        self._signature = info["signature"]
        self._timestamp = info["timestamp"]
        self._nonce = info["nonce"]
        self._key = info['key']

    def check(self):
        x_list = [self._key, self._timestamp, self._nonce]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        hashcode = sha1.hexdigest()
        if hashcode == self._signature:
            return True
        else:
            return False


def check_sig(f):
    @functools.wraps(f)
    def wrapper(obj):
        session = web.config.get('_session')
        if session and session.is_login:
            return f(obj)
        sig = Signature()
        if sig.check():
            return f(obj)
        else:
            info = {"success": False, "info": "身份验证未通过!"}
            return json.dumps(info)
    return wrapper


class CreateSignature(object):
    def __init__(self, client):
        self.client = client
        self.nonce = str(random.randint(1000, 10000))
        self.timestamp = str(int(time.time()))
        # self.key = Config().rcf.get("toollab", "authkey")
        self.key = Config().get_webauth("tool_lab","authkey")
        self.signature = self.get_signature()

    def get_signature(self):
        x_list = [self.key, self.timestamp, self.nonce]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        return sha1.hexdigest()

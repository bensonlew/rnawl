# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
# modified by hd @20200821 增加创新平台的身份验证
import hashlib
import web
from mainapp.models.clientkey import ClientKey
import datetime
import functools
from netaddr import IPNetwork, IPAddress
import random
import time
import json


class Signature(object):

    def __init__(self):
        data = web.input()
        info = {
            "signature": data.signature if hasattr(data, "signature") else web.ctx.env.get('HTTP_SIGNATURE'),
            "timestamp": data.timestamp if hasattr(data, "timestamp") else web.ctx.env.get('HTTP_TIMESTAMP'),
            "nonce": data.nonce if hasattr(data, "nonce") else web.ctx.env.get('HTTP_NONCE'),
            "client": data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT'),
        }

        for data in info.values():
            if not data:
                raise web.badrequest

        self._signature = info["signature"]
        self._timestamp = info["timestamp"]
        self._nonce = info["nonce"]
        self._client = info["client"]
        self._ip = web.ctx.ip

    def check(self):
        # diff = datetime.datetime.now() - datetime.datetime.fromtimestamp(int(self._timestamp))
        client = ClientKey(self._client)
        # if client.timelimit and abs(diff.seconds) > client.timelimit:
        #     info = {"success": False, "info": "验证时间超时!"}
        #     return json.dumps(info)
        # if client.ipmask:
        #     ip_list = client.ipmask.split(";")
        #     in_range = False
        #     for r in ip_list:
        #         if IPAddress(self._ip) in IPNetwork(r):
        #             in_range = True
        #     if in_range is False:
        #         # raise web.unauthorized
        #         info = {"success": False, "info": "IP不在允许范围内!"}
        #         return json.dumps(info)

        x_list = [client.key, self._timestamp, self._nonce]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        hashcode = sha1.hexdigest()
        if hashcode == self._signature:
            return True
        else:
            return False

    def checkV2(self):
        """
        创新平台的接口验证
        :return:
        """
        client = ClientKey(self._client)
        x_list = [client.n_key, self._timestamp, self._nonce]
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
            # raise web.unauthorized
            info = {"success": False, "info": "身份验证未通过!"}
            return json.dumps(info)
    return wrapper


def check_sig_v2(f):
    @functools.wraps(f)
    def wrapper(obj):
        session = web.config.get('_session')
        if session and session.is_login:
            return f(obj)
        sig = Signature()
        if sig.checkV2():
            return f(obj)
        else:
            # raise web.unauthorized
            info = {"success": False, "info": "身份验证未通过!"}
            return json.dumps(info)
    return wrapper


class CreateSignature(object):
    def __init__(self, client):
        self.client = client
        self.nonce = str(random.randint(1000, 10000))
        self.timestamp = str(int(time.time()))
        self.signature = self.get_signature()

    def get_signature(self):
        data = ClientKey(self.client)
        x_list = [data.key, self.timestamp, self.nonce]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        return sha1.hexdigest()


class CreateSignatureV2(object):
    def __init__(self, client):
        self.client = client
        self.nonce = str(random.randint(1000, 10000))
        self.timestamp = str(int(time.time()))
        self.signature = self.get_signature()

    def get_signature(self):
        data = ClientKey(self.client)
        x_list = [data.n_key, self.timestamp, self.nonce]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        return sha1.hexdigest()

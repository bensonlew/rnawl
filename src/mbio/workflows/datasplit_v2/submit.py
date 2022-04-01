# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190311

"""
接口投递，用于将拆分的文库拆分、样本拆分、样本质控独立开
"""

import re
import os
import time
import json
import random
import urllib
import urllib2
import hashlib
import httplib
import datetime
import argparse
from socket import error as SocketError
# from mainapp.libs.signature import CreateSignature
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError

class Submit(object):
    """投递接口的对象"""
    def __init__(self, params, url, types=None):
        if types:
            if types == "sanger":
                # self.url = "http://bcl.i-sanger.com"
                self.url = "http://wpm.i-sanger.com"
                self.client = "client01"
                self.mysql_client_key = "1ZYw71APsQ"
            elif types == "tsanger":
                # self.url = "http://bcl.tsanger.com"
                self.url = "http://wpm.i-sanger.com"
                # self.client = "client03"
                self.client = "client01"
                # self.mysql_client_key = "hM4uZcGs9d"
                self.mysql_client_key = "1ZYw71APsQ"
            elif types == "nbtsg":
                self.url = "http://wpm2.sanger.com"
                self.client = "client03"
                self.mysql_client_key = "hM4uZcGs9d"
            else:
                # self.url = "http://10.101.203.193:9090"
                # self.url = "http://192.168.12.101:9090"
                self.url = "http://bcl.tsg.com"
                self.client = "client03"
                self.mysql_client_key = 'hM4uZcGs9d'
        else:
            # self.url = "http://bcl.i-sanger.com"
            self.url = "http://wpm.i-sanger.com"
            self.client = "client01"
            self.mysql_client_key = "1ZYw71APsQ"
        self.api = url
        self.params = params

    def webapitest(self):
        """用于投递接口"""
        httpHandler = urllib2.HTTPHandler(debuglevel=1)
        httpsHandler = urllib2.HTTPSHandler(debuglevel=1)
        opener = urllib2.build_opener(httpHandler, httpsHandler)
        urllib2.install_opener(opener)
        dataa = urllib.urlencode(self.params)
        url = "%s?%s" % (self.url + self.api, self.signature())
        request = urllib2.Request(url, dataa)
        for i in range(3):
            try:
                response = urllib2.urlopen(request)
            except (urllib2.HTTPError, urllib2.URLError, httplib.HTTPException, SocketError) as e:
                print "接口投递失败, 第{}次(最多三次尝试), 错误:{}".format(i + 1, e)
            else:
                the_page = response.read()
                print "接口投递成功"
                return json.loads(the_page)
        # return {"info": "拆分接口投递失败！", "success": False}

    def signature(self):
        timestamp = str(int(time.time()))
        nonce = str(random.randint(1000, 10000))
        web_key = self.mysql_client_key
        sha1 = hashlib.sha1()
        key_list = [web_key, nonce, timestamp]
        key_list.sort()
        map(sha1.update, key_list)
        hashkey = sha1.hexdigest()
        signature = {
            "client": self.client,
            "nonce": nonce,
            "timestamp": timestamp,
            "signature": hashkey
        }
        return urllib.urlencode(signature)

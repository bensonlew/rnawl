# -*- coding: utf-8 -*-
# __author__ = 'xuting'
from __future__ import division
import web
import json
import datetime
import os
from biocluster.config import Config
from mainapp.models.data_exchange import Identity


class UploadTask(object):
    def POST(self):
        data = web.input()
        print "收到请求，请求的内容为:"
        print data
        code = data.identity
        ip = data.ip
        user = data.user
        mode = data.mode
        test = True if data.mode == 'tsanger' else False
        info = Identity(test=test).get_target_path(code)
        if not info["success"]:
            info["rel_path"] = ""
            return json.dumps(info)
        else:
            upload_info = dict()
            upload_info["code"] = code
            upload_info["request_time"] = datetime.datetime.now()
            upload_info["ip"] = ip
            upload_info["user"] = user
            Identity().add_upload_record(upload_info)
            self.config = Config().get_netdata_config(mode)
            info["abs_path"] = os.path.join(self.config[mode + "_path"], info["rel_path"])
            return json.dumps(info)

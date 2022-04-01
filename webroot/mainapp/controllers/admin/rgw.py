# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
import web
from mainapp.libs.signature import check_sig
import json
from mainapp.config.db import Config


class RgwAction(object):

    @check_sig
    def POST(self):
        data = web.input()
        if not (hasattr(data, "region") and hasattr(data, "bucket")):
            info = {"success": False, "info": "缺少参数!"}
            return json.dumps(info)
        return json.dumps(Config().get_rgw_conn(data.region, data.bucket))
# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import web
import json
import datetime
import os
from biocluster.config import Config
from mainapp.models.mongo.metaasv import Metaasv
from bson.objectid import ObjectId
from types import StringTypes
import re
from mainapp.libs.signature import check_sig
from mainapp.libs.input_check import meta_check
from mainapp.controllers.project.metaasv_controller import MetaasvController


class ConvertLevelAction(MetaasvController):
    """
    此功能用于asv分类学水平实时的调取数据，生成接口
    避免后台存入比较多的数据导致冗余
    """
    def __init__(self):
        super(ConvertLevelAction, self).__init__(instant=True)
        self.metaasv = Metaasv()

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        postArgs = ['level_id', 'submit_location', "asv_id", "task_type"]
        for arg in postArgs:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'parameters missing:%s' % arg}
                return json.dumps(info)
        otu_path = os.path.join(Config().WORK_DIR, "tmp", datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3] + ".otu.xls")
        self.metaasv.export_otu_table_by_level(data.asv_id, otu_path, data.level_id)
        self.metaasv.add_otu_detail(otu_path, data.asv_id, data.level_id)
        info = dict()
        info["success"] = True
        info["info"] = "已成功完成计算"
        return json.dumps(info)
# -*- coding: utf-8 -*-
# __author__ = 'shijin'
import web
import json
import datetime
from mainapp.controllers.project.ref_rna_controller import RefRnaController
from mbio.api.to_file.ref_rna import *
from mainapp.libs.signature import check_sig

class GenesetDeleteAction(RefRnaController):
    def __init__(self):
        super(GenesetDeleteAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        #print data
        default_argu = ['geneset_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        #print(data.geneset_id)
        if self.ref_rna.delete_geneset(data.geneset_id):
            task_info = {"success": True, "info": "任务运行成功."}
            task_info['content'] = {
                'ids': {
                    'id': "",
                    'name': ""
                    }}
            return json.dumps(task_info)
        else:
            task_info = {"success": False, "info": "任务运行失败."}
            return json.dumps(task_info)

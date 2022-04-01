# -*- coding: utf-8 -*-
# __author__ = 'shijin'

import json

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class GenesetDeleteAction(WholeTranscriptomeController):
    def __init__(self):
        super(GenesetDeleteAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        # print data
        default_argu = ['geneset_id', 'task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                variables = []
                variables.append(argu)
                info = {'success': False, 'info': '%s参数缺少!' % argu, 'code': 'C2901201', 'variables': variables}
                return json.dumps(info)
        # print(data.geneset_id)
        if self.whole_transcriptome.delete_geneset(data.geneset_id, data.task_id):
            task_info = {"success": True, "info": "任务运行成功."}
            task_info['content'] = {
                'ids': {
                    'id': "",
                    'name': ""
                }}
            return json.dumps(task_info)
        else:
            task_info = {"success": False, "info": "任务运行失败.", 'code': 'C2901202', 'variables': ''}
            return json.dumps(task_info)

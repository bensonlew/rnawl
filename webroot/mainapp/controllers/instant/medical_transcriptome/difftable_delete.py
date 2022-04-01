# -*- coding: utf-8 -*-
# __author__ = 'shijin'
import web
import json
import datetime
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mbio.api.to_file.medical_transcriptome import *
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId


class DifftableDeleteAction(MedicalTranscriptomeController):
    def __init__(self):
        super(DifftableDeleteAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        self.input_data = data
        #print data
        default_argu = ['diff_id', 'task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                variables = []
                variables.append(argu)
                info = {'success': False, 'info': '%s参数缺少!' % argu, 'code': 'C2901201', 'variables': variables}
                return json.dumps(info)
        params = self.pack_params()
        diff_infos = self.medical_transcriptome.get_diff_table_name(data.diff_id,data.task_id)
        diff_name = diff_infos["name"]
        time_now = datetime.datetime.now()
        name  = "Delete_" + diff_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            task_id=data.task_id,
            name=name,
            version="v1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='cancel interaction task',
            params=params,
            deleted_geneset  = diff_infos["genesets"],
            status="start")
        main_id = self.medical_transcriptome.insert_main_table('sg_difftable_delete', main_info)

        #print(data.geneset_id)
        if self.medical_transcriptome.delete_difftable(data.diff_id, data.task_id):
            query_dict ={'_id': ObjectId(main_id)}
            update_dict = {'status': "end", 'desc': "删除成功"}
            self.medical_transcriptome.update_db_record('sg_difftable_delete', query_dict=query_dict, insert_dict=update_dict)
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

    def pack_params(self):
        params_dict = dict()
        for each in self.input_data:
            if each == "task_type":
                params_dict[each] = int(self.input_data[each])
            else:
                params_dict[each] = self.input_data[each]
        params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
        return params

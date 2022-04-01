# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
import datetime
import string
from random import choice
from bson import ObjectId
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig


class PanCoreAction(MetaasvController):
    """
    Metaasv Pan_core分析
    """
    def __init__(self):
        super(PanCoreAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        postArgs = ['asv_id', 'group_id', 'level_id', "group_detail"]
        for arg in postArgs:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'parameters missing:%s' % arg}
                return json.dumps(info)
        task_name = 'metaasv.report.pan_core'
        task_type = 'workflow'
        time_now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        level_name = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV"]
        main_table_name = 'PanCore' + level_name[int(data.level_id) - 1] + "_" + time_now
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        group_detal_dict = json.loads(data.group_detail)
        specimen_ids = list()
        for v in group_detal_dict.values():
            if len(v) < 2:
                info = {'success': False, 'info': '每个组别至少应该有两个样本！'}
                return json.dumps(info)
            for tmp in v:
                specimen_ids.append(tmp)
        specimen_ids = ",".join(specimen_ids)
        params_json = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'submit_location': data.submit_location,
            'task_type': str(data.task_type)
        }

        unique_id = self.get_unique()
        mongo_data = {
            'project_sn': task_info['project_sn'],
            'task_id': task_info['task_id'],
            'asv_id': ObjectId(data.asv_id),
            'status': 'start',
            'desc': '正在计算Pan OTU',
            'unique_id': unique_id,
            'name': 'Pan_' + time_now,
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "level_id": int(data.level_id),
            "params": json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        }
        main_pan_table_id = self.metaasv.insert_main_table(
            'pan_core', mongo_data)
        update_info = {str(main_pan_table_id): 'pan_core', str(
            main_pan_table_id): 'pan_core'}
        options = {
            "in_otu_table": data.asv_id,
            "group_table": data.group_id,
            'update_info': json.dumps(update_info),
            "group_detail": data.group_detail,
            "samples": self.metaasv.sampleIdToName(task_info['task_id'], specimen_ids),
            "level": int(data.level_id),
            "main_id": str(main_pan_table_id),
        }
        to_file = [
            "metaasv.export_otu_table_by_level(in_otu_table)", "metaasv.export_group_table_by_detail(group_table)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="PanCore/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(PanCoreAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_pan_table_id),
                'name': 'Pan_' + time_now
            }
        }
        return json.dumps(task_info)

    def get_unique(self):
        chars = string.ascii_letters + string.digits
        unique_id = "".join([choice(chars) for i in range(10)])
        collection = self.metaasv.db["pan_core"]
        result = collection.find_one({"unique_id": unique_id})
        if result:
            return self.metaasv.get_unique()
        return unique_id

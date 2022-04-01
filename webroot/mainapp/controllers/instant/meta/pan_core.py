# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import web
import json
import datetime
import string
from random import choice
from bson import ObjectId
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
# from mainapp.models.mongo.public.meta.meta import Meta


class PanCore(MetaController):
    def __init__(self):
        super(PanCore, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        postArgs = ['group_id', 'level_id', 'submit_location', "group_detail"]
        for arg in postArgs:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'parameters missing:%s' % arg}
                return json.dumps(info)
        task_name = 'meta.report.pan_core'
        task_type = 'workflow'
        time_now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        level_name = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"]
        main_table_name = 'PanCore' + level_name[int(data.level_id) - 1] + "_" + time_now  # modified by hongdongxuan 20170323
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!", 'code':'C2202701'}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])
        group_detal_dict = json.loads(data.group_detail)
        specimen_ids = list()
        for v in group_detal_dict.values():
            if len(v) < 2:
                info = {'success': False, 'info': '每个组别至少应该有两个样本！', 'code':'C2202702'}
                return json.dumps(info)
            for tmp in v:
                specimen_ids.append(tmp)
        #if len(specimen_ids) < 5:
        #    info = {'success': False, 'info': '样本数量少于5个,请选择更多的样本！'}
        #    return json.dumps(info)
        specimen_ids = ",".join(specimen_ids)
        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'submit_location': data.submit_location,
            'task_type': data.task_type
        }
        to_file = [
            "meta.export_otu_table_by_level(in_otu_table)", "meta.export_group_table_by_detail(group_table)"]
        unique_id = self.get_unique()
        mongo_data = {
            'project_sn': task_info['project_sn'],
            'task_id': task_info['task_id'],
            'otu_id': ObjectId(data.otu_id),
            'type': 1,  # pan
            'status': 'start',
            'desc': '正在计算Pan OTU',
            'unique_id': unique_id,
            'name': 'Pan_' + time_now,
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "level_id": int(data.level_id),
            "params": json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        }
        main_pan_table_id = self.meta.insert_main_table(
            'sg_otu_pan_core', mongo_data)
        mongo_data['type'] = 2
        mongo_data['desc'] = '正在计算Core OTU'
        mongo_data['name'] = 'Core_' + time_now
        main_core_table_id = self.meta.insert_main_table(
            'sg_otu_pan_core', mongo_data)
        update_info = {str(main_pan_table_id): 'sg_otu_pan_core', str(
            main_core_table_id): 'sg_otu_pan_core'}
        options = {
            "in_otu_table": data.otu_id,
            "group_table": data.group_id,
            'update_info': json.dumps(update_info),
            "group_detail": data.group_detail,
            "samples": self.meta.sampleIdToName(specimen_ids),
            "level": int(data.level_id),
            "main_pan_id": str(main_pan_table_id),
            "main_core_id": str(main_core_table_id)
        }
        self.set_sheet_data(name=task_name, options=options, main_table_name="PanCore/" + main_table_name,
                            module_type=task_type, to_file=to_file)  # modified by hongdongxuan 20170322 在main_table_name前面加上文件输出的文件夹名
        task_info = super(PanCore, self).POST()
        task_info['content'] = {
            'ids': [{
                'id': str(main_pan_table_id),
                'name': 'Pan_' + time_now
            }, {
                'id': str(main_core_table_id),
                'name': 'Core_' + time_now
            }]
        }
        return json.dumps(task_info)

    def get_unique(self):
        chars = string.ascii_letters + string.digits
        unique_id = "".join([choice(chars) for i in range(10)])
        collection = self.meta.db["sg_otu_pan_core"]
        result = collection.find_one({"unique_id": unique_id})
        if result:
            return self.meta.get_unique()
        return unique_id

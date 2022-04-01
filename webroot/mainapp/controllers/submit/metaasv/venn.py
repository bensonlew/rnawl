# -*- coding: utf-8 -*-
# __author__ = 'xuting'  last modify by qindanhua 20170110
import web
import json
import datetime
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON


class VennAction(MetaasvController):
    """
    metaasv Venn图
    """
    def __init__(self):
        super(VennAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        postArgs = ['group_id', 'level_id', "group_detail", 'submit_location']
        for arg in postArgs:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'parameters missing: %s' , "variables":[ arg], "code" : "C2203903"}
                return json.dumps(info)

        task_name = 'metaasv.report.venn'
        task_type = 'workflow'

        my_param = dict()
        my_param['asv_id'] = data.asv_id
        my_param['level_id'] = int(data.level_id)
        my_param['group_id'] = data.group_id
        my_param['group_detail'] = group_detail_sort(data.group_detail)
        my_param["submit_location"] = data.submit_location
        my_param["task_type"] = data.task_type

        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!"}
            return json.dumps(info)

        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        main_table_name = 'Venn_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        group_detal_dict = json.loads(data.group_detail)
        if len(group_detal_dict) < 2:
            info = {"success": False, "info": "进行Venn分析，分组方案的分组类别必须大于等于2且小于等于6！"}
            return json.dumps(info)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ('status', 'start'),
            ('name', main_table_name),
            ('desc', 'venn分析正在计算中'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(my_param, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.metaasv.insert_none_table('venn')
        update_info = {str(main_table_id): 'venn'}

        specimen_ids = list()
        for v in group_detal_dict.values():
            for tmp in v:
                specimen_ids.append(tmp)
        specimen_ids = ",".join(specimen_ids)
        options = {
            "in_otu_table": data.asv_id,
            "update_info": json.dumps(update_info),
            "group_detail": data.group_detail,
            "group_table": data.group_id,
            "samples": self.metaasv.sampleIdToName(task_info['task_id'], specimen_ids),
            "level": data.level_id,
            "asv_id": str(data.asv_id),
            "main_id": str(main_table_id),
            "main_table_data": SON(mongo_data)
        }
        to_file = ["metaasv.export_otu_table_by_level(in_otu_table)", "metaasv.export_group_table_by_detail(group_table)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="Venn/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(VennAction, self).POST()
        if task_info['success']:
            task_info['content'] = {
                'ids': {
                    'id': str(main_table_id),
                    'name': main_table_name
                }}
        return json.dumps(task_info)

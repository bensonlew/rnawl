# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
from mainapp.controllers.project.metaasv_controller import MetaasvController
import datetime
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON


class RankAbundanceAction(MetaasvController):
    """
    Metaasv Rank_abundance曲线分析
    """
    def __init__(self):
        super(RankAbundanceAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        params_name = ['asv_id', 'level_id', 'group_id', 'group_detail']
        for param in params_name:
            if not hasattr(data, param):
                info = {"success": False, "info": "parameters missing:%s" % param}
                return json.dumps(info)
        if int(data.level_id) not in range(1, 10):
            variables = []
            variables.append(data.level_id)
            info = {"success": False, "info": "level{}不在规定范围内!".format(data.level_id), }
            return json.dumps(info)
        my_param = dict()
        my_param['asv_id'] = data.asv_id
        my_param['level_id'] = int(data.level_id)
        my_param['submit_location'] = data.submit_location
        my_param['task_type'] = str(data.task_type)
        my_param['group_detail'] = group_detail_sort(data.group_detail)
        my_param['group_id'] = data.group_id
        task_name = 'metaasv.report.rank_abundance'
        task_type = 'workflow'

        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        main_table_name = 'Rank_abundance' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ("level_id", int(data.level_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(my_param, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.metaasv.insert_none_table('rank_abundance')
        update_info = {str(main_table_id): 'rank_abundance'}

        options = {
            "update_info": json.dumps(update_info),
            "otu_table": data.asv_id,
            "level": int(data.level_id),
            "main_id": str(main_table_id),
            "group_table": data.group_id,
            'group_detail': data.group_detail,
            "main_table_data": SON(mongo_data)
        }
        to_file = ["metaasv.export_otu_table_by_level(otu_table)", "metaasv.export_group_table_by_detail(group_table)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="Rank_abundance/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(RankAbundanceAction, self).POST()
        if task_info['success']:
            task_info['content'] = {
                'ids': {
                    'id': str(main_table_id),
                    'name': main_table_name
                }}
        return json.dumps(task_info)

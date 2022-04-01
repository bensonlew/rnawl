# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
from mainapp.controllers.project.metaasv_controller import MetaasvController
import datetime
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON


class RarefactionAction(MetaasvController):
    """
    metaasv 稀释曲线分析
    """
    ESTIMATORS = ['ace', 'chao', 'coverage', 'pd', 'shannon', 'shannoneven', 'simpson', 'simpsoneven', 'sobs']

    def __init__(self):
        super(RarefactionAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        params_name = ['asv_id', 'level_id', 'index_type', 'submit_location', 'group_id', 'group_detail']
        for param in params_name:
            if not hasattr(data, param):
                info = {"success": False, "info": "parameters missing:%s" % param}
                return json.dumps(info)
        for index in data.index_type.split(','):
            if index not in self.ESTIMATORS:
                variables = []
                variables.append(index)
                info = {"success": False, "info": "指数类型不正确{}".format(index)}
                return json.dumps(info)
        if int(data.level_id) not in range(1, 10):
            variables = []
            variables.append(data.level_id)
            info = {"success": False, "info": "level{}不在规定范围内!".format(data.level_id)}
            return json.dumps(info)
        my_param = dict()
        my_param['asv_id'] = data.asv_id
        my_param['level_id'] = int(data.level_id)
        my_param['index_type'] = data.index_type
        my_param['submit_location'] = data.submit_location
        my_param['task_type'] = str(data.task_type)
        my_param['group_detail'] = group_detail_sort(data.group_detail)
        my_param['group_id'] = data.group_id
        task_name = 'metaasv.report.rarefaction'
        task_type = 'workflow'

        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        main_table_name = 'Rarefaction_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ("index_type", str(data.index_type)),
            ("level_id", int(data.level_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(my_param, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.metaasv.insert_none_table('rarefaction')
        update_info = {str(main_table_id): 'rarefaction'}

        options = {
            "update_info": json.dumps(update_info),
            "asv_id": data.asv_id,
            "otu_table": data.asv_id,
            "indices": data.index_type,
            "level": data.level_id,
            "rare_id": str(main_table_id),
            "group_detail": data.group_detail,
            'group': data.group_id,
            "main_table_data": SON(mongo_data)
        }
        if hasattr(data, "freq"):
            options["freq"] = data.freq

        to_file = ["metaasv.export_otu_table_by_detail(otu_table)", "metaasv.export_group_table_by_detail(group)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="Rarefaction/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(RarefactionAction, self).POST()
        if task_info['success']:
            task_info['content'] = {
                'ids': {
                    'id': str(main_table_id),
                    'name': main_table_name
                }}
        return json.dumps(task_info)

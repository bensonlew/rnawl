# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
from mainapp.controllers.project.meta_controller import MetaController
# from mainapp.libs.signature import check_sig
# from mainapp.models.workflow import Workflow
# from mainapp.models.mongo.meta import Meta
# from mainapp.models.mongo.estimator import Estimator
import datetime
# from mainapp.libs.param_pack import GetUploadInfo
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON


class RarefactionAction(MetaController):
    """

    """
    ESTIMATORS = ['ace', 'bootstrap', 'chao', 'coverage', 'default', 'heip', 'invsimpson', 'jack', 'npshannon',
                  'shannon', 'shannoneven', 'simpson', 'simpsoneven', 'smithwilson', 'sobs']

    def __init__(self):
        super(RarefactionAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        params_name = ['otu_id', 'level_id', 'index_type', 'freq', 'submit_location', 'group_id', 'group_detail']
        for param in params_name:
            if not hasattr(data, param):
                info = {"success": False, "info": "parameters missing:%s" % param}
                return json.dumps(info)
        for index in data.index_type.split(','):
            if index not in self.ESTIMATORS:
                variables = []
                variables.append(index)
                info = {"success": False, "info": "指数类型不正确{}".format(index), 'code':'C2203301', 'variables':variables}
                return json.dumps(info)
        if int(data.level_id) not in range(1, 10):
            variables = []
            variables.append(data.level_id)
            info = {"success": False, "info": "level{}不在规定范围内!".format(data.level_id), 'code':'C2203302', 'variables':variables}
            return json.dumps(info)
        my_param = dict()
        my_param['otu_id'] = data.otu_id
        my_param['level_id'] = int(data.level_id)
        # my_param['indices'] = data.index_type
        my_param['freq'] = int(data.freq)
        #sort_index = data.index_type.split(',')
        #sort_index.sort()
        #sort_index = ','.join(sort_index)
        #my_param['index_type'] = sort_index
        my_param['index_type'] = data.index_type
        my_param['submit_location'] = data.submit_location
        my_param['task_type'] = data.task_type
        my_param['group_detail'] = group_detail_sort(data.group_detail)
        my_param['group_id'] = data.group_id
        # params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        # otu_info = self.meta.get_otu_table_info(data.otu_id)
        task_name = 'meta.report.rarefaction'
        task_type = 'workflow'

        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!", 'code':'C2203303'}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])
        main_table_name = 'Rarefaction_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),
            ("level_id", int(data.level_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(my_param, sort_keys=True, separators=(',', ':')))
        ]
        # main_table_id = Estimator().add_rare_collection(data.level_id, params, data.otu_id, main_table_name)
        main_table_id = self.meta.insert_none_table('sg_alpha_rarefaction_curve')
        update_info = {str(main_table_id): 'sg_alpha_rarefaction_curve'}

        options = {
            "update_info": json.dumps(update_info),
            "otu_id": data.otu_id,
            "otu_table": data.otu_id,
            "indices": data.index_type,
            "level": data.level_id,
            "freq": data.freq,
            "rare_id": str(main_table_id),
            "group_detail": data.group_detail,
            "main_table_data": SON(mongo_data)
        }

        to_file = "meta.export_otu_table_by_detail(otu_table)"
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

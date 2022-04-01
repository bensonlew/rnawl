# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
import web
import json
import random
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId
from mainapp.libs.param_pack import *
from biocluster.config import Config
from mainapp.models.mongo.submit.denovo_rna.denovo_express import DenovoExpress
import types
from mainapp.models.mongo.meta import Meta
from mainapp.models.workflow import Workflow
from mainapp.controllers.project.denovo_controller import DenovoController


class DiffExpress(DenovoController):
    def __init__(self):
        super(DiffExpress, self).__init__()

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        #print data
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        my_param = dict()
        my_param['express_id'] = data.express_id
        my_param['group_detail'] = group_detail_sort(data.group_detail)
        my_param['group_id'] = data.group_id
        my_param['control_id'] = data.control_id
        my_param['ci'] = data.ci
        my_param['submit_location'] = data.submit_location
        params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        express_info = Meta(db=self.mongodb).get_main_info(data.express_id, 'sg_denovo_express')
        if express_info:
            main_table_name = "DiffExpress_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            task_id = express_info["task_id"]
            project_sn = express_info["project_sn"]
            diff_express_id = DenovoExpress().add_express_diff(params=params, samples=None, compare_column=None, project_sn=project_sn, task_id=task_id, express_id=data.express_id)
            update_info = {str(diff_express_id): "sg_denovo_express_diff"}
            update_info = json.dumps(update_info)
            options = {
                "express_file": data.express_id,
                "update_info": update_info,
                "group_id": data.group_id,
                "control_file": data.control_id,
                "ci": data.ci,
                "diff_express_id": str(diff_express_id)
            }
            to_file = ["denovo.export_express_matrix(express_file)", "denovo.export_control_file(control_file)"]
            if data.group_id != 'all':
                to_file.append("denovo.export_group_table_by_detail(group_file)")
                options.update({
                    "group_file": data.group_id,
                    "group_detail": data.group_detail,
                })
            self.set_sheet_data(name='denovo_rna.report.diff_express', options=options, main_table_name=main_table_name, module_type='workflow', to_file=to_file, main_id=diff_express_id, collection_name="sg_denovo_express_diff")
            task_info = super(DiffExpress, self).POST()
            task_info['content'] = {'ids': {'id': str(diff_express_id), 'name': main_table_name}}
            #print task_info
            return json.dumps(task_info)
        else:
            info = {"success": False, "info": "express_id不存在，请确认参数是否正确！!"}
            return json.dumps(info)

    def check_options(self, data):
        """
        检查网页端传进来的参数是否正确
        """
        params_name = ['express_id', 'ci', 'group_detail', 'group_id', 'control_id', 'submit_location']
        success = []
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数!")
        if float(data.ci) >= 1 or float(data.ci) <= 0:
            success.append("显著性水平ci不在范围内")
        for ids in [data.express_id, data.group_id, data.control_id]:
            ids = str(ids)
            #print type(ids)
            if not isinstance(ids, ObjectId) and not isinstance(ids, types.StringTypes):
                success.append("传入的id：{}不是一个ObjectId对象或字符串类型".format(ids))
        return success

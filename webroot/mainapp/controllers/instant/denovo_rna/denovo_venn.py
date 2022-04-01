# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
import web
import json
from mainapp.controllers.core.basic import Basic
from mainapp.libs.param_pack import group_detail_sort
from mainapp.models.mongo.meta import Meta
from mainapp.controllers.project.denovo_controller import DenovoController


class DenovoVenn(DenovoController):
    def __init__(self):
        super(DenovoVenn, self).__init__(instant=True)

    def POST(self):
        data = web.input()
        self.data = data
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        denovo = Meta(db=self.mongodb)  # 需指定mongodb; self.mongodb为tsanger_rna/sanger_rna
        express_info = denovo.get_main_info(data.express_id, 'sg_denovo_express')
        if not express_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_id = express_info["task_id"]
        project_sn = express_info["project_sn"]
        task_name = 'denovo_rna.report.venn'
        my_param = dict()
        my_param['express_id'] = data.express_id
        my_param['group_id'] = data.group_id
        my_param['group_detail'] = group_detail_sort(data.group_detail)
        my_param["submit_location"] = data.submit_location
        my_param["task_type"] = data.task_type
        self.params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        main_table_name = "Venn_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = [
            ('project_sn', project_sn),
            ('task_id', task_id),
            ('express_id', data.express_id),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(my_param, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = denovo.insert_main_table('sg_denovo_venn', mongo_data)
        update_info = {str(main_table_id): 'sg_denovo_venn'}
        self.options = {
            "express_file": data.express_id,
            "group_table": data.group_id,
            "express_id": str(self._mainTableId),
            "group_detail": data.group_detail,
            "update_info": json.dumps(update_info),
            "main_id": str(main_table_id)
        }
        to_file = ["denovo.export_express_matrix(express_file)", "denovo.export_group_table_by_detail(group_table)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type='workflow', to_file=to_file, main_id=main_table_id, collection_name="sg_denovo_venn")
        task_info = super(DenovoVenn, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        #print self.return_msg
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传进来的参数是否正确
        """
        params_name = ['group_id', 'express_id', "group_detail", 'submit_location']
        success = []
        table_dict = json.loads(data.group_detail)
        #print "收到请求, 请求的内容为："
        #print data
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数!")
        if not isinstance(table_dict, dict):
            success.append("传入的table_dict不是一个字典")
        return success

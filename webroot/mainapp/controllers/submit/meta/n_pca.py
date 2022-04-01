# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'
import web
import json
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.models.mongo.meta import Meta
from mainapp.libs.signature import check_sig
from bson import ObjectId
import datetime



class NPcaAction(MetaController):
    def __init__(self):
        super(NPcaAction, self).__init__(instant=False)
        # super(NPca, self).__init__()

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['otu_id', 'level_id', 'submit_location', 'group_detail', 'group_id', 'second_group_id', 'second_group_detail', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'meta.report.n_pca'
        task_type = 'workflow'
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!", 'code':'C2202401'}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])

        if isinstance(data.second_group_detail, dict):
            table_dict = data.second_group_detail
        else:
            table_dict = json.loads(data.second_group_detail)
        if not isinstance(table_dict, dict):
            raise Exception("请选择二级分组！二级分组与一级分组选择的样本请保持一致，至少选择6个分组，每组至少3个样本,参数%s必须是字典格式" %("second_group_detail"))

        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'second_group_detail': group_detail_sort(data.second_group_detail),
            'second_group_id': data.second_group_id,
            'submit_location': data.submit_location,
            'task_type': data.task_type
            }
        main_table_name = 'NPca_' + data.level_id + \
            '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        # main_table_name = 'NPca_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),
            #('table_type', 'dist'),
            #('tree_type', 'cluster'),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.meta.insert_main_table('sg_npca', mongo_data)
        update_info = {str(main_table_id): 'sg_npca'}
        options = {
            'otu_table': data.otu_id,
            'otu_id': data.otu_id,
            'level': int(data.level_id),
            'second_group_table': data.second_group_id,
            'second_group_detail': data.second_group_detail,
            'group_table': data.group_id,
            'group_detail': data.group_detail,
            'update_info': json.dumps(update_info),
            #'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'n_pca_id': str(main_table_id)
            }
        #to_file = 'meta.export_otu_table_by_detail(otu_table)'
        to_file = ["meta.export_otu_table_by_detail(otu_table)", "meta.export_cascading_table_by_detail(second_group_table)","meta.export_group_table_by_detail(group_table)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="NPca/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(NPcaAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        print(self.return_msg)
        return json.dumps(task_info)
        # return json.dumps({'success': True, 'info': 'shenghe log'})

# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import SON


class TreeExtendAction(MetaasvController):
    """
    Metaasv 个性化系统发育树分析
    """
    def __init__(self):
        super(TreeExtendAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['asv_id', 'level_id', 'group_id', 'group_detail', 'submit_location','method','select_spe','upload_spe', 'file_dir_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' , "variables":[ argu]}
                return json.dumps(info)

        task_name = 'metaasv.report.external_tree'
        task_type = 'workflow'  # 可以不配置
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        params = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'submit_location': data.submit_location,
            'task_type': data.task_type,
            'method' : data.method,
            'select_spe' : data.select_spe,
            'upload_spe' : data.upload_spe,
            'file_dir_id' : data.file_dir_id
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'PlotTree_' + \
            '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params),
            ("tree_type", 'extend')

        ]
        main_table_id = self.metaasv.insert_none_table('personal_phylo_tree')
        update_info = {str(main_table_id): 'personal_phylo_tree'}
        options = {
                   'asv_id': data.asv_id,
                   'level': int(data.level_id),
                   'group_id': data.group_id,
                   'update_info': json.dumps(update_info),
                   'group_detail': data.group_detail,
                   # 'params': params,
                   'main_id': str(main_table_id),
                   'main_table_data': SON(mongo_data),
                   'upload_spe' : data.upload_spe ,
                   'method' : data.method,
                   'select_spe_fasta' : data.select_spe  #根据select_spe 筛选出的序列
                   }
        to_file = ['metaasv.export_rep_seq_by_samples_abund(select_spe_fasta)']
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="PlotTree/" + main_table_name,
                            module_type=task_type,
                            to_file=to_file)
        task_info = super(TreeExtendAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

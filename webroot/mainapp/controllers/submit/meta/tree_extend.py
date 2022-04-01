# -*- coding: utf-8 -*-

import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import SON


class TreeExtendAction(MetaController):
    def __init__(self):
        super(TreeExtendAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['otu_id', 'level_id', 'group_id', 'group_detail', 'submit_location','method','select_spe','upload_spe', 'file_dir_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' , "variables":[ argu], "code" : "C2204401"}
                return json.dumps(info)
        # if int(data.level_id) < int(data.color_level_id):
        #     info = {'success': False, 'info': '颜色设置水平必须高于进化树绘制水平!', 'code':'C2203101'}
        #     return json.dumps(info)
        # if int(data.topn) < 2:
        #     if int(data.topn) != 0:
        #         info = {'success': False, 'info': '至少选择丰度高的物种2个及以上', 'code':'C2203102'}
        #         return json.dumps(info)
        task_name = 'meta.report.external_tree'
        task_type = 'workflow'  # 可以不配置
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!", 'code':'C2204402'}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])
        params = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            #'color_level_id': int(data.color_level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'submit_location': data.submit_location,
            'task_type': '',
            'method' : data.method,
            'select_spe' : data.select_spe,
            'upload_spe' : data.upload_spe,
            'file_dir_id': data.file_dir_id,
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'Personality_PlotTree_' + \
            '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params),
            ("tree_type", 'extend')

        ]
        main_table_id = self.meta.insert_none_table('sg_phylo_tree')
        update_info = {str(main_table_id): 'sg_phylo_tree'}
        options = {
                    #'otu_table': data.otu_id,
                   'otu_id': data.otu_id,
                   'level': int(data.level_id),
                   #'color_level_id': int(data.color_level_id),
                   #'sample_group': data.group_id,
                   'group_id': data.group_id,
                   'update_info': json.dumps(update_info),
                   'group_detail': data.group_detail,
                   'params': params,
                   'main_id': str(main_table_id),
                   #'topN': int(data.topn),
                   'main_table_data': SON(mongo_data),
                   'upload_spe' : data.upload_spe ,
                   'method' : data.method,
                   'select_spe_fasta' : data.select_spe  #根据select_spe 筛选出的序列
                   }
        #to_file = ['meta.export_otu_table_by_detail(otu_table)', 'meta.export_group_table_by_detail(sample_group)']
        to_file = ['meta.export_rep_seq_by_samples_abund(select_spe_fasta)']
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="Personality_PlotTree/" + main_table_name,
                            module_type=task_type,
                            to_file=to_file)
        task_info = super(TreeExtendAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

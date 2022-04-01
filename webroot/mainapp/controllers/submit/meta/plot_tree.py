# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import SON


class PlotTreeAction(MetaController):
    def __init__(self):
        super(PlotTreeAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['otu_id', 'level_id', 'color_level_id', 'group_id', 'group_detail', 'submit_location', 'topn', 'method']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' , "variables":[ argu], "code" : "C2203104"}
                return json.dumps(info)
        if int(data.level_id) < int(data.color_level_id):
            info = {'success': False, 'info': '颜色设置水平必须高于进化树绘制水平!', 'code':'C2203101'}
            return json.dumps(info)
        if int(data.topn) < 2:
            if int(data.topn) != 0:
                info = {'success': False, 'info': '至少选择丰度高的物种2个及以上', 'code':'C2203102'}
                return json.dumps(info)
        task_name = 'meta.report.plot_tree'
        task_type = 'workflow'  # 可以不配置
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!", 'code':'C2203103'}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])
        task_id = otu_info['task_id']
        method = data.method   #'NJ'
        if method not in ['NJ', 'ML', 'MP']:
            info = {'success': False, 'info': "Value of  method parameters not in ['NJ', 'ML', 'MP']" , "code" : "C2203105"}
            return json.dumps(info)

        if method in ['MP']:
            if int(data.topn) > 500:
                info = {'success': False, 'info': "top num can not more than 500", "code" : "C2203109"}
                return  json.dumps(info)
        elif method in ['NJ', 'ML']:
            if int(data.topn) > 500:
                info = {'success': False, 'info': "top num can not more than 500", "code" : "C2203110"}
                return  json.dumps(info)

        run_all_tree_condition = False
        if run_all_tree_condition :
            # ##判断如果没有总树，则多传2个参数用来先建总树
            tree_type_map = {
                'NJ' : 'phylo',
                'ML' : 'phylo_ml',
                'MP' : 'phylo_mp'
            }

            tree_type = tree_type_map[method]

            # 判断是否要运行总树
            run_all_tree = ''
            otu_table = self.meta.get_mongo_common('sg_otu',{"task_id":task_id, 'name':"OTU_Taxon_Origin"})
            ori_otu_table_id= ''
            if otu_table:
                ori_otu_table_id =  otu_table['_id']
                search = {"task_id":task_id, 'table_id': ori_otu_table_id, "table_type": 'otu', "tree_type": tree_type}
                if  self.meta.get_mongo_common('sg_newick_tree',search):
                   run_all_tree = 'end'

            if run_all_tree =='':
                search_is_run = {'task_id':task_id, "run_all":'yes'}
                res_search_is_run = self.meta.get_mongo_common('sg_phylo_tree',search_is_run)
                if res_search_is_run :
                    if res_search_is_run['status'] in ['end']:
                        run_all_tree = 'end'
                    elif res_search_is_run['status'] in ['start']:
                        run_all_tree = 'start'
                    else:
                        run_all_tree = 'not_start'
                else:
                    run_all_tree = 'not_start'

        # run_all_tree为start值时，将重头建树
        else:
            run_all_tree = 'start'


        params = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'color_level_id': int(data.color_level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'submit_location': data.submit_location,
            'task_type': data.task_type,
            'topn': data.topn,
            'method' : method
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'PlotTree_' + \
            '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_id),
            ('otu_id', ObjectId(data.otu_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]

        if method == 'ML':
            mongo_data.append(("software","IQ-TREE"))
        else:
            mongo_data.append(("software","MEGA"))

        if run_all_tree == 'not_start':
            mongo_data.append(("run_all","yes"))

        main_table_id = self.meta.insert_none_table('sg_phylo_tree')
        update_info = {str(main_table_id): 'sg_phylo_tree'}
        options = {
                   'otu_id': data.otu_id,
                   'level': int(data.level_id),
                   'color_level_id': int(data.color_level_id),
                   'group_id': data.group_id,
                   'update_info': json.dumps(update_info),
                   'group_detail': data.group_detail,
                   'params': params,
                   'main_id': str(main_table_id),
                   'topN': int(data.topn),
                   'main_table_data': SON(mongo_data),
                   'task_id' : task_info['task_id'],
                   'method' : method
                   }

        # >>>zouguanqing 201910 V4
        if run_all_tree == 'not_start':
            if ori_otu_table_id == '':
                info = {'success': False, 'info': "原始的otu表没有找到，无法重头建总树" , "code" : "C2203106"}
                return json.dumps(info)

            options['ori_otu_id'] = str(ori_otu_table_id)
            options['otu_table'] = data.otu_id
            options['sample_group'] = data.group_id
            options['run_tree'] = 'all'
            options['otu_seq'] = ori_otu_table_id
            to_file = ['meta.export_otu_table_by_detail(otu_table)', 'meta.export_group_table_by_detail(sample_group)']
            to_file.append('meta.export_otu_seqs(otu_seq)')

        elif run_all_tree == 'start':
            options["run_tree"] = 'part'
            options['otu_seq'] = data.otu_id
            to_file = ['meta.get_seq_by_select_samples_pick_pre_nums(otu_seq)']
        else:
            options['otu_table'] = data.otu_id
            options['sample_group'] = data.group_id
            options["run_tree"] = 'pick'
            to_file = ['meta.export_otu_table_by_detail(otu_table)', 'meta.export_group_table_by_detail(sample_group)']

        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="PlotTree/" + main_table_name,
                            module_type=task_type,
                            to_file=to_file)
        task_info = super(PlotTreeAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.models.mongo.metagenomic import Metagenomic
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from mainapp.libs.signature import check_sig
from .comm_creat import CommCreatAction

class TwoGroupAction(MetagenomicController):
    def __init__(self):
        super(TwoGroupAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['anno_type', 'ci', 'correction', 'correction_ci', 'coverage', 'geneset_id', 'group_detail', 'group_id', 'method', 'submit_location', 'test', 'tail_type', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % argu}
                return json.dumps(info)
        if float(data.ci) > 1 or float(data.ci) < 0 :
            info = {"success": False, "info": "参数错误：显著性水平参数值错误, 应该在[0,1]", "code":"C2402401", "variables":""}
            return json.dumps(info)
        if float(data.correction_ci) > 1 or float(data.correction_ci) < 0 :
            info = {"success": False, "info": "参数错误：多重检验校正显著性水平参数值错误, 应该在[0,1]", "code":"C2402402", "variables":""}
            return json.dumps(info)
        if data.correction not in  ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of correction (%s)" % data.correction}
            return json.dumps(info)
        if data.tail_type not in ["two.side", "greater", "less"]:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of tail_type (%s)" % data.tail_type }
            return json.dumps(info)
        if data.test not in ["chi", "fisher", "kru_H", "mann", "anova", "student", "welch", "signal"]:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of test (%s)" % data.test }
            return json.dumps(info)
        if float(data.coverage) not in [0.90, 0.95, 0.98, 0.99, 0.999]:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of coverage (%s)" % data.coverage}
        add_personal = CommCreatAction()
        info = add_personal.judge_database(data)
        if info != "":
            return json.dumps(info)
        table_dict = json.loads(data.group_detail)
        if len(table_dict) != 2 or data.group_id == 'all':
            info = {"success": False, "info": "分析只适用于两个分组的情况！", "code":"C2402403", "variables":""}
            return  json.dumps(info)
        if not isinstance(table_dict, dict):
            info = {"success": False, "info": "PARAMETERS ERROR: wrong type of group_detail, dict expected!"}
            return  json.dumps(info)
        task_name = 'metagenomic.report.two_group'
        module_type = 'workflow'
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        task_info = metagenomic.get_task_info(geneset_info['task_id'])
        params_json = {
            'anno_type': data.anno_type,
            'geneset_id': data.geneset_id,
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'method': data.method,
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            'ci': float(data.ci),
            'correction': data.correction,
            # 'correction_ci': float(data.correction_ci),
            'correction_ci': float(data.correction_ci) if float(data.correction_ci) != 1 else 1,  # 为1时取整数 by ghd @ 20181015
            'tail_type': data.tail_type,
            'test': data.test,
            'coverage': float(data.coverage),
        }
        #group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('desc', 'processing'),
            ('anno_type', data.anno_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('test_type', 'two_group'),
        ]
        level = ""
        level_name = ""
        if data.anno_type not in ["gene"]:
            params_json['anno_id'] = data.anno_id
            params_json['level_id'] = int(data.level_id)  # 此处是level缩写，需要再转义成真实level name
            level = data.level_id
            level_name = self.level_id(data.level_id)
            if hasattr(data, 'second_level'):
                params_json['second_level'] = data.second_level
            #main_table_name = 'DiffTwoGroup_' + data.anno_type.upper() + '_' + level_name + '_' + \
            #              datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        else:
            main_table_name = 'DiffTwoGroup_GENE_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        #mongo_data.append(('name', main_table_name))
        #mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        #main_table_id = self.metagenomic.insert_main_table('metastat', mongo_data)
        #update_info = {str(main_table_id): 'metastat'}
        [geneset_table, gene_list] = metagenomic.export_geneset_table(data.geneset_id, data.method)
        options = {
            #'update_info': json.dumps(update_info),
            # 'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            # 'group_id': data.group_id,
            # 'geneset_id': data.geneset_id,
            'group_detail': data.group_detail,
            # 'anno_type': data.anno_type,
            'group_file': data.geneset_id,  # 为修改id转样品功能
            'test': data.test,
            'correction': data.correction,
            'ci': float(data.ci),
            'type': data.tail_type,
            'coverage': data.coverage,
            #'main_id': str(main_table_id),
            # 'group_name': metagenomic.get_group_name(data.group_id),
            'group_name': "group_name",
            # 'group_name': ','.join(table_dict.keys()),
            # 'main_table_data': SON(mongo_data),
            'geneset_table': geneset_table,
        }
        to_file =[]
        #个性化数据库调用
        add_personal = CommCreatAction()
        data,options, params_json,info, to_file = add_personal.add_personal(data,options,params_json,geneset_table,to_file)
        if info!= "":
            return json.dumps(info)
        #个性化数据库用调用
        if data.anno_type not in ['gene']:
            options['level_id'] = self.level_id(data.level_id)
            anno_info = metagenomic.get_anno_info(data.anno_id, data.anno_type)
            # options['anno_id'] = data.anno_id
            options['anno_table'] = self.use_s3(anno_info['anno_file'])
            if data.anno_type not in 'nr':
                options['lowest_level'] = self.use_s3(anno_info['lowest_level'])
            if hasattr(data, 'second_level'):
                if data.second_level != "":
                    options['level_type_name'] = self.level_convert(data.second_level, data.level_id)
            main_table_name = 'DiffTwoGroup_' + data.anno_type.upper() + '_' + level_name + '_' + \
                          datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            main_table_name = main_table_name.replace(" ","_")
        else:
            options['gene_list'] = gene_list
        if data.test == "signal":
            options['group_file'] = ':'.join([data.geneset_id, data.group_id])
            to_file.append('metagenomic.export_group_table_for_signal(group_file)')
        else:
            to_file.append('metagenomic.export_group_table_by_detail(group_file)')
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.metagenomic.insert_main_table('metastat', mongo_data)
        update_info = {str(main_table_id): 'metastat'}
        options['update_info'] = json.dumps(update_info)
        options['main_id'] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            task_id=task_info['task_id'],
                            project_sn=task_info['project_sn'], module_type=module_type, params=params_json,
                            to_file=to_file)
        task_info = super(TwoGroupAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

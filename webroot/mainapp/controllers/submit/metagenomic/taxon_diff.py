# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'

import web
import json
import datetime
from bson import ObjectId
from mbio.packages.metagenomic.id_convert import name2id
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.metagenomic_controller import MetagenomicController


class TaxonDiffAction(MetagenomicController):
    def __init__(self):
        super(TaxonDiffAction, self).__init__(instant=False)

    def POST(self):
        data = web.input()
        print data
        default_argu = ['anno_type', 'anno_id', 'level_id', 'submit_location',
                        'group_id', 'group_detail', 'method', 'correction']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s!" % argu, "code": "C2401103"}
                return json.dumps(info)
        if hasattr(data, "post_hoc"):
            test_type = "multiple"
            if data.method not in ["kru_H", "anova"]:
                info = {"success": False, "info": "PARAMETERS ERROR: wrong value of method (%s)" % data.method}
        else:
            test_type = "two_group"
            tail_type = data.tail_type
            if tail_type == "two_side":
                tail_type = "two.side"
        #     if data.tail_type not in ["two.side", "greater", "less"]:
        #         info = {"success": False, "info": "PARAMETERS ERROR: wrong value of tail_type (%s)" % data.tail_type}
        #     if float(data.ci) > 1 or float(data.ci) < 0 :
        #         info = {"success": False, "info": "参数错误：显著性水平参数值错误, 应该在[0,1]", "code":"C2402401", "variables":""}
        #         return json.dumps(info)
        #     if float(data.ci_levle) > 1 or float(data.correction_ci) < 0 :
        #         info = {"success": False, "info": "参数错误：多重检验校正显著性水平参数值错误, 应该在[0,1]", "code":"C2402402", "variables":""}
        #         return json.dumps(info)
        #     if data.correction not in  ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]:
        #         info = {"success": False, "info": "PARAMETERS ERROR: wrong value of correction (%s)" % data.correction}
        #         return json.dumps(info)
        #     if data.tail_type not in ["two.side", "greater", "less"]:
        #         info = {"success": False, "info": "PARAMETERS ERROR: wrong value of tail_type (%s)" % data.tail_type }
        #         return json.dumps(info)
        #     if data.method not in ["chi", "fisher", "kru_H", "mann", "anova", "student", "welch", "signal"]:
        #         info = {"success": False, "info": "PARAMETERS ERROR: wrong value of method (%s)" % data.method }
        #         return json.dumps(info)
        #     if float(data.ci_level) not in [0.90, 0.95, 0.98, 0.99, 0.999]:
        #         info = {"success": False, "info": "PARAMETERS ERROR: wrong value of ci_level (%s)" % data.coverage}
        #         return json.dumps(info)

        task_name = 'metagenomic.report.taxon_diff'
        module_type = 'workflow'
        task_type = 2
        anno_info = self.metagenomic.common_find_one("anno_" + data.anno_type,
                                                     {"_id": ObjectId(data.anno_id)})
        name_id = name2id(anno_info["task_id"], type="task")
        name_id = {k: v for k, v in name_id.items() if v in data.group_detail}
        params_json = {
            'anno_type': data.anno_type,
            'anno_id': data.anno_id,
            'level_id': data.level_id,
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'method': data.method,
            'correction': data.correction,
            'task_type': task_type,
            'submit_location': data.submit_location
        }
        if test_type == "two_group":
            params_json["ci_level"] = data.ci_level
            params_json["tail_type"] = data.tail_type
            if hasattr(data, "pair_id"):
                params_json["pair_id"] = data.pair_id
        else:
            params_json["post_hoc"] = data.post_hoc
            params_json["post_hoc_level"] = data.post_hoc_level
        mongo_data = [
            ('project_sn', anno_info['project_sn']),
            ('task_id', anno_info['task_id']),
            ('status', 'start'),
            ('desc', 'processing'),
            ('anno_type', data.anno_type),
            ('test_type', test_type),
            ('specimen_list', ','.join(name_id.values())),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        level_lab = ['d', 'k', 'p', 'c', 'o', 'f', 'g', 's']
        col = level_lab[int(data.level_id) - 1] + '__'
        to_file = ['metagenomic.export_group_by_detail_taxon(group)']
        options = {
            'table': anno_info["anno_file"],
            'level_id': data.level_id,
            'name2id': json.dumps(name_id),
            'col': col,
            'method': data.method,
            'correction': data.correction,
            'group': data.group_detail,
            'group_detail': data.group_detail,
            'test_type': test_type,
        }
        if test_type == "two_group":
            a_name = "DiffTwo"
            options["ci_level"] = data.ci_level
            options["tail_type"] = tail_type
            if hasattr(data, "pair_id"):
                options["paired"] = data.pair_id
                to_file.append('metagenomic.export_signal_pair_sample(paired)')
        else:
            a_name = "DiffMul"
            options["post_hoc"] = data.post_hoc
            options["post_hoc_level"] = data.post_hoc_level
        main_table_name = a_name + "_{}_{}_{}".format(data.anno_type.title(), col[0],
                                                      datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.metagenomic.insert_main_table('taxon_diff', mongo_data)
        update_info = {str(main_table_id): 'taxon_diff'}
        options['main_id'] = str(main_table_id)
        options['main_col'] = "taxon_diff"
        options['update_info'] = json.dumps(update_info)
        # self.set_sheet_data(name=task_name, options=options,
        #                     main_table_name=a_name + "_reads/" + main_table_name.strip().split("_")[0] + '/' + main_table_name,
        #                     task_id=anno_info['task_id'], project_sn=anno_info['project_sn'],
        #                     module_type=module_type, params=params_json, to_file=to_file)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=a_name + "_reads/" + main_table_name,
                            task_id=anno_info['task_id'], project_sn=anno_info['project_sn'],
                            module_type=module_type, params=params_json, to_file=to_file)
        task_info = super(TaxonDiffAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

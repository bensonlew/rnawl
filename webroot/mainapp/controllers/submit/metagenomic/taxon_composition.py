# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'

import web
import json
import datetime
from bson import ObjectId
from mbio.packages.metagenomic.id_convert import name2id
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.metagenomic_controller import MetagenomicController


class TaxonCompositionAction(MetagenomicController):
    def __init__(self):
        super(TaxonCompositionAction, self).__init__(instant=False)

    def POST(self):
        data = web.input()
        print data
        default_argu = ['anno_type', 'anno_id', 'submit_location', 'group_id', 'group_detail',
                        'group_method', 'graphic_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s!" % argu, "code" : "C2401103"}
                return json.dumps(info)
        if data.group_method not in ["", "sum", "average", "middle"]:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of group_method (%s)" % data.group_method, "code" : "C2401104"}
            return json.dumps(info)
        task_name = 'metagenomic.report.taxon_composition'
        module_type = 'workflow'
        task_type = 2
        anno_info = self.metagenomic.common_find_one("anno_" + data.anno_type,
                                                     {"_id": ObjectId(data.anno_id)})
        name_id = name2id(anno_info["task_id"], type="task")
        name_id = {k: v for k, v in name_id.items() if v in data.group_detail}
        params_json = {
            'anno_type': data.anno_type,
            'anno_id': data.anno_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'group_method': data.group_method,
            'combine_value': float(data.combine_value),
            'graphic_type': data.graphic_type,
            'task_type': task_type,
            'submit_location': data.submit_location
        }
        mongo_data = [
            ('project_sn', anno_info['project_sn']),
            ('task_id', anno_info['task_id']),
            ('status', 'start'),
            ('desc', 'processing'),
            ('anno_type', data.anno_type),
            ('graphic_type', data.graphic_type),
            ('specimen_list', ','.join(name_id.values())),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        if data.graphic_type in ["bubble"]:
            if hasattr(data, "color_level"):
                params_json["color_level"] = int(data.color_level)
        level_lab = ['d', 'k', 'p', 'c', 'o', 'f', 'g', 's']
        col = level_lab[int(data.level_id) - 1] + '__'
        options = {
            'graphic_type': data.graphic_type,
            'table': anno_info["anno_file"],
            'level_id': data.level_id,
            'col': col,
            'group_method': data.group_method,
            'group': data.group_detail,
            'name2id': json.dumps(name_id),
            'others': float(data.combine_value)
        }
        to_file = ['metagenomic.export_group_by_detail_taxon(group)']
        main_table_name = "BarPie_{}_{}_{}".format(data.anno_type.title(), col[0],
                                                   datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.metagenomic.insert_main_table('taxon_composition', mongo_data)
        update_info = {str(main_table_id): 'taxon_composition'}
        options['main_id'] = str(main_table_id)
        options['main_col'] = "taxon_compositon"
        options['update_info'] = json.dumps(update_info)
        # self.set_sheet_data(name=task_name, options=options,
        #                     main_table_name="BarPie_reads/" + main_table_name.strip().split("_")[0] + '/' + main_table_name,
        #                     task_id=anno_info['task_id'], project_sn=anno_info['project_sn'],
        #                     module_type=module_type, params=params_json, to_file=to_file)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name="BarPie_reads/" + main_table_name,
                            task_id=anno_info['task_id'], project_sn=anno_info['project_sn'],
                            module_type=module_type, params=params_json, to_file=to_file)
        task_info = super(TaxonCompositionAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.metagbin_controller import MetagbinController
from mainapp.models.mongo.metagbin import Metagbin
from mainapp.libs.signature import check_sig
from bson import SON


class BinTreeAction(MetagbinController):
    def __init__(self):
        super(BinTreeAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'submit_location', 'task_type', "tree_type", 'method', 'num']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'metagbin.report.bin_tree'
        module_type = 'workflow'
        if int(data.num) <=100:
            if data.tree_type not in ['ML', 'MP', "NJ"]:
                info = {'success': False, 'info': 'parameters missing: tree_type is wrong'}
                return json.dumps(info)
        elif int(data.num) >100:
            if data.tree_type not in ['ML']:
                info = {'success': False, 'info': 'parameters missing: tree_type is not ML'}
                return json.dumps(info)
        if data.method not in ['corgene','16s']:
            info = {'success': False, 'info': 'parameters missing: method is not corgene or 16s'}
            return json.dumps(info)
        else:
            if data.method in ['corgene']:
                if not hasattr(data, "bins"):
                    info = {'success': False, 'info': 'parameters missing: bins is not exsits'}
                    return json.dumps(info)
                if not hasattr(data, "bin_list"):
                    info = {'success': False, 'info': 'parameters missing: bin_list is not exsits'}
                    return json.dumps(info)
            if data.method in ['16s']:
                if hasattr(data, "bins") and hasattr(data, "bin_list"):
                    pass
                elif not hasattr(data, "bins") and not hasattr(data, "bin_list"):
                    pass
                else:
                    info = {'success': False, 'info': 'parameters missing: bins and bin_list must be exsits'}
                    return json.dumps(info)
                if hasattr(data, "genomes") and hasattr(data, "genome_list"):
                    pass
                elif not hasattr(data, "genomes") and not hasattr(data, "genome_list"):
                    pass
                else:
                    info = {'success': False, 'info': 'parameters missing: genomes and genome_list must be exsits'}
                    return json.dumps(info)
                if hasattr(data, "ref") and hasattr(data, "ref_id"):
                    pass
                elif not hasattr(data, "ref") and not hasattr(data, "ref_id"):
                    pass
                else:
                    info = {'success': False, 'info': 'parameters missing: ref and ref_id must be exsits'}
                    return json.dumps(info)
                if not hasattr(data, "bins") and not hasattr(data, "genomes") and not hasattr(data, "ref"):
                    info = {'success': False, 'info': 'parameters missing: bins and genome and ref have one'}
                    return json.dumps(info)
        project_sn = self.metagbin.get_projectsn(data.task_id)
        params = {
            'method': data.method,
            'submit_location': data.submit_location,
            'tree_type': data.tree_type,
            'task_type': int(data.task_type),
            'num': int(data.num)
        }
        main_table_name = 'PlotTree_' + data.method + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
        if data.method in ['16s']:
            if hasattr(data, "ref_id"):
                params['ref_id'] = data.ref_id
            if hasattr(data, "bins"):
                params['bins'] = data.bins
            if hasattr(data, "genomes"):
                params['genomes'] = data.genomes
        elif data.method in ['corgene']:
            params['bins'] = data.bins
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_id = self.metagbin.insert_main_table('tree',mongo_data)
        Metagbin().insert_main_table_new("tree", str(main_id), {"main_id": main_id})
        update_info[str(main_id)] = 'tree'
        options = {
            'update_info': json.dumps(update_info),
            'main_id': str(main_id),
            "num": data.num,
            "method": data.method,
            "tree_type":data.tree_type
            }
        if data.method in ['16s']:
            if hasattr(data, "ref"):
                options['ref'] = data.ref
            if hasattr(data, "bin_list"):
                options['bin_list'] = data.bin_list
                options['bins'] = data.bins
            if hasattr(data, "genome_list"):
                options['genome_list'] = data.genome_list
                options['genomes'] = data.genomes
        elif data.method in ['corgene']:
            options['bin_list'] = data.bin_list
            options['bins'] = data.bins
        self.set_sheet_data(name=task_name,options=options,
                            main_table_name="PlotTree/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(BinTreeAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)
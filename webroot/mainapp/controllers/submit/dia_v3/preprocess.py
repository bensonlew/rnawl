# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.dia_controller import DiaController
from mbio.api.to_file.dia import *
from collections import OrderedDict
import unittest
import os

class PreprocessAction(DiaController):
    def __init__(self):
        super(PreprocessAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['all_eliminate', 'all_percent', 'if_group', 'fill_type',
                        'fillna', 'submit_location', 'task_id', 'task_type']
        task_info = self.dia.get_task_info(data.task_id)

        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "parameters missing:%s" % argu}
                return json.dumps(info)

        if data.fillna not in ["bpca", 'grr', 'impseqrob', 'seqknn', 'impseq', 'rowmedian', 'knnmethod',
                               "min", "mean", "rf", "none"]:
            variables = list()
            variables.append(data.fillna)
            info = {"success": False, "info": "缺失值方法错误：%s" % data.fillna, 'code':'C2301701', 'variables':variables}
            return json.dumps(info)

        if data.if_group == 'yes' and data.group_percent == '':
            info = {"success": False, "info": "请输入阈值标准值"}
            return json.dumps(info)

        task_name = 'dia_v3.report.preprocess'
        module_type = 'workflow'

        params_json = {
            'submit_location': data.submit_location,
            'task_id': data.task_id,
            "task_type": int(data.task_type),
            # 'software': "DIA",
            'all_eliminate': data.all_eliminate,
            'all_percent': data.all_percent,
            'if_group': data.if_group,
            'fill_type': data.fill_type,
            'fillna': data.fillna,
        }

        if hasattr(data, "group_specific"):
            params_json.update({"group_specific": data.group_specific})
        if hasattr(data, "group_percent"):
            params_json.update({"group_percent": data.group_percent})
        # time_now = datetime.datetime.now()
        # name = "ProteinTable_Origin_" + time_now.strftime("%Y%m%d_%H%M%S")
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', data.task_id),
            ('status', 'start'),
            ('desc', 'preprocessed exp profile'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('type', 'origin'),
            ('version', 'v3')
        ]

        main_table_name = "ProteinTable_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))

        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.dia.insert_main_table('sg_express', mongo_data)
        update_info = {str(main_table_id): 'sg_express'}

        options = {
            'update_info': json.dumps(update_info),
            'group_table': data.task_id,
            'fill_type': data.fill_type,
            'fillna': data.fillna,
            'main_table_id': str(main_table_id),
            'name': main_table_name,
            'if_group': data.if_group,
            'all_eliminate': data.all_eliminate,
            'all_percent': float(data.all_percent),
            'searchdb': data.task_id,
        }

        if hasattr(data, "group_specific"):
            options.update({"group_specific": data.group_specific})
        if data.if_group == 'yes' and hasattr(data, "group_percent"):
            options.update({"group_percent": float(data.group_percent)})

        express = self.dia.get_main_info_by_record('sg_express', task_id=data.task_id, type='raw')
        # double check how to find the raw exp table
        if express:
            raw_exp_id = express["main_id"]
        else:
            info = {'success': False, 'info': '该项目没有raw表'}
            return json.dumps(info)
        options['raw_exp_id'] = str(raw_exp_id)
        options['raw_path'] = str(raw_exp_id)

        to_file = ['dia.export_group_origin(group_table)',
                   'dia.export_exp_matrix_all(raw_path)',
                   'dia.export_searchdb(searchdb)']
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name='Preprocess/' + main_table_name,
                            task_id=data.task_id,
                            project_sn=task_info['project_sn'],
                            module_type=module_type, params=params_json, to_file=to_file)

        task_info = super(PreprocessAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }
        }

        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.dia.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.dia.update_group_compare_is_use(data.task_id, data.control_id)

        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/dia_v3/preprocess "
        # cmd += "-b http://192.168.12.102:9090 "
        cmd += '-b http://bcl.tsg.com '
        args = dict(
            task_id="dia_v3",
            task_type="2",
            submit_location="diapreprocess",
            express_id='5fc71df617b2bf494ba82c10',
            all_eliminate='all',
            all_percent='90',
            if_group='yes',
            group_specific='any',
            group_percent='50',
            fill_type='group',
            fillna='min',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

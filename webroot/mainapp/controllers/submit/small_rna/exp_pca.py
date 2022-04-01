# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, qinjincheng'

from mainapp.controllers.project.small_rna_controller import SmallRnaController
from mainapp.libs.signature import check_sig
import web
import json
import datetime
from bson.objectid import ObjectId
from collections import OrderedDict
import unittest
import os

class ExpPcaAction(SmallRnaController):
    def __init__(self):
        super(ExpPcaAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()

        # check arg
        basic_args = ['task_id', 'submit_location', 'task_type']
        basic_args.extend(['group_id', 'group_dict', 'exp_id'])
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Lack argument: {}'.format(arg)}
                return json.dumps(info)

        # set variables
        task_id = data.task_id
        submit_location = data.submit_location
        task_type = data.task_type
        group_id = str(data.group_id)
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        exp_id = data.exp_id

        # prepare main_info for inserting a new document in sg_exp_pca
        exp_info = self.small_rna.get_exp_params_info(exp_id, task_id)
        project_sn = exp_info['project_sn']
        time_now = datetime.datetime.now()
        name = 'ExpPca_{}_{}'.format('all', time_now.strftime("%Y%m%d_%H%M%S"))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        desc = 'Expression pca main table built at controller'
        group_id = data.group_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        eps = 0
        #判断跑不跑置信圈
        for a in group_dict:
            if len(group_dict[a]) <3:
                eps += 1
        # for sample_infos in group_dict.values():
        #     if len(sample_infos) < 3:
        #         eps += 1
            else:
                continue
        if eps == 0:
            if str(data.group_id).lower() != 'all':
                ellipse = 'yes'
            else:
                ellipse = 'no'

        else:
            ellipse = 'no'

        params = {
            'task_id': task_id,
            'submit_location': submit_location,
            'task_type': int(task_type),
            'group_id': group_id,
            'group_dict': group_dict,
            'exp_id': exp_id,
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = {
            # essential keys
            'task_id': task_id,
            'project_sn': project_sn,
            'created_ts': created_ts,
            'name': name,
            'desc': desc,
            'params': params,
            # alterable keys
            'exp_id': ObjectId(exp_id),
            # status of process
            'status': 'start'
        }

        if ellipse == 'yes':
            main_info.update({
                'ellipse': 'yes'
            })
        if ellipse == 'no':
            main_info.update({
                'ellipse': 'no'
            })

        main_id = self.small_rna.insert_main_table('sg_exp_pca', main_info)

        # prepare options for workflow
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        options = {
            'ellipse': ellipse,
            'exp_matrix': exp_id, # to_file require
            'group_dict': json.dumps(group_dict), # to_file require, no use at workflow
            "pca_main_id": str(main_id),
            "group_table": json.dumps(group_dict),
            'update_info': json.dumps({str(main_id): 'sg_exp_pca'})
        }
        # prepare to file
        to_files = ['small_rna.export_exp_matrix(exp_matrix)',
                    'small_rna.export_group(group_table)']
        # prepare sheet data for workflow
        task_name = 'small_rna.report.exp_pca'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,
                            module_type='workflow',
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=task_id)

        # run workflow and obtain return value
        task_info = super(ExpPcaAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': name}}
        task_info['group_dict'] = group_dict

        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.small_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.small_rna.update_group_compare_is_use(data.task_id, data.control_id)

        return json.dumps(task_info)

class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''

    def test(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += '-fr no '
        cmd += '-c {} '.format('client03')
        cmd += 's/small_rna/exp_pca '
        cmd += '-b http://192.168.12.102:9090 '
        args = {
            'task_id': 'small_rna',
            'submit_location': 'exp_pca',
            'task_type': '2',

            'group_id': '5be24abac82c5131355dcd00',
            'group_dict': json.dumps({'GH':['GH1','GH2'], 'NFAs':['NFAs1','NFAs2'], 'Normal':['Normal1','Normal2'], 'PRL': ['PRL1', 'PRL2']}).replace('"', '\\"'),
            'exp_id': '5bfb8c0ea4e1af5ee61ea9b6',
        }
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()

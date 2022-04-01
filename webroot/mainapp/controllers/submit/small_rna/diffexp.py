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

class DiffexpAction(SmallRnaController):
    def __init__(self):
        super(DiffexpAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()

        # check arg
        basic_args = ['task_id', 'submit_location', 'task_type']
        basic_args.extend(['group_id', 'group_dict', 'control_id', 'exp_id'])
        # diff_method in ['DESeq2', 'edegR', 'DEGseq']
        # stat_type in ['padjust', 'pvalue']
        # stat_cutoff ~ 0.05
        # float(fc) > 1.0
        basic_args.extend(['diff_method', 'fc'])
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Lack argument: {}'.format(arg)}
                return json.dumps(info)

        # set variables
        task_id = data.task_id
        submit_location = data.submit_location
        task_type = int(data.task_type)
        group_id = str(data.group_id)
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        control_id = data.control_id
        exp_id = data.exp_id
        diff_method = data.diff_method
        # stat_type = data.stat_type
        # stat_cutoff = data.stat_cutoff
        # fc = data.fc
        seq_type = 'all'

        # check diff_method
        exp_info = self.small_rna.get_exp_params_info(exp_id, task_id)
        project_sn = exp_info['project_sn']
        if group_id.lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        group_size = list()
        message = list()
        for key in group_dict:
            group_size.append(len(group_dict[key]))
        group_size.sort()
        if group_size[0] == group_size[-1] == 1:
            if diff_method == 'DESeq2':
                # message.append('DESeq2 can not handle single sample comparison, please select DEGseq or edgeR')
                message.append('请选择DEGseq、edgeR、NOIseq或者Limma进行单样本组与单样本组的组间比较')
        elif group_size[0] == 1 and group_size[1] >= 2:
            if diff_method == 'DESeq2' or diff_method == 'DEGseq':
                # message.append('if you want to make a comparison betwenn single sample and multiple samples, please select edgeR')
                message.append('请选择edgeR进行单样本组与多样本组的组间比较')
        # elif group_size[0] >= 2:
        #     if diff_method == 'DEGseq':
        #         # message.append('if you want to make a comparison betwenn multiple samples and multiple samples, we recommand DESeq2 or edgeR')
        #         message.append('请选择DESeq2或edgeR进行多样本组与多样本组的组间比较')
        if message:
            info = {'success': False, 'info': message[0]}
            return json.dumps(info)

        # prepare main_info for inserting a new document in sg_diff
        params = {
            'task_id': task_id,
            'submit_location': submit_location,
            'task_type': int(task_type),
            'group_id': group_id,
            'group_dict': group_dict,
            'control_id': control_id,
            'exp_id': exp_id,
            'diff_method': diff_method,
            # 'stat_type': stat_type,
            # 'stat_cutoff': stat_cutoff,
            'fc': str(float(data.fc)),
        }
        # for special args
        if hasattr(data, 'correct_method'):
            correct_method = data.correct_method
            params.update(dict(correct_method=correct_method))
        else:
            correct_method = 'BH'
        if hasattr(data, 'stat_type'):
            stat_type = data.stat_type
            params.update(dict(stat_type=stat_type))
        if hasattr(data, 'stat_cutoff'):
            stat_cutoff = data.stat_cutoff
            params.update(dict(stat_cutoff=stat_cutoff))
        if hasattr(data, "prob"):
            params.update({
                "prob": data.prob
            })
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))

        # create main table record
        time_now = datetime.datetime.now()
        name = 'Diff_{}_{}'.format(diff_method, time_now.strftime("%Y%m%d_%H%M%S"))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        desc = 'Differential expression main table built at controller'
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
            # 'group_id': ObjectId(group_id),
            # 'group_dict': group_dict,
            # 'control_id': ObjectId(control_id),
            # 'diff_method': str(diff_method),
            # 'stat_type': str(stat_type),
            # 'stat_cutoff': float(stat_cutoff),
            # 'fc': float(fc),
            # status of process
            'status': 'start'
        }
        main_id = self.small_rna.insert_main_table('sg_diff', main_info)

        # prepare options for workflow
        options = {
            'exp_matrix': exp_id, # to_file require
            'group_dict': json.dumps(group_dict), # to_file require, no use at workflow
            'seq_type': seq_type,  # to_file require, no use at workflow
            'count': self.small_rna.get_count_exp_id(exp_id),
            'exp_id': str(exp_id),
            'group': json.dumps(group_dict),
            'cmp': control_id,
            'diff_main_id': str(main_id),
            'method': diff_method,
            'fc': float(data.fc),
            # 'pvalue': float(data.stat_cutoff),
            # 'pvalue_padjust': stat_type,
            # 'padjust_way': correct_method,
            'update_info': json.dumps({str(main_id): "sg_diff"})
        }
        if hasattr(data, 'stat_type'):
            options.update({
                'pvalue_padjust': data.stat_type
            })
        if hasattr(data, 'stat_cutoff'):
            options.update({
                'pvalue':float(data.stat_cutoff)
            })
        if hasattr(data, 'correct_method'):
            options.update({
                'padjust_way': data.correct_method
            })
        if hasattr(data, 'prob'):
            options.update({
                'prob': data.prob
            })

        # prepare to file
        to_files = ['small_rna.export_exp_matrix(exp_matrix)',
                    'small_rna.export_exp_matrix(count)',
                    'small_rna.export_group(group)',
                    'small_rna.export_compare(cmp)']

        # prepare sheet data for workflow
        task_name = 'small_rna_v2.report.diffexp'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,
                            module_type='workflow',
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=task_id)

        # run workflow and obtain return value
        task_info = super(DiffexpAction, self).POST()
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

    def test_DESeq2_padjust_BH(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += '-fr no '
        cmd += '-c {} '.format('client03')
        cmd += 's/small_rna/diffexp '
        cmd += '-b http://bcl.tsg.com '
        args = {
            'task_id': 'tsg_37304',
            'submit_location': 'diffexp',
            'task_type': '2',

            'group_id': '5eba05b317b2bf24ccac5ff1',
            'group_dict': json.dumps({"ABA":["ABA_1","ABA_2","ABA_3"],"BR":["BR_1","BR_2","BR_3"],"GA3":["GA3_1","GA3_2","GA3_3"],"IAA":["IAA_1","IAA_2","IAA_3"],"Ler":["Ler_1","Ler_2","Ler_3"],"NDGA":["NDGA_1","NDGA_2","NDGA_3"],"TIBA":["TIBA_1","TIBA_2","TIBA_3"],"ag11":["ag11_1","ag11_2","ag11_3"]}).replace('"', '\\"'),
            'control_id': '5eba05b317b2bf24ccac5ff2',
            'exp_id': '5eba071517b2bf24ccaf063a',

            'diff_method': 'DESeq2',
            'stat_type': 'padjust',
            'stat_cutoff': '0.05',
            'correct_method': 'BH',
            # 'prob': '0.8',
            'fc': '2.0',
        }
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print cmd
        os.system(cmd)

    # def test_DEGseq_padjust_BH(self):
    #     cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
    #     cmd += 'post '
    #     cmd += '-fr no '
    #     cmd += '-c {} '.format('client03')
    #     cmd += 's/small_rna/diffexp '
    #     cmd += '-b http://192.168.12.102:9090 '
    #     args = {
    #         'task_id': 'small_rna',
    #         'submit_location': 'diffexp',
    #         'task_type': '2',
    #
    #         'group_id': '5be24abac82c5131355dcd00',
    #         'group_dict': json.dumps({'GH':['GH1'], 'NFAs':['NFAs1'], 'Normal':['Normal1'], 'PRL': ['PRL1']}).replace('"', '\\"'),
    #         'control_id': '5be24c18c82c5131355ddc5f',
    #         'exp_id': '5bfb8c0ea4e1af5ee61ea9b6',
    #
    #         'diff_method': 'DEGseq',
    #         'stat_type': 'padjust',
    #         'stat_cutoff': '0.05',
    #         'correct_method': 'BH',
    #         'fc': '2.0',
    #     }
    #     arg_names, arg_values = args.keys(), args.values()
    #     cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
    #     print cmd
    #     os.system(cmd)
    #
    # def test_edgeR_padjust_BH(self):
    #     cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
    #     cmd += 'post '
    #     cmd += '-fr no '
    #     cmd += '-c {} '.format('client03')
    #     cmd += 's/small_rna/diffexp '
    #     cmd += '-b http://192.168.12.102:9090 '
    #     args = {
    #         'task_id': 'small_rna',
    #         'submit_location': 'diffexp',
    #         'task_type': '2',
    #
    #         'group_id': '5be24abac82c5131355dcd00',
    #         'group_dict': json.dumps({'GH':['GH1','GH2'], 'NFAs':['NFAs1','NFAs2'], 'Normal':['Normal1','Normal2'], 'PRL': ['PRL1', 'PRL2']}).replace('"', '\\"'),
    #         'control_id': '5be24c18c82c5131355ddc5f',
    #         'exp_id': '5bfb8c0ea4e1af5ee61ea9b6',
    #
    #         'diff_method': 'edgeR',
    #         'stat_type': 'padjust',
    #         'stat_cutoff': '0.05',
    #         'correct_method': 'BH',
    #         'fc': '2.0',
    #     }
    #     arg_names, arg_values = args.keys(), args.values()
    #     cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
    #     print cmd
    #     os.system(cmd)

if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os
from collections import OrderedDict
from bson.objectid import ObjectId


class ExpDistributionAction(ProkRNAController):
    def __init__(self):
        super(ExpDistributionAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'exp_id', 'group_dict', 'exp_level']
        # basic_args += ['exp_level']  # exp_level分为三种mrna, srna, mrnasrna
        # v3.2 exp_level 更新为ref_gene, novel_gene, sRNA
        # check arg
        '''
                Array
        (
            [exp_id] => 5b4fe64cf6b9e46b1e6ddc50
            [exp_level] => G
            [group_dict] => {"A1":["A1_1","A1_2","A1_3"],"A2":["A2_1","A2_2","A2_3"]}
            [group_id] => 5b4dabdfffec6004d1584a10
            [submit_location] => expgraph
            [task_id] => tsanger_31035
            [task_type] => 2
            [type] => all
        )
        '''

        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        exp_info = self.prok_rna.get_exp_params_info(data.exp_id, data.task_id)  # 返回sg_exp这次的分析结果
        project_sn = exp_info["project_sn"]
        task_id = data.task_id

        if data.group_id.lower() == 'all':
            group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        else:
            collection = self.db['sg_specimen_group']
            result = collection.find_one({"_id": ObjectId(data.group_id)})
            group_names = result["category_names"]
            group_dict_in = json.loads(data.group_dict)
            group_dict = OrderedDict()
            for each in group_names:
                if each in group_dict_in.keys():
                    group_dict[each] = group_dict_in[each]

        # create main table record
        exp_level = data.exp_level
        quant_method = exp_info['method']
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            group_id=data.group_id,
            group_dict=group_dict,
            exp_level=data.exp_level,
            # type=data.type,
        )
        name = "ExpDistribution" + '_' + data.exp_level + '_' + quant_method + '_'
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='density and box and volin plot main table',
            params=params,
            status="start",
            exp_level=exp_level,
        )
        main_id = self.prok_rna.insert_main_table('sg_exp_graph', main_info)

        # prepare option for workflow
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        options = {
            "exp_matrix": data.exp_id,
            "group_dict": json.dumps(group_dict),
            "graph_main_id": str(main_id),
            "exp_level": data.exp_level,
            # "type": data.type,  # 这个type和最后的道标会产生关系，需要认真考虑
            "update_info": json.dumps({str(main_id): "sg_exp_graph"})  # to update sg_status
        }

        # prepare to file
        to_files = ["prok_rna.export_exp_matrix_prok(exp_matrix)",]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'prok_rna.report.exp_distribution'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ExpDistributionAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.prok_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.prok_rna.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/prok_rna/exp_distribution "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="prok_rna_srna",
            task_type="2",
            submit_location="expgraph",
            exp_id="5b877be2a4e1af3eb5578ab9",
            group_id="5b72454477b3f3b113361971",
            exp_level='mRNA',
            group_dict=r'{"WT":["WT_1", "WT_2", "WT_3"],"rcsBKO": [ "rcsBKO_1", "rcsBKO_2", "rcsBKO_3"]}'.replace('"', '\\"'),
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)

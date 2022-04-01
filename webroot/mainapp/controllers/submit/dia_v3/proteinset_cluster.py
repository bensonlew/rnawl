# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.dia_controller import DiaController
from mbio.api.to_file.dia import *
from mainapp.libs.signature import check_sig
import unittest
import os


class ProteinsetClusterAction(DiaController):
    def __init__(self):
        super(ProteinsetClusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'use_group',
                       'proteinset_id', 'sct', 'gct', ]
        basic_args += ["express_id"]
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        # if data.exp_level != data.proteinset_type:
        #     info = "exp_level is not consistent with proteinset_type!"
        #     return json.dumps(info)
        if data.group_id.lower() == 'all' and not data.use_group.lower() in ['none', 'no']:
            info = {'success': False, 'info': "group 设置为all时,不能选择以分组均值分析"}
            return json.dumps(info)
        # exp_info = self.dia.get_exp_params_info(data.task_id)
        exp_info = self.dia.get_main_info(main_id=data.express_id, collection_name="sg_express", task_id=data.task_id)
        # express_id = str(exp_info['main_id'])
        print(data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)

        # create main table record
        # exp_level = exp_info['exp_level']
        # exp_type, quant_method = exp_info['exp_type'], exp_info['method']
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            express_id=data.express_id,
            group_id=data.group_id,
            group_dict=group_dict,
            # quant_method=quant_method,
            use_group=data.use_group,
            proteinset_id=data.proteinset_id,
            sct=data.sct,
            gct=data.gct,
        )

        n_clusters = 10
        scm = "complete"
        scd = "correlation"
        gcm = "average"
        gcd = "euclidean"
        if data.sct == 'no':
            pass
        elif data.sct == 'kmeans':
            params.update(dict(scd=data.scd, n_clusters=int(data.n_clusters), ))
            scd = data.scd
            n_clusters = int(data.n_clusters)
        else:
            params.update(dict(scm=data.scm, scd=data.scd))
            scm = data.scm
            scd = data.scd
        if data.gct.startswith('h'):
            params.update(dict(gcm=data.gcm, gcd=data.gcd, n_clusters=int(data.n_clusters)))
            gcm = data.gcm
            gcd = data.gcd
            n_clusters = int(data.n_clusters)
        elif data.gct.startswith('k'):
            params.update(dict(gcd=data.gcd, n_clusters=int(data.n_clusters)))
            gcd = data.gcd
            n_clusters = int(data.n_clusters)
        else:
            pass

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Pset_Cluster" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            group_dict=json.loads(data.group_dict, object_pairs_hook=OrderedDict),
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='proteinset cluster main table',
            params=params,
            status="start",
            version='v3',
        )
        main_id = self.dia.insert_main_table('sg_proteinset_cluster', main_info)

        print("++++ exp_id protein_set_id")
        print(data.express_id + ";" + data.proteinset_id)
        # prepare option for workflow
        options = dict(
            exp_matrix=data.express_id + ";" + data.proteinset_id,
            group_dict=data.group_dict,
            group_id=data.group_id,
            group=data.group_dict,
            cluster_main_id=str(main_id),
            n_clusters=n_clusters,
            use_group=data.use_group,
            scm=scm,
            scd=scd,
            sct=data.sct,
            gcm=gcm,
            gct=data.gct,
            gcd=gcd,
            update_info=json.dumps({str(main_id): "sg_proteinset_cluster"})  # to update sg_status
        )

        # prepare to file
        to_files = ["dia.export_proteinset_exp_matrix_new(exp_matrix)",
                    "dia.export_group(group)", ]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'dia_v3.report.proteinset_cluster'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ProteinsetClusterAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # 更新基因集的使用信息

        self.dia.insert_proteinset_info(data.proteinset_id, name, str(main_id))
        ##
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.dia.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.dia.update_group_compare_is_use(data.task_id, data.control_id)
        ##
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
        cmd += "s/dia_v3/proteinset_cluster "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="dia_test",
            task_type="2",
            submit_location="proteinsetcluster",
            group_id="5fc71d5217b2bf494ba7ddef",
            group_dict=json.dumps({"C_3": ["C_3_1", "C_3_2", "C_3_3", "C_3_4", "C_3_5"],
                                   "R_3": ["R_3_1", "R_3_2", "R_3_3", "R_3_4", "R_3_5"]}).replace('"', '\\"'),
            express_id="5fc737f517b2bf708ed8252e",
            proteinset_id="5fc5f40c17b2bf1d5aaa2c28",
            use_group='no',
            sct="hierarchy",
            gct="hierarchy",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mbio.api.to_file.denovo_rna_v2 import *
from mainapp.libs.signature import check_sig
import unittest
import os


class GenesetClusterAction(DenovoRnaV2Controller):
    def __init__(self):
        super(GenesetClusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'exp_id', 'exp_level', 'use_group',
                       'geneset_id', 'geneset_type', 'sct', 'gct',
                       # 'scm', 'gcm',
                       # 'scd', 'gcd',
                       # 'n_clusters',
                       ]
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                var = []
                var.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % (arg), "code": 'C1601001', "variables": var}
                return json.dumps(info)
        if data.exp_level != data.geneset_type:
            info = {'success': False, 'info': "exp_level is not consistent with geneset_type!", "code": 'C1601002', "variables": ''}
            # info = "exp_level is not consistent with geneset_type!"
            return json.dumps(info)
        if data.group_id.lower() == 'all' and data.use_group.lower() == 'yes':
            info = {'success': False, 'info': "group 设置为all时,不能选择以分组均值分析", "code": 'C1601003', "variables": ''}
            return json.dumps(info)
        exp_info = self.denovo_rna_v2.get_exp_params_info(data.exp_id, data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)

        # create main table record
        exp_level = exp_info['exp_level']
        # exp_type, quant_method = exp_info['exp_type'], exp_info['method']
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            group_id=data.group_id,
            exp_level=exp_level,
            group_dict=group_dict,
            # quant_method=quant_method,
            use_group=data.use_group,
            geneset_id=data.geneset_id,
            geneset_type=data.geneset_type,
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
        name = "GeneSet_Cluster" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            version="v2",
            task_id=task_id,
            name=name,
            group_dict=json.loads(data.group_dict, object_pairs_hook=OrderedDict),
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='geneset cluster main table',
            params=params,
            status="start",
            type=exp_level,
        )
        main_id = self.denovo_rna_v2.insert_main_table('sg_geneset_cluster', main_info)

        # prepare option for workflow
        options = dict(
            exp_matrix=data.exp_id + ";" + data.geneset_id,
            group_dict=data.group_dict,
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
            update_info=json.dumps({str(main_id): "sg_geneset_cluster"})  # to update sg_status
        )

        # prepare to file
        to_files = ["denovo_rna_v2.export_geneset_exp_matrix(exp_matrix)",
                    "denovo_rna_v2.export_group(group)",]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'denovo_rna_v2.report.geneset_cluster'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(GenesetClusterAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # 更新基因集的使用信息
        self.denovo_rna_v2.insert_geneset_info(data.geneset_id, 'sg_geneset_cluster', str(main_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.denovo_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.denovo_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/denovo_rna_v2/geneset_cluster "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="denovo_rna_v2_upgrade",
            task_type="2",
            submit_location="genesetcluster",
            group_id="5d47f30b17b2bf08c619bf78",
            group_dict=json.dumps({"A":["A1","A2","A3"],"B":["B1","B2","B3"],"C":["C1","C2","C3"]}).replace('"', '\\"'),
            exp_id="5d41458917b2bf10c8b3ffb0",
            exp_level="G",
            geneset_id="5cb64c0d17b2bf6187b11fbd",
            geneset_type="G",
            quant_method='RSEM',
            use_group='no',
            n_clusters='10',
            sct="hierarchy",
            scm="complete",
            scd="euclidean",
            gct="no",
            gcm="average",
            gcd="euclidean",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

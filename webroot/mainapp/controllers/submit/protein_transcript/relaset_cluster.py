# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.protein_transcript_controller import ProteinTranscriptController
from mbio.api.to_file.itraq_tmt import *
from mainapp.libs.signature import check_sig
import unittest
import os


class RelasetClusterAction(ProteinTranscriptController):
    def __init__(self):
        super(RelasetClusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "protein_group_dict", 'rna_group_dict', 'rna_exp_id']
        basic_args += ['relaset_id', 'protein_group_id', 'rna_group_id']
        basic_args += ['gct', 'use_group']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        if not data.protein_group_dict or not data.rna_group_dict:
            info = {'success': False, 'info': "样品分组为空"}
            return json.dumps(info)
        relation_info = self.itraq_tmt.get_relation_condition(data.task_id)
        if not relation_info:
            info = {'success': False, 'info': "关联信息为空，请确保改项目已经跟转录组项目做好关联"}
            return json.dumps(info)

        if (data.rna_group_id.lower() == 'all' or data.protein_group_id.lower() == 'all') and data.use_group.lower() != 'no':
            info = {'success': False, 'info': "group 设置为all时,不能选择以分组均值分析"}
            return json.dumps(info)

        project_sn = relation_info["project_sn"]
        rna_type = relation_info['rna_type']
        rna_task_id = relation_info['rna_task_id']

        protein_group_dict = json.loads(data.protein_group_dict, object_pairs_hook=OrderedDict)
        rna_group_dict = json.loads(data.rna_group_dict, object_pairs_hook=OrderedDict)
        if (len(protein_group_dict.keys()) < 2 or len(rna_group_dict.keys()) < 2) and data.use_group.lower() != "no":
            info = {'success': False, 'info': "使用分组样本计算(求中值或中位数),需选择两个或两个以上的组别"}
            return json.dumps(info)

        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            rna_exp_id=data.rna_exp_id,
            # protein_group_dict=data.protein_group_dict,
            protein_group_dict=protein_group_dict,
            # rna_group_dict=data.rna_group_dict,
            rna_group_dict=rna_group_dict,
            protein_group_id=data.protein_group_id,
            rna_group_id=data.rna_group_id,
            gct=data.gct,
            use_group=data.use_group,
            relaset_id=data.relaset_id,
        )

        n_clusters = 10
        gcm = "average"
        gcd = "euclidean"
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
        name = "relaset_cluster" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            rna_type=rna_type,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='the cluster analyse of one relaset',
            params=params,
            status="start",
            version='v2.2'
        )
        main_id = self.itraq_tmt.insert_main_table('sg_relaset_cluster', main_info)

        # prepare option for workflow
        options = dict(
            protein_matrix=data.task_id,
            rna_matrix=str(data.rna_exp_id),
            rna_group_dict=data.rna_group_dict,
            group=data.rna_group_dict,
            rna_group_id=data.rna_group_id,
            protein_group_dict=data.protein_group_dict,
            protein_group_id=data.protein_group_id,
            relaset_list=data.relaset_id,
            rna_type=rna_type,
            n_clusters=n_clusters,
            use_group=data.use_group,
            gcm=gcm,
            gct=data.gct,
            gcd=gcd,
            main_id=str(main_id),
            update_info=json.dumps({str(main_id): "sg_relaset_cluster"})  # to update sg_status
        )

        # prepare to file
        to_files = ["protein_transcript.export_exp_matrix_itraq(protein_matrix)",
                    "protein_transcript.export_relaset_list(relaset_list)",
                    "protein_transcript.export_rna_protein_group(group)",
                    ]
        if rna_type == 'prok_rna':
            to_files.append("protein_transcript.export_exp_matrix_prok(rna_matrix)")
        if rna_type == 'ref_rna_v1':
            to_files.append("protein_transcript.export_exp_matrix_refrnav1(rna_matrix)")
        if rna_type == 'ref_rna_v2':
            to_files.append("protein_transcript.export_exp_matrix_refrnav2(rna_matrix)")
        if rna_type == 'denovo_rna_v2':
            to_files.append("protein_transcript.export_exp_matrix_denovo(rna_matrix)")
        if rna_type == 'whole_transcriptome':
            to_files.append("protein_transcript.export_exp_matrix_whole_transcriptome(rna_matrix)")
        if rna_type == 'medical_transcriptome':
            to_files.append("protein_transcript.export_exp_matrix_medical_transcriptome(rna_matrix)")
        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'protein_transcript.report.relaset_cluster'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(RelasetClusterAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # 更新基因集的使用信息
        self.itraq_tmt.insert_relaset_info(data.relaset_id, name, str(main_id))
        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/isanger/sanger_bioinfo/bin/webapitest_new.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/protein_transcript/relaset_cluster "
        cmd += "-b http://bcl.i-sanger.com "
        cmd += "-dbversion 1 "
        args = dict(
            task_id="majorbio_318246",
            task_type="2",
            submit_location="relasetcluster",
            rna_exp_id="5ff775e0ec02cc0b1b91af05",
            rna_group_dict=json.dumps({"M":["M1","M2","M3"],"W":["W1","W2","W3"]}).replace('"', '\\"'),
            protein_group_dict=json.dumps({"Wild_type":["W1","W2","W3"],"Mutant":["M1","M2","M3"]}).replace('"', '\\"'),
            rna_group_id='5ff77534ec02cc0b1b6c7438',
            protein_group_id='6010d2aaf6b9e46cb0837349',
            use_group='no',
            n_clusters='10',
            gct="hierarchy",
            gcm="average",
            gcd="euclidean",
            relaset_id="6011331759c37d855d8b458f"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.protein_transcript_dia_v3_controller import ProteinTranscriptDiaV3Controller
from mbio.api.to_file.dia import *
from mainapp.libs.signature import check_sig
import unittest
import os


class ExpressPrCorrAction(ProteinTranscriptDiaV3Controller):
    def __init__(self):
        super(ExpressPrCorrAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "protein_group_dict", 'rna_group_dict', 'rna_exp_id']
        basic_args += ['corr_method', 'geneset_id', 'proteinset_id']
        basic_args += ['cor_cutoff', 'padjust_way', 'sig_type']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        if not data.protein_group_dict or not data.rna_group_dict:
            info = {'success': False, 'info': "样品分组为空"}
            return json.dumps(info)
        relation_info = self.dia.get_relation_condition(data.task_id)
        if not relation_info:
            info = {'success': False, 'info': "关联信息为空，请确保改项目已经跟转录组项目做好关联"}
            return json.dumps(info)
        if not hasattr(data, 'corr_method') or data.corr_method.lower() not in ['spearman', 'pearson', 'kendall']:
            info = {'success': False, 'info': "没有您输入的相关性方法"}
            return json.dumps(info)

        project_sn = relation_info["project_sn"]
        rna_type = relation_info['rna_type']
        rna_task_id = relation_info['rna_task_id']

        protein_group_dict = json.loads(data.protein_group_dict, object_pairs_hook=OrderedDict)
        rna_group_dict = json.loads(data.rna_group_dict, object_pairs_hook=OrderedDict)
        if len(protein_group_dict) != len(rna_group_dict):
            info = {'success': False, 'info': "蛋白和转录的组别不一样，无法做相关性分析"}
            return json.dumps(info)
        for n, samples in enumerate(protein_group_dict.values()):
            if len(samples) != len(rna_group_dict.values()[n]):
                info = {'success': False, 'info': "蛋白和转录某组别的样本数量不一样，无法做相关性分析"}
                return json.dumps(info)

        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            protein_group_dict=protein_group_dict,
            rna_group_dict=rna_group_dict,
            protein_group_id=data.protein_group_id,
            rna_group_id=data.rna_group_id,
            corr_method=data.corr_method,
            geneset_id=data.geneset_id,
            proteinset_id=data.proteinset_id,
            cor_cutoff=data.cor_cutoff,
            padjust_way=data.padjust_way,
            sig_type=data.sig_type,
            rna_exp_id=data.rna_exp_id
        )

        if "qvalue" in data:
            params.update({'qvalue': data.qvalue})
        if "pvalue" in data:
            params.update({'pvalue': data.pvalue})
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "express_pr_corr" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            rna_type=rna_type,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='the corr analyse of the express of proteinset and geneset',
            params=params,
            status="start"
        )
        main_id = self.dia.insert_main_table('sg_express_pr_corr', main_info)

        # prepare option for workflow
        options = dict(
            protein_matrix=data.task_id,
            rna_matrix=str(data.rna_exp_id),
            rna_group_dict=data.rna_group_dict,
            protein_group_dict=data.protein_group_dict,
            geneset_list=str(data.geneset_id) + ',' + rna_type,
            proteinset_list=str(data.proteinset_id),
            group_rna=data.rna_group_dict,
            group_protein=data.protein_group_dict,
            corr_method=data.corr_method,
            rna_type=rna_type,
            cor_cutoff=data.cor_cutoff,
            padjust_way=data.padjust_way,
            main_id=str(main_id),
            update_info=json.dumps({str(main_id): "sg_express_pr_corr"})  # to update sg_status
        )

        # if "qvalue" in data:
        if data.sig_type == 'qvalue':  # modified by zhangyitong on 20211118
            options.update({'sig_type':1, 'qvalue_cutoff': data.qvalue})
        # if "pvalue" in data:
        elif data.sig_type == "pvalue":
            options.update({'sig_type':0, 'pvalue_cutoff': data.pvalue})

        # prepare to file
        to_files = ["protein_transcript_dia_v3.export_exp_matrix_labelfree(protein_matrix)",
                    "protein_transcript_dia_v3.export_proteinset_list(proteinset_list)",
                    "protein_transcript_dia_v3.export_geneset_list(geneset_list)",
                    "protein_transcript_dia_v3.export_rna_group(group_rna)",
                    "protein_transcript_dia_v3.export_protein_group(group_protein)",
                    ]
        if rna_type == 'prok_rna':
            to_files.append("protein_transcript_dia_v3.export_exp_matrix_prok(rna_matrix)")
        if rna_type == 'ref_rna_v1':
            to_files.append("protein_transcript_dia_v3.export_exp_matrix_refrnav1(rna_matrix)")
        if rna_type == 'ref_rna_v2':
            to_files.append("protein_transcript_dia_v3.export_exp_matrix_refrnav2(rna_matrix)")
        if rna_type == 'denovo_rna_v2':
            to_files.append("protein_transcript_dia_v3.export_exp_matrix_denovo(rna_matrix)")
        if rna_type == 'whole_transcriptome':
            to_files.append("protein_transcript_dia_v3.export_exp_matrix_whole_transcriptome(rna_matrix)")
        if rna_type == 'medical_transcriptome':
            to_files.append("protein_transcript_dia_v3.export_exp_matrix_medical_transcriptome(rna_matrix)")
        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'protein_transcript_dia_v3.report.express_pr_corr'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ExpressPrCorrAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # 更新基因集的使用信息
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
        cmd += "s/protein_transcript_labelfree/express_pr_corr "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_34739",
            task_type="2",
            submit_location="expressprcorr",
            rna_exp_id="5c218d70a4e1af619b59ed03",
            rna_group_dict=json.dumps({"Con":["Con1","Con2","Con3"],"Drug":["Drug1","Drug2","Drug3"]}).replace('"', '\\"'),
            protein_group_dict=r'{"F":["F_1","F_2","F_3"],"S":["S_1","S_2","S_3"]}'.replace(
                '"', '\\"'),
            corr_method="spearman",
            proteinset_id="5d22a83b17b2bf1c689e2823",
            # proteinset_id="5c2103e0a4e1af61602f1f91",
            geneset_id="5c218d9fa4e1af619b5b6b3c",
            # geneset_id="5c218d9fa4e1af619b5b6b3c",
            pvalue="0.01",
            qvalue="0.01",
            cor_cutoff="0.8",
            sig_type="qvalue",
            padjust_way="fdr_bh",
            protein_group_id="5d22a82717b2bf1c689db738",
            rna_group_id="5c218916a4e1af619b2de3c4"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()

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


class PrKeggClassAction(ProteinTranscriptDiaV3Controller):
    def __init__(self):
        super(PrKeggClassAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "protein_kegg_id", 'rna_kegg_id', ]
        basic_args += ['geneset_id', 'proteinset_id']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        relation_info = self.dia.get_relation_condition(data.task_id)
        task_info = self.dia.get_task_info(data['task_id'])
        if not relation_info:
            info = {'success': False, 'info': "关联信息为空，请确保改项目已经跟转录组项目做好关联"}
            return json.dumps(info)

        project_sn = relation_info["project_sn"]
        rna_type = relation_info['rna_type']
        rna_task_id = relation_info['rna_task_id']

        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            protein_kegg_id=str(data.protein_kegg_id),
            rna_kegg_id=str(data.rna_kegg_id),
            geneset_id=str(data.geneset_id),
            proteinset_id=str(data.proteinset_id),
        )

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "relate_kegg_class" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            rna_type=rna_type,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='the kegg_anno_results class of proteinsets and genesets',
            params=params,
            status="start"
        )
        main_id = self.dia.insert_main_table('sg_pr_kegg_class', main_info)

        # prepare option for workflow
        options = dict(
            gene_kegg_table=rna_type + ',' + str(data.rna_kegg_id),
            gene_kegg_level_table=rna_type + ',' + str(data.rna_kegg_id),
            protein_kegg_table=str(data.protein_kegg_id),
            protein_kegg_level_table=str(data.protein_kegg_id),
            protein_kegg_id=str(data.protein_kegg_id),
            pr_list=rna_type + '___' + str(data.geneset_id) + '___' + str(data.proteinset_id),
            add_info_rna=rna_type + ',' + str(data.rna_kegg_id),
            add_info_protein=str(data.protein_kegg_id),
            task_id=str(data.task_id),
            rna_type=rna_type,
            main_id=str(main_id),
            update_info=json.dumps({str(main_id): "sg_pr_kegg_class"})  # to update sg_status
        )
        if "database_version" in  task_info:
            kegg_version = task_info["database_version"].get("kegg", "2017").split("_")[0]
            options.update({"kegg_version": kegg_version})

        # prepare to file
        to_files = ["protein_transcript_dia_v3.export_add_info_protein(add_info_protein)",
                    "protein_transcript_dia_v3.export_add_info_rna(add_info_rna)",
                    "protein_transcript_dia_v3.export_kegg_table_protein(protein_kegg_table)",
                    "protein_transcript_dia_v3.export_kegg_table_rna(gene_kegg_table)",
                    "protein_transcript_dia_v3.export_kegg_level_table_protein(protein_kegg_level_table)",
                    "protein_transcript_dia_v3.export_kegg_level_table_rna(gene_kegg_level_table)",
                    "protein_transcript_dia_v3.export_genepro_list(pr_list)",
                    ]
        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'protein_transcript_dia_v3.report.pr_kegg_class'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(PrKeggClassAction, self).POST()
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
        cmd += "s/protein_transcript_labelfree/pr_kegg_class "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_34739",
            # task_id="tsg_33293",
            task_type="2",
            submit_location="prkeggclass",
            # protein_kegg_id="5c2103d5a4e1af61602e9c67",
            protein_kegg_id="5d22a83917b2bf1c689e1264",
            rna_kegg_id="5c218caaa4e1af619b4fd964",
            # proteinset_id="5c2103e0a4e1af61602f1f93",
            proteinset_id="5d22a83b17b2bf1c689e2825,5d22a83b17b2bf1c689e2827",
            # geneset_id="5c218da2a4e1af619b5b6b3e",
            geneset_id="5c218d9fa4e1af619b5b6b3c,5c218da2a4e1af619b5b6b3e",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()

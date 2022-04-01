# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.protein_transcript_labelfree_controller import ProteinTranscriptLabelfreeController
from mbio.api.to_file.labelfree import *
from mainapp.libs.signature import check_sig
import unittest
import os


class RelateDiffAction(ProteinTranscriptLabelfreeController):
    def __init__(self):
        super(RelateDiffAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "gene_diff_id", 'gene_compare_group', 'protein_diff_id', 'protein_compare_group', 'species']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        if not data.gene_compare_group or not data.protein_compare_group:
            info = {'success': False, 'info': "差异分组为空"}
            return json.dumps(info)
        if not data.gene_diff_id or not data.protein_diff_id:
            info = {'success': False, 'info': "请选择要分析的差异表"}
            return json.dumps(info)
        relation_info = self.labelfree.get_relation_condition(data.task_id)
        if not relation_info:
            info = {'success': False, 'info': "关联信息为空，请确保改项目已经跟转录组项目做好关联"}
            return json.dumps(info)

        project_sn = relation_info["project_sn"]
        rna_type = relation_info['rna_type']
        rna_task_id = relation_info['rna_task_id']

        task_info = self.labelfree.get_task_info(data.task_id)
        if task_info.has_key('protein_fa'):
            seq = task_info['protein_fa']
        else:
            info = {'success': False, 'info': "找不到蛋白的序列文件"}
            return json.dumps(info)

        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            gene_diff_id=data.gene_diff_id,
            protein_diff_id=data.protein_diff_id,
            gene_compare_group=data.gene_compare_group,
            protein_compare_group=data.protein_compare_group,
            combine_score=data.combine_score,
            taxon=data.taxon,
            species=data.species
        )

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "relate_diff" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            rna_type=rna_type,
            name=name,
            gene_diff_id=data.gene_diff_id,
            protein_diff_id=data.protein_diff_id,
            gene_compare_group=data.gene_compare_group,
            protein_compare_group=data.protein_compare_group,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='the combine situation of the protein and gene',
            params=params,
            seq=seq,
            species=data.species,
            status="start",
            version='v2.2'
        )
        main_id = self.labelfree.insert_main_table('sg_p2g_diff', main_info)

        # prepare option for workflow
        options = dict(
            fc=','.join([data.protein_diff_id, data.gene_diff_id, rna_type]),
            diff_protein=data.protein_diff_id,
            diff_rna=rna_type,
            gene_diff_id=data.gene_diff_id,
            gene_compare_group=data.gene_compare_group,
            protein_compare_group=data.protein_compare_group,
            relation_tab=data.task_id,
            rna_type=rna_type,
            species=data.species,
            combine_score=data.combine_score,
            main_id=str(main_id),
            seq=seq,
            update_info=json.dumps({str(main_id): "sg_p2g_diff"})  # to update sg_status
        )

        # prepare to file
        to_files = ["protein_transcript_labelfree.export_diff_labelfree(diff_protein)",
                    "protein_transcript_labelfree.export_diff_rna(diff_rna)",
                    "protein_transcript_labelfree.export_relation_tab(relation_tab)",
                    "protein_transcript_labelfree.export_fc_info(fc)"
                    ]

                # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'protein_transcript_labelfree.report.relate_diff'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(RelateDiffAction, self).POST()
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
        cmd += "s/protein_transcript_labelfree/relate_diff "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_34739",
            task_type="2",
            submit_location="diffprogene",
            gene_compare_group="Con|Drug",
            protein_compare_group="S|F",
            gene_diff_id="5c218da2a4e1af619b5b6b40",
            protein_diff_id="5d22a83917b2bf1c689e18ee",
            species="9606",
            combine_score='300',
            taxon='all'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

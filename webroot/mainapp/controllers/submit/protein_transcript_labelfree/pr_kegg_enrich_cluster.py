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


class PrKeggEnrichClusterAction(ProteinTranscriptLabelfreeController):
    def __init__(self):
        super(PrKeggEnrichClusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "protein_kegg_id", 'rna_kegg_id', ]
        basic_args += ['geneset_id', 'proteinset_id', 'gct', 'padjust_way', 'cluster_col']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        relation_info = self.labelfree.get_relation_condition(data.task_id)
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
            gct=data.gct,
            cluster_col=data.cluster_col,
            padjust_way=data.padjust_way
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
        name = "relate_kegg_enrich_cluster" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            rna_type=rna_type,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='the kegg enrich cluster of proteinsets and genesets',
            params=params,
            status="start"
        )
        main_id = self.labelfree.insert_main_table('sg_pr_kegg_enrich_cluster', main_info)

        # prepare option for workflow
        options = dict(
            gene_list=rna_type + ';' + str(data.geneset_id),
            protein_list=str(data.proteinset_id),
            protein_kegg=str(data.protein_kegg_id),
            protein_info=str(data.protein_kegg_id),
            gene_kegg=rna_type + ',' + str(data.rna_kegg_id),
            gene_info=rna_type + ',' + str(data.rna_kegg_id),
            rna_type=rna_type,
            method=data.padjust_way,
            main_table_id=str(main_id),
            n_clusters=n_clusters,
            gcm=gcm,
            gct=data.gct,
            gcd=gcd,
            cluster_col=data.cluster_col,
            update_info=json.dumps({str(main_id): "sg_pr_kegg_enrich_cluster"})  # to update sg_status
        )


        # prepare to file
        to_files = ["protein_transcript_labelfree.export_proteinsets_list(protein_list)",
                    "protein_transcript_labelfree.export_genesets_list(gene_list)",
                    "protein_transcript_labelfree.export_kegg_table_rna(gene_kegg)",
                    "protein_transcript_labelfree.export_add_info_rna(gene_info)",
                    "protein_transcript_labelfree.export_kegg_table_protein(protein_kegg)",
                    "protein_transcript_labelfree.export_add_info_protein(protein_info)",
                    ]
        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'protein_transcript_labelfree.report.pr_kegg_enrich_cluster'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(PrKeggEnrichClusterAction, self).POST()
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
        cmd += "s/protein_transcript_labelfree/pr_kegg_enrich_cluster "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_34739",
            task_type="2",
            submit_location="prkeggenrichcluster",
            protein_kegg_id="5d22a83917b2bf1c689e1264",
            rna_kegg_id="5c218caaa4e1af619b4fd964",
            # proteinset_id="5c2103e0a4e1af61602f1f93",
            proteinset_id="5d22a83b17b2bf1c689e2825,5d22a83b17b2bf1c689e2827",
            # geneset_id="5c218da2a4e1af619b5b6b3e",
            geneset_id="5c218d9fa4e1af619b5b6b3c,5c218da2a4e1af619b5b6b3e",
            padjust_way='BH',
            n_clusters='10',
            gct="hierarchy",
            gcm="average",
            gcd="euclidean",
            cluster_col='Ratio_in_study'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
# __author__ = "fengyitong 2018-12-11"

from biocluster.workflow import Workflow
import re
from collections import defaultdict
import unittest
from biocluster.wpm.client import *
import datetime
import copy
import os
import pandas as pd
from collections import OrderedDict
from mbio.packages.labelfree.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json
import glob


class RelasetClusterWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RelasetClusterWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='protein_matrix', type='string'),
            dict(name='rna_matrix', type='string'),
            dict(name='relaset_list', type='string'),
            dict(name='rna_group_dict', type='string'),
            dict(name='protein_group_dict', type='string'),
            dict(name='rna_group_id', type='string'),
            dict(name='protein_group_id', type='string'),
            dict(name='diff_rna', type='string'),
            dict(name="n_clusters", type='int'),
            dict(name="use_group", type="string"),
            dict(name="group", type="string"),
            dict(name="gct", type="string"),
            dict(name="gcm", type="string"),
            dict(name="gcd", type="string"),
            dict(name="rna_type", type='string', default="ref_rna_v2"),
            dict(name="main_id", type='string'),
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/7_relate/02_relate_exp/02_relate_cluster')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(RelasetClusterWorkflow, self).send_log(data)

    def merge_matrix(self):
        def get_rel_col(inlist, rel_list, n):
            rel_copy = copy.copy(rel_list)
            outlist = list()
            for i in inlist:
                for rel in rel_copy:
                    if i == rel.split('|')[n]:
                        outlist.append(rel)
                        rel_copy.remove(rel)
                        break
            return outlist
        #由于图片展示一边全是蓝一边全是红，所以加入了线下的归一化的方法
        def normal(exp_file,out):
            with open(exp_file) as er, open(out,'w') as ow:
                ow.write(er.readline())
                for line in er:
                    line = line.strip().split('\t')
                    if line:
                        acc = line[0]
                        line = [float(x) for x in line[1:]]
                        norm = [str(x / (sum(line) + 0.01) * len(line) + 0.1) for x in line]
                        ow.write(acc + '\t' + '\t'.join(norm) + '\n')
            return out
        with open(self.option('relaset_list'),'r') as rel_r:
            rel=rel_r.readlines()
            proteins = [x.strip().split('|')[0] for x in rel]
            transcripts = [x.strip().split('|')[1] for x in rel]
            rels = [x.strip() for x in rel]
        self.merge_out = os.path.join(self.output_dir,'merge_matrix')

        p_normaled = normal(self.option('protein_matrix'), os.path.join(self.work_dir, 'protein_normaled.txt'))
        r_normaled = normal(self.option('rna_matrix'), os.path.join(self.work_dir, 'rna_normaled.txt'))
        # protein_matrix=pd.read_table(self.option('protein_matrix') ,index_col=0, dtype={0:str})
        protein_matrix = pd.read_table(p_normaled, index_col=0, dtype={0: str})
        protein_matrix = protein_matrix.loc[proteins, :]
        # rna_matrix = pd.read_table(self.option('rna_matrix'), index_col=0)
        rna_matrix = pd.read_table(r_normaled, index_col=0)
        rna_matrix = rna_matrix.loc[list(set(transcripts)), :]
        rename_protein = OrderedDict()
        rename_rna = OrderedDict()
        for col in protein_matrix.columns:
            rename_protein[col] = col + '_protein'
        protein_matrix.rename(columns=rename_protein, inplace=True)
        for col in rna_matrix.columns:
            rename_rna[col] = col + '_rna'
        rna_matrix.rename(columns=rename_rna, inplace=True)
        merge_matrix = pd.DataFrame(columns=protein_matrix.columns.tolist() + rna_matrix.columns.tolist())
        for rel in rels:
            rel_dict = protein_matrix.loc[rel.split('|')[0],].to_dict()
            rel_dict.update(rna_matrix.loc[rel.split('|')[1],].to_dict())
            rel_dict['rela_id'] = rel
            merge_matrix = merge_matrix.append(rel_dict, ignore_index=True)
        merge_matrix = merge_matrix.set_index('rela_id')


        # p_normaled = normal(self.option('protein_matrix'), os.path.join(self.work_dir,'protein_normaled.txt'))
        # r_normaled = normal(self.option('rna_matrix'), os.path.join(self.work_dir,'rna_normaled.txt'))
        # # protein_matrix=pd.read_table(self.option('protein_matrix') ,index_col=0, dtype={0:str})
        # protein_matrix=pd.read_table(p_normaled, index_col=0, dtype={0:str})
        # protein_matrix=protein_matrix.loc[proteins,:]
        # # rna_matrix = pd.read_table(self.option('rna_matrix'), index_col=0)
        # rna_matrix = pd.read_table(r_normaled, index_col=0)
        # rna_matrix = rna_matrix.loc[transcripts, :]
        # protein_matrix['rela_id'] = get_rel_col(protein_matrix.index.tolist(), rels, 0)
        # rna_matrix['rela_id'] = get_rel_col(rna_matrix.index.tolist(), protein_matrix['rela_id'].tolist(),1)
        # protein_matrix = protein_matrix.set_index('rela_id')
        # rna_matrix = rna_matrix.set_index('rela_id')
        # rename_protein = OrderedDict()
        # rename_rna=OrderedDict()
        # for col in protein_matrix.columns:
        #     rename_protein[col] = col + '_protein'
        # protein_matrix.rename(columns=rename_protein, inplace=True)
        # for col in rna_matrix.columns:
        #     rename_rna[col] = col + '_rna'
        # rna_matrix.rename(columns=rename_rna, inplace=True)
        # merge_matrix = protein_matrix.join(rna_matrix, how='inner')
        merge_matrix.to_csv(self.merge_out, sep='\t', header=True, index=True)

    def run_tool(self):
        self.tool = self.add_tool("labelfree.exp_cluster")
        self.tool.on("end", self.set_db)
        options = dict(
            exp=self.merge_out,
            group=self.option('group'),
            n_clusters=int(self.option('n_clusters')),
            sct='no',
            gct=self.option('gct'),
            gcm=self.option('gcm'),
            gcd=self.option('gcd'),
            use_group=self.option('use_group'),
        )
        self.tool.set_options(options)
        self.tool.run()

    def run(self):
        self.merge_matrix()
        self.run_tool()
        super(RelasetClusterWorkflow, self).run()

    def set_db(self):
        cluster = self.api.api("protein_transcript_labelfree.relaset")
        # add result info
        cluster.add_relaset_cluster(self.tool.output_dir, main_id=self.option('main_id'))
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        pt_expression_matrix = os.path.join(self.tool.output_dir, "expression_matrix.xls")
        pt_seq_tree = os.path.join(self.tool.output_dir, "seq.cluster_tree.txt")
        if os.path.exists(pt_expression_matrix) and os.path.exists(pt_seq_tree):
            chart.chart_pt_proteinsetcluster(pt_expression_matrix, pt_seq_tree)
            chart.to_pdf()
            # move pdf to result dir
            pdf_file = os.path.join(self.work_dir, "pt_proteinsetcluster.heatmap_new.pdf")
            if os.path.exists(pdf_file):
                os.link(pdf_file, os.path.join(self.tool.output_dir, "heat.pdf"))
            pdf_files = glob.glob(self.work_dir + "/pt_proteinsetcluster__*.pdf")
            for pdf_file in pdf_files:
                os.link(pdf_file, os.path.join(self.tool.output_dir, "sub"+os.path.basename(pdf_file)[22:].split('.')[0].split('_')[1]+'.pdf'))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["7_relate", "", "蛋白组与转录组关联分析",0],
            ["7_relate/02_relate_exp", "", "表达量信息", 0],
            ["7_relate/02_relate_exp/02_relate_cluster", "", "关联数据表达量聚类分析", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "关联集聚类分析结果目录"],
            ["./heat.pdf", "", "热图"],
            ["./*sub*.pdf", "", "子聚类图"],
        ])
        super(RelasetClusterWorkflow, self).end()

class TestFunction(unittest.TestCase):
    def test(self):
        worker = worker_client()
        id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        data = {
            "id": "protein_transcript_labelfree_" + id,
            "type": "workflow",
            "name": "protein_transcript_labelfree.protein_transcript_labelfree",
            "options": dict(
                rna_type="ref_rna_v1",
                gff="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/biomart/Homo_sapiens.GRCh37.biomart",
                pep="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/cds/Homo_sapiens.GRCh37.pep.fa",
                protein_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/protein.list",
                transcript_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/gene.list",
                protein_faa="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/DB.fasta",
                )
        }

        info = worker.add_task(data)
        print(info)


    # def test(self):
    #     import random
    #     from mbio.workflows.single import SingleWorkflow
    #     from biocluster.wsheet import Sheet
    #     data = {
    #         "id": "protein_transcript_labelfree",
    #         "type": "workflow",
    #         "name": "protein_transcript_labelfree.protein_transcript_labelfree",
    #         "options": dict(
    #             rna_type="ref_rna_v1",
    #             gff="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/biomart/Homo_sapiens.GRCh37.biomart",
    #             pep="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/cds/Homo_sapiens.GRCh37.pep.fa",
    #             protein_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/protein.list",
    #             transcript_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/gene.list",
    #             protein_faa="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/DB.fasta",
    #             )
    #     }
    #     # data['options']['method'] = 'rsem'
    #     # wsheet = Sheet(data=data)
    #     # wf = SingleWorkflow(wsheet)
    #     # wf.run()
    #     # #
    #     # data['id'] += '1'
    #     # data['options']['method'] = 'salmon'
    #     # wsheet = Sheet(data=data)
    #     # wf = SingleWorkflow(wsheet)
    #     # wf.run()
    #     #
    #     data['id'] += '_fyt'
    #     wsheet = Sheet(data=data)
    #     wf = SingleWorkflow(wsheet)
    #     wf.run()


if __name__ == '__main__':
    unittest.main()

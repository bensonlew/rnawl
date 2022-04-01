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
import json
from mbio.packages.dia_v3.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder



class SampleProteinGeneWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SampleProteinGeneWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='protein_matrix', type='string'),
            dict(name='rna_matrix', type='string'),
            dict(name='relation_tab', type='string'),
            dict(name="rna_group_dict", type='string'),
            dict(name="protein_group_dict", type='string'),
            dict(name="value_type", type='string', default="fpkm"),
            dict(name="rna_cutoff", type='string', default="1"),
            dict(name="protein_cutoff", type='string', default="0"),
            dict(name="rna_type", type='string', default="ref_rna_v2"),
            dict(name="main_id", type='string'),
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.filter_exp_matrix(self.option('rna_matrix'), cutoff = self.option('rna_cutoff'))
        self.filter_exp_matrix(self.option('protein_matrix'), cutoff = self.option('protein_cutoff'))
        with open(self.option('protein_matrix'), 'r') as pro, \
                open(self.option('rna_matrix'), 'r') as trans, \
                open(self.option('relation_tab'), 'r') as rel_t:
            _ = pro.readline()
            _ = trans.readline()
            self.protein_list = [x.split('\t')[0].strip() for x in pro.readlines()]
            self.transcript_list = [x.split('\t')[0].strip() for x in trans.readlines()]
            for line in rel_t:
                line = line.strip().split('\t')
                if line[0] == 'related' and len(line) >1:
                    self.related_list = line[1].split(';')
        self.copy_protein = copy.copy(self.protein_list)
        print(len(self.protein_list))
        print(len(self.transcript_list))
        self.remove_trans = list()
        self.relation_list = list()
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/7_relate/01_relate_sample')
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
        super(SampleProteinGeneWorkflow, self).send_log(data)

    def run(self):
        self.start_listener()
        self.fire("start")
        self.protein_gene()
        # super(SampleProteinGeneWorkflow, self).run()

    def set_db(self):
        with open(self.output_dir + '/relationship_list', 'w') as rel_w:
            rel_w.write('related' + '\t' + ';'.join(self.relation_list) + '\n')
            rel_w.write('failed_proteins' + '\t' + ';'.join(self.protein_list) + '\n')
            rel_w.write('failed_transcripts' + '\t' + ';'.join(self.transcript_list) + '\n')
        relation = self.api.api("protein_transcript_dia_v3.protein_transcript")
        params = dict(
            protein_group_dict=self.option('protein_group_dict'),
            rna_group_dict=self.option('rna_group_dict'),
            value_type=self.option('value_type'),
            rna_cutoff=self.option('rna_cutoff'),
            protein_cutoff=self.option('protein_cutoff'),
        )
        print(len(self.protein_list))
        print(len(self.transcript_list))
        relation.add_relation(main_id = self.option('main_id'), params = params, relation_file= self.output_dir + '/relationship_list', table = 'sg_sample_protein_gene', detail_id= 'p2g_sample_id')
        """
        保存结果表到mongo数据库中
        """
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        chart.output_dir = self.work_dir + '/'

        protein_transcript_baseinfo = os.path.join(self.output_dir, "relationship_list")
        if os.path.exists(protein_transcript_baseinfo):
            chart.chart_protein_transcript_baseinfo(protein_transcript_baseinfo)
            chart.to_pdf()

            # move pdf to result dir
            pdf_file = os.path.join(self.work_dir, "protein_transcript_baseinfo.venn2.pdf")
            if os.path.exists(pdf_file):
                os.link(pdf_file, os.path.join(self.output_dir, "venn.pdf"))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["7_relate", "", "蛋白组与转录组关联分析", 0],
            ["7_relate/01_relate_sample", "", "关联基本信息", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "联合分析样本数据关联结果"],
            ["relationship_list", "txt", "联合分析样本数据关联结果表"],
            ["venn.pdf", "", "关联基本信息venn图"],
        ])
        super(SampleProteinGeneWorkflow, self).end()

    def protein_gene(self):
        for protein in self.copy_protein:
            for rela in self.related_list:
                if rela.split('|')[0] == protein:
                    gene = rela.split('|')[1]
                    if gene in self.transcript_list or gene in self.remove_trans:
                        self.relation_list.append(rela)
                        self.protein_list.remove(protein)
                        try:
                            self.transcript_list.remove(gene)
                            self.remove_trans.append(gene)
                        except:
                            self.remove_trans.append(self.transcript_list.pop())
        self.set_db()

    def filter_exp_matrix(self, matrix, cutoff):
        filter_list = list()
        with open(matrix, 'r') as mat_r:
            header = mat_r.readline()
            for line in mat_r:
                tmp = line.strip().split('\t')
                cut = 0
                for i in tmp[1:]:
                    if float(i) < float(cutoff):
                        cut += 1
                if not cut > float(len(tmp))/2:
                    filter_list.append(line)
        with open(matrix, 'w') as mat_w:
            mat_w.write(header)
            for line in filter_list:
                mat_w.write(line)

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


if __name__ == '__main__':
    unittest.main()

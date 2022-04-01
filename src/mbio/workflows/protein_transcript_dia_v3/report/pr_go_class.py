# -*- coding: utf-8 -*-
# __author__ = "fengyitong 2019-01-14"

from biocluster.workflow import Workflow
import unittest
from biocluster.wpm.client import *
import datetime
import re
import os
from mbio.packages.dia_v3.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import glob
import json

# from collections import OrderedDict


class PrGoClassWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PrGoClassWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='gene_go', type='string'),
            dict(name='protein_go', type='string'),
            dict(name='protein_go_id', type='string'),
            dict(name='gene_go_id', type='string'),
            dict(name="rna_type", type='string', default="ref_rna_v2"),
            dict(name="main_id", type='string'),
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/7_relate/03_relate_anno/01_relate_go')
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
        super(PrGoClassWorkflow, self).send_log(data)

    def merge_pro_gene_go(self):
        merge_dict = dict()
        term_list = ["biological_process", "cellular_component", "molecular_function"]
        for term in term_list:
            merge_dict[term] = dict()
            merge_dict[term]['2'] = dict()
            merge_dict[term]['3'] = dict()
            merge_dict[term]['4'] = dict()
        with open(self.option('protein_go'), 'r') as pro_go:
            _ = pro_go.readline()
            for line in pro_go:
                line = line.strip().split('\t')
                term_type, term, go_id, level, protein_num, protein_percent, protein_list = line
                if go_id not in  merge_dict[term_type][level]:
                    merge_dict[term_type][level][go_id] = dict(
                        term_type = term_type,
                        term = term,
                        go_id = go_id,
                        level = level,
                        protein_num = protein_num,
                        protein_percent = protein_percent,
                        protein_list = protein_list,
                        gene_num = '0',
                        gene_percent = '0',
                        gene_list = ''
                    )
                info_dict = merge_dict[term_type][level][go_id]
                info_dict['protein_num'] = protein_num
                info_dict['protein_percent'] = protein_percent
                info_dict['protein_list'] = protein_list

        with open(self.option('gene_go'), 'r') as gene_go:
            _ = gene_go.readline()
            for line in gene_go:
                line = line.strip().split('\t')
                term_type, term, go_id, level, gene_num, gene_percent, gene_list = line
                if go_id not in  merge_dict[term_type][level]:
                    merge_dict[term_type][level][go_id] = dict(
                        term_type = term_type,
                        term = term,
                        go_id = go_id,
                        level = level,
                        protein_num = '0',
                        protein_percent = '0',
                        protein_list = '',
                        gene_num = gene_num,
                        gene_percent = gene_percent,
                        gene_list = gene_list
                    )
                info_dict = merge_dict[term_type][level][go_id]
                info_dict['gene_num'] = gene_num
                info_dict['gene_percent'] = gene_percent
                info_dict['gene_list'] = gene_list
        self.merge_go_file = os.path.join(self.output_dir, 'pro_gene_go_class.xls')
        with open(self.merge_go_file, 'w') as go_w:
            go_w.write("Term type\tTerm\tGO\tLevel\t" + "protein_num\tprotein_percent\tprotein_list\t" + "gene_num\tgene_percent\tgene_list" + "\n")
            for term in term_list:
                for level in merge_dict[term]:
                    for go in merge_dict[term][level]:
                        tmp_dict = merge_dict[term][level][go]
                        tmp_str = '\t'.join(['{' + key + '}' for key in ['term_type', 'term', 'go_id',
                                                                         'level', 'protein_num', 'protein_percent',
                                                                         'protein_list', 'gene_num', 'gene_percent',
                                                                         'gene_list']])
                        go_w.write(tmp_str.format(**tmp_dict) + '\n')
        self.set_db()

    def run(self):
        self.start_listener()
        self.fire("start")
        self.merge_pro_gene_go()
        # super(PrGoClassWorkflow, self).run()

    def set_db(self):
        go = self.api.api("protein_transcript_dia_v3.relaset")
        # add result info
        go.add_pr_go_class(self.merge_go_file, main_id=self.option('main_id'))
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        pt_goclass = os.path.join(self.output_dir, "pro_gene_go_class.xls")
        if os.path.exists(pt_goclass):
            chart.chart_pt_goclass(pt_goclass)
            chart.to_pdf()
            # move pdf to result dir
            pdf_file = os.path.join(self.work_dir, "pt_goclass.go_bar.pdf")
            if os.path.exists(pdf_file):
                os.link(pdf_file, os.path.join(self.output_dir, "bar.pdf"))
            pdf_file = os.path.join(self.work_dir, "pt_goclass.double_nest.pdf")
            if os.path.exists(pdf_file):
                os.link(pdf_file, os.path.join(self.output_dir, "pie.pdf"))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["7_relate", "", "蛋白组与转录组关联分析",0],
            ["7_relate/03_relate_anno", "", "功能注释信息", 0],
            ["7_relate/03_relate_anno/01_relate_go", "", "GO注释", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "转录蛋白结合展示GO分类目录"],
            ["bar.pdf", "", "直方图"],
            ["pie.pdf", "", "饼图"],
        ])
        super(PrGoClassWorkflow, self).end()

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

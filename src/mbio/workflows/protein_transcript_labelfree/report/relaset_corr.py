# -*- coding: utf-8 -*-
# __author__ = "fengyitong 2018-12-11"

from biocluster.workflow import Workflow
import unittest
from biocluster.wpm.client import *
import datetime
import re
from mbio.packages.labelfree.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import os
import json
import glob


class RelasetCorrWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RelasetCorrWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='protein_matrix', type='string'),
            dict(name='rna_matrix', type='string'),
            dict(name='relaset_list', type='string'),
            dict(name='rna_group_dict', type='string'),
            dict(name='protein_group_dict', type='string'),
            dict(name="group", type="string"),
            dict(name="rna_type", type='string', default="ref_rna_v2"),
            dict(name="corr_method", type='string', default="spearman"),
            dict(name="use_group", type='string', default="no"),
            dict(name="exp_deal", type='string', default="no"),
            dict(name="main_id", type='string'),
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/7_relate/02_relate_exp/03_relate_corr')
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
        super(RelasetCorrWorkflow, self).send_log(data)

    def run_tool(self):
        self.tool = self.add_tool("protein_transcript_labelfree.pro_gene_corr")
        self.tool.on("end", self.set_db)
        options = dict(
            group=self.option('group'),
            protein_matrix=self.option('protein_matrix'),
            rna_matrix=self.option('rna_matrix'),
            relaset_list=self.option('relaset_list'),
            corr_method=self.option('corr_method'),
            use_group=self.option('use_group'),
            exp_deal=self.option('exp_deal'),
        )
        self.tool.set_options(options)
        self.tool.run()

    def run(self):
        self.run_tool()
        super(RelasetCorrWorkflow, self).run()

    def set_db(self):
        corr = self.api.api("protein_transcript_labelfree.relaset")
        # add result info
        corr.add_relaset_corr(self.tool.output_dir, main_id=self.option('main_id'),workflow_out = self.get_workflow_output_dir())
        self.end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if re.match(r'tsanger:', workflow_output):
            workflow_output = workflow_output.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        elif re.match(r'^\w+://\S+/.+$', workflow_output):
            workflow_output = workflow_output
        else:
            workflow_output = workflow_output.replace('sanger:', '/mnt/ilustre/data/')
        return workflow_output

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        pt_corrScatter = os.path.join(self.tool.output_dir, "protein_rna_density_matrix.xls")
        pt_corrScatter_pvalue = os.path.join(self.tool.output_dir, "corr_info")
        if os.path.exists(pt_corrScatter) and os.path.exists(pt_corrScatter_pvalue):
            chart.chart_pt_proteinsetcluster_corr(pt_corrScatter, pt_corrScatter_pvalue)
            chart.to_pdf()
            # move pdf to result dir
            pdf_file = os.path.join(self.work_dir, "pt_proteinsetcluster_corr.scatter_new.pdf")
            if os.path.exists(pdf_file):
                os.link(pdf_file, os.path.join(self.tool.output_dir, "single.pdf"))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["7_relate", "", "蛋白组与转录组关联分析",0],
            ["7_relate/02_relate_exp", "", "表达量信息", 0],
            ["7_relate/02_relate_exp/03_relate_corr", "", "关联数据表达相关性分析", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "关联集相关性分析结果目录"],
            ["./single.pdf", "", "单组相关性图"],
            ["./multily_corr.pdf", "", "多组相关性图"],
        ])
        super(RelasetCorrWorkflow, self).end()

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

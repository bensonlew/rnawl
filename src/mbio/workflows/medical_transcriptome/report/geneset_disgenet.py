# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.workflow import Workflow
import os
import pandas as pd
import unittest
import re
import json
import glob
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import types
from bson.objectid import ObjectId
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.chart.chart_geneset import ChartGeneset


class GenesetDisgenetWorkflow(Workflow):
    """
    基因集DisGeNET富集分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetDisgenetWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gene_list", "type": "infile", "format": "medical_transcriptome.gene_list"},
            {"name": "entrez_list", "type": "infile", "format": "medical_transcriptome.common"},
            {"name": "padjust_method", "type": "string", "default": "BH"},
            {"name": "dsi", "type": "string", "default": None},
            {"name": "dpi", "type": "string", "default": None},
            {"name": "score", "type": "string", "default": "0.5"},
            {"name": "el", "type": "string", "default": "Definitive"},   # comma between terms
            {"name": "ei", "type": "string", "default": None},
            {"name": "main_table_data", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "main_id", "type": "string", "default": None},
            {"name": "task_type", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "geneset_name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.entrez_covert = self.add_tool("medical_transcriptome.geneset.gene2entrez")
        self.disgenet_enrich = self.add_tool("medical_transcriptome.geneset.disgenet_enrich")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/03 Enrich/05 DisGeNET')
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
        super(GenesetDisgenetWorkflow, self).send_log(data)

    def run(self):
        self.entrez_covert.on("end", self.run_disgenet_enrich)
        self.get_run_log()
        self.run_entrez_convert()
        super(GenesetDisgenetWorkflow, self).run()

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        disgenet_enrich_table = self.disgenet_enrich.output_dir + '/DisGeNET_enrich_stat.xls'
        with open(disgenet_enrich_table, "r") as f:
            lines = f.readlines()
            if len(lines) >= 2:
                chart.chart_geneset_enrich_disgenet(disgenet_enrich_table, geneset_name=self.option('geneset_name'))
                chart.to_pdf()
                # move pdf to result dir
                pdf_file = glob.glob(self.work_dir + "/*.pdf")
                for p in pdf_file:
                    os.link(p, self.disgenet_enrich.output_dir + "/" + os.path.basename(p))

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_geneset_disgenet_enrich", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_disgenet = self.api.api('medical_transcriptome.geneset_disgenet')
        self.logger.info("开始进行disgenet_enrich的导表")
        output_file = os.path.join(self.disgenet_enrich.work_dir, "enrichment_result.txt")
        g2e_path = os.path.join(self.entrez_covert.work_dir, "gene2entrez.list")
        if not os.path.exists(output_file):
            api_disgenet.update_main_table_status(disgenet_id=self.option("main_id"))
        if os.path.exists(output_file):
            api_disgenet.add_annotation_detail(enrich_path=output_file,
                                               g2e_path=g2e_path,
                                               disgenet_id=self.option("main_id"))
        self.end()

    def run_entrez_convert(self):
        opts = {
            "gene_list": self.option("gene_list"),
            "entrez_list": self.option("entrez_list"),
            "geneset_id": self.option("geneset_id"),
        }
        self.entrez_covert.set_options(opts)
        self.entrez_covert.run()

    def run_disgenet_enrich(self):
        gene_list = os.path.join(self.entrez_covert.work_dir, "{}_entrez.list".format(self.option('geneset_id')))
        g2e_path = os.path.join(self.entrez_covert.work_dir, "gene2entrez.list")
        opts = {
            "gene_list": gene_list,
            "dsi": self.option("dsi"),
            "dpi": self.option("dpi"),
            "score": self.option("score"),
            "el": self.option("el"),
            "ei": self.option("ei"),
            "padjust_method": self.option("padjust_method"),
            'g2e_path': g2e_path,
        }
        self.disgenet_enrich.set_options(opts)
        self.disgenet_enrich.on("end", self.set_db)
        self.disgenet_enrich.run()

    def end(self):
        enrich_path = self.disgenet_enrich.output_dir + "/enrichment_result.txt"
        if os.path.exists(enrich_path):
            enrichment = pd.read_table(enrich_path, sep="\t", header=0)
            enrichment = enrichment.drop(["qvalue"], axis=1)
            enrichment = enrichment.rename(columns={"GeneRatio": "ratio_in_study", "BgRatio": "ratio_in_pop"})
            enrichment.to_csv(self.disgenet_enrich.output_dir + "/DisGeNET_enrich_stat.xls",
                              sep="\t",  header=True, index=False,)
            os.remove(enrich_path)
        self.chart()
        # set_output
        if os.path.exists(os.path.join(self.disgenet_enrich.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.disgenet_enrich.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.disgenet_enrich.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.disgenet_enrich.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录", 0],
            ["04 GeneSet/03 Enrich", "", "基因集功能富集", 0],
            ["04 GeneSet/03 Enrich/05 DisGeNET", "", "基因集DisGeNET功能富集", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "基因集DisGeNET富集分析文件", 0, ],
            ["./DisGeNET_enrich_stat.xls", "xls", "DisGeNET富集分析统计表", 0, ],
            ["./DisGeNET_enrich_detail.xls", 'xls', '基因/转录本对应DisGeNET富集详情表', 0, ],
            ["*bar.pdf", "pdf", "DisGeNET富集分析柱形图", 0],
            ["*bar_line.pdf", "pdf", "DisGeNET富集分析柱形图(带折线)", 0],
            ["*buble.pdf", "pdf", "DisGeNET富集分析气泡图", 0],
            ["*buble2.pdf", "pdf", "DisGeNET富集分析气泡图(分散型)", 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        super(GenesetDisgenetWorkflow, self).end()




class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.medical_transcriptome.report.geneset_disgenet import GenesetDisgenetWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'geneset_disgenet_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'medical_transcriptome.report.geneset_disgenet',
            'options': {
                "gene_list": "/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/gene_test.list",
                "padjust_method": "BH",
                "dsi": "0.1",
                "dpi": "0.1",
                "score": "0.1",
                "el": "definitive,strong,limited,moderate",
                "geneset_id": "5f45c70417b2bf78d9c9c174",
                "task_id": "tsg_medical_transcriptome",
                "project_sn": "188_5d01dede4f911",
            }
        }
        wsheet_object = Sheet(data=data)
        wf = GenesetDisgenetWorkflow(wsheet_object)
        wf.sheet.id = 'tsg_medical_transcriptome'
        wf.sheet.project_sn = '188_5d01dede4f911'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()



if __name__ == '__main__':
    unittest.main()

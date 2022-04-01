# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# last update liubinxu
from biocluster.workflow import Workflow
from collections import defaultdict
import os
import time
import re
from itertools import chain
from mbio.packages.prok_rna.chart import Chart
import glob
import json

class GenesetClassWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetClassWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "geneset_go", "type": "string"},
            {"name": "gene_go", "type": "string", "default":None},
            {"name": "geneset_cog", "type": "string"},
            # {"name": "geneset_go2", "type": "string"},
            # {"name": "geneset_kegg", "type": "string"},
            # {"name": "kegg_table", "type": "infile", "format": "denovo_rna_v2.kegg_table"},
            {"name": "anno_type", "type": "string"},
            {"name": "geneset_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "version", "type": "string"},
            {"name": "type", "type": "string"},
            {"name": "geneset_name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.start_listener()
        self.fire("start")
        self.set_db()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('prok_rna.geneset')
        self.logger.info("注释分类分析结果入库")
        if self.option("anno_type") == "cog":
            time.sleep(30)
            output_file = self.option("geneset_cog")
            os.link(output_file, self.output_dir + "/cog_class_table.xls")
            api_geneset.add_geneset_cog_detail(output_file, self.option("main_table_id"))
        elif self.option("anno_type") == "go":
            time.sleep(30)
            output_file = self.option("geneset_go")
            os.link(output_file, self.output_dir + "/go_class_table.xls")
            api_geneset.add_go_regulate_detail(output_file, self.option("main_table_id"), gene_go=self.option("gene_go"))
        else:
            output_file = self.output_dir + '/kegg_stat.xls'
            pathway_file = self.output_dir + '/pathways'
            api_geneset.add_kegg_regulate_detail(self.option("main_table_id"), output_file)
            api_geneset.add_kegg_regulate_pathway(pathway_file, self.option("main_table_id"))
        # os.link(output_file, self.output_dir + "/" + os.path.basename(output_file))
        print(output_file)
        self.end()

    def move_pdf(self, origin, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(origin):
            os.link(origin, new)

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        if self.option("anno_type") == "go":
            geneset_go_class = os.path.join(self.output_dir, 'go_class_table.xls')
            chart.prok_geneset_goclass(geneset_go_class)
        elif self.option("anno_type") == "cog":
            geneset_list = self.option('geneset_name').split(',')
            geneset_cog = os.path.join(self.output_dir, 'cog_class_table.xls')
            chart.chart_geneset_class_cog(geneset_cog, geneset_list)
        chart.to_pdf()

        # move pdf
        target_dir = self.output_dir
        if self.option("anno_type") == "go":
            pdf = glob.glob(os.path.join(self.work_dir, '*go_bar_geneset.go_bar.pdf'))
            for i in pdf:
                if os.path.basename(i).split('go_bar_geneset.go_bar.pdf')[0] == "":
                    prefix = 'all'
                else:
                    prefix = os.path.basename(i).split('go_bar_geneset.go_bar.pdf')[0]
                self.move_pdf(i, os.path.join(target_dir, prefix + '_go_class.pdf'))
        elif self.option("anno_type") == "cog":
            pdf = os.path.join(self.work_dir, 'cog_annot.gene_set.column.pdf')
            self.move_pdf(pdf, os.path.join(target_dir, 'cog_class.pdf'))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("anno_type") == "go":
            result_dir.add_relpath_rules([
                [".", "", "基因集GO分类结果目录"],
                ["go_class_table.xls", "", "GO分类统计表"],
                # ["./estimators.xls", "xls", "alpha多样性指数表"]
            ])
            result_dir.add_regexp_rules([
                [r".*go_class\.pdf", "", "GO分类统计柱形图"],
            ])
        elif self.option("anno_type") == "cog":
            result_dir.add_relpath_rules([
                [".", "", "基因集COG分类结果目录"],
                ["cog_class_table.xls", "", "COG分类统计表"],
                ["cog_class.pdf", "", "COG分类统计柱状图"],
            ])
        # print self.get_upload_files()
        # self.set_end()
        # self.fire('end')
        # self.end_unfinish_job()
        # self._upload_result()
        # self._import_report_data()
        # self._update("set_end")
        # self.step.finish()
        # self.step.update()
        # self.logger.info("运行结束!")
        # self._save_report_data()
        super(GenesetClassWorkflow, self).end()

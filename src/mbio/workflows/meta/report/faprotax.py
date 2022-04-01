# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
import glob
import os
from bson import ObjectId
import re
import json
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class FaprotaxWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FaprotaxWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile","format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "group_table", "type": "infile","format": "meta.otu.group_table"},
            {"name": "otu_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "otu_id", "type": 'string'},
            {"name": "group_detail", "type": "string"},
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.faprotax = self.add_tool('meta.faprotax_predict')
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False


    def run_faprotax_predict(self):
        with open(self.option("group_table").prop['path'],"r") as f,open(self.work_dir +"/group.txt","w") as t:
            data = f.readlines()
            t.write("#sample\t#group\n")
            group_detail = json.loads(self.option("group_detail"))
            if len(group_detail) == 1 and "All" in group_detail:
                for i in data[1:]:
                    t.write(i.split("\t")[0] + "\t" + "All" + "\n")
            else:
                for i in data[1:]:
                    t.write(i)
        mapping_file = self.work_dir + "/mapping_file"
        self.faprotax.set_options({
            "otu_table": self.option("otu_table")
            })
        self.faprotax.on("end", self.set_db)
        self.faprotax.run()

    def run(self):
        self.run_faprotax_predict()
        super(FaprotaxWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        link_dir(self.faprotax.output_dir,self.output_dir)
        self.logger.info("正在写入mongo数据库")
        self.api_faprotax = self.api.api("faprotax")
        self.api_faprotax.add_Faprotax(prediction_file=self.output_dir + "/function_prediction.txt",group_file=self.work_dir +"/group.txt",main_id=self.option("main_id"))
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "sg_faprotax")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "faprotax",
                "interaction": 1,
                "main_table": "sg_faprotax",
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Faprotax功能预测结果目录",0],
            ["./function_prediction.txt", "txt", "Faprotax功能预测结果表",0],
            ["./report.txt", "txt", "Faprotax功能预测报告文件",0],
            ["./FAPROTAX功能Heatmap图.pdf", "pdf", "FAPROTAX功能Heatmap图", 0],
        ])
        super(FaprotaxWorkflow, self).end()
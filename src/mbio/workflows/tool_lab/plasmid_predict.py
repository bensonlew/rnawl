# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.tool_lab.common_function import fasta_dir_sort
import os
import subprocess
import gevent

class PlasmidPredictWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PlasmidPredictWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "file_style", "type": "string"},
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.tools = []

    def check_options(self):
        if (not self.option("fasta").is_set) and (not self.option("fasta_dir").is_set):
            raise OptionError("请传入序列文件或序列文件夹！", code="")

    def run(self):
        if self.option("fasta_dir").is_set:
            fasta_dir_sort(self.option("fasta_dir").prop['path'],self.work_dir + "/fasta_dir_sort")
        self.run_anno()
        super(PlasmidPredictWorkflow, self).run()

    def run_anno(self):
        """
        质粒预测
        :return:
        """
        if self.option("fasta").is_set:
            self.plasmid_predict = self.add_tool("tool_lab.plasmid_predict")
            opts = {
                "fasta": self.option("fasta"),
                "sample_name": self.option("fasta").prop['path'].split("/")[-1].split(".")[0]
            }
            self.plasmid_predict.set_options(opts)
            self.plasmid_predict.on('end', self.set_output)
            self.plasmid_predict.run()
        else:
            files = os.listdir(os.path.join(self.work_dir, "fasta_dir_sort"))
            for file in files:
                self.plasmid_predict = self.add_tool("tool_lab.plasmid_predict")
                opts= {
                    "fasta": os.path.join(os.path.join(self.work_dir,"fasta_dir_sort"),file),
                    "sample_name": file.rstrip(".fasta")
                }
                self.plasmid_predict.set_options(opts)
                self.tools.append(self.plasmid_predict)
            self.logger.info(self.tools)
            if self.tools:
                if len(self.tools) > 1:
                    self.on_rely(self.tools, self.set_output)
                elif len(self.tools) == 1:
                    self.tools[0].on('end', self.set_output)
                for tool in self.tools:
                    self.logger.info(tool)
                    gevent.sleep(1)
                    tool.run()

    def set_output(self):
        if self.option("fasta").is_set:
            link_dir(self.plasmid_predict.output_dir, self.output_dir)
        else:
            for tool in self.tools:
                self.logger.info(tool)
                for file in os.listdir(tool.output_dir):
                    if os.path.exists(self.output_dir + "/" + file):
                        os.remove(self.output_dir + "/" + file)
                    os.link(tool.output_dir + "/" + file,self.output_dir + "/" + file)
        self.set_db()

    def set_db(self):
        api_main = self.api.api("tool_lab.plasmid_predict")
        if self.option('main_id'):
            mainid =  self.option('main_id')
        else:
            mainid = api_main.add_plasmid_main()
        self.logger.info(mainid)
        if self.option("fasta").is_set:
            api_main.add_plasmid_detail(main_id=mainid, file_path=self.output_dir)
        else:
            for tool in self.tools:
                api_main.add_plasmid_detail(main_id=mainid, file_path=tool.output_dir)
        api_main.add_plasmid_main(main_id=mainid)
        for file in os.listdir(self.output_dir):
            if file.endswith(".log"):
                os.remove(self.output_dir + '/' + file)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "样品流程分析结果目录", 0],
            ["./*plasflow_predictions.tsv_plasmids.fasta", "fasta", "鉴定为plasmid的序列", 0],
            ["./*plasflow_predictions.tsv_chromosomes.fasta", "fasta", "鉴定为chromosome的序列件", 0],
            ["./*plasflow_predictions.tsv_unclassified.fasta", "fasta", "鉴定为unclassified的序列件", 0],
            ["./*plasflow_predictions.tsv", "tsv", "质粒鉴定统计表", 0],
            ["./*plasmid_blast.txt", "txt", "质粒与PLSDB数据库blast的m8格式文件件", 0],
            ["./*plasmid_plsdb.tsv", "tsv", "质粒注释详情表", 0],
        ])
        super(PlasmidPredictWorkflow, self).end()
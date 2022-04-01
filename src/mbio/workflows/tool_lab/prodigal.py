# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.tool_lab.common_function import fasta_dir_sort
import os
import subprocess
import gevent

class ProdigalWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProdigalWorkflow, self).__init__(wsheet_object)
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
        #self.prodigal = self.add_module('tool_lab.prodigal')
        self.tools = []

    def check_options(self):
        if (not self.option("fasta").is_set) and (not self.option("fasta_dir").is_set):
            raise OptionError("请传入序列文件或序列文件夹！", code="")

    def run(self):
        if self.option("fasta_dir").is_set:
            fasta_dir_sort(self.option("fasta_dir").prop['path'],self.work_dir + "/fasta_dir_sort")
        self.run_anno()
        super(ProdigalWorkflow, self).run()

    def run_anno(self):
        """
        prodial 注释
        :return:
        """
        if self.option("fasta").is_set:
            self.prodigal = self.add_module('tool_lab.prodigal')
            opts = {
                "genome": self.option("fasta"),
                "sample": self.option("fasta").prop['path'].split("/")[-1].split(".")[0],
            }
            self.prodigal.set_options(opts)
            self.prodigal.on('end', self.set_output)
            self.prodigal.run()
        else:
            files = os.listdir(os.path.join(self.work_dir, "fasta_dir_sort"))
            for file in files:
                self.prodigal = self.add_module('tool_lab.prodigal')
                opts= {
                    "genome": os.path.join(os.path.join(self.work_dir,"fasta_dir_sort"),file),
                    "sample": file.rstrip(".fasta")
                }
                self.prodigal.set_options(opts)
                self.tools.append(self.prodigal)
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
            link_dir(self.prodigal.output_dir + '/' +self.option("fasta").prop['path'].split("/")[-1].split(".")[0],self.output_dir + '/' +self.option("fasta").prop['path'].split("/")[-1].split(".")[0])
        else:
            for tool in self.tools:
                self.logger.info(tool)
                for file in os.listdir(tool.output_dir):
                    if os.path.exists(self.output_dir + "/" + file):
                        os.removedirs(self.output_dir + "/" + file)
                    link_dir(tool.output_dir + "/" + file, self.output_dir + "/" + file)
        self.set_db()

    def set_db(self):
        api_main = self.api.api("tool_lab.prodigal")
        if self.option('main_id'):
            mainid =  self.option('main_id')
        else:
            mainid = api_main.add_prodigal_main()
        self.logger.info(mainid)
        if self.option("fasta").is_set:
            api_main.add_prodigal_detail(main_id=mainid, file_path=self.output_dir)
        else:
            for tool in self.tools:
                if os.listdir(tool.output_dir):
                    api_main.add_prodigal_detail(main_id=mainid, file_path=tool.output_dir)
        api_main.add_prodigal_main(main_id=mainid)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "Prodigal编码基因预测流程分析结果目录", 0],
            ["./*", "dir", "编码基因预测结果文件夹", 0],
            ["./*/Prodigal_Gene prediction.faa", "faa", "编码基因预测蛋白文件", 0],
            ["./*/Prodigal_Gene prediction.fna", "fna", "编码基因预测核酸文件", 0],
            ["./*/Prodigal_Gene prediction.gff3", "gff", "编码基因预测信息gff文件", 0],
            ["./*/Prodigal_Gene prediction_detail.xls", "xls", "编码基因预测详情表", 0],
            ["./*/Prodigal_Gene prediction_stat.xls", "gff", "编码基因预测统计表", 0]
        ])
        super(ProdigalWorkflow, self).end()
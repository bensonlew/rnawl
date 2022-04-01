# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.tool_lab.common_function import fasta_dir_sort
import os
import subprocess
import gevent

class SrnaPredictWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SrnaPredictWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "file_style", "type": "string"},
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},
            {"name": "species", "type": "string", "default": "Bacteria"}, # 物种 Bacteria,Archaea,Viruses
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
        self.run_predict()
        super(SrnaPredictWorkflow, self).run()

    def run_predict(self):
        """
        sRNA预测
        :return:
        """
        if self.option("fasta").is_set:
            self.predict = self.add_tool("tool_lab.srna_predict")
            opts = {
                "fasta": self.option("fasta"),
                "sample_name": self.option("fasta").prop['path'].split("/")[-1].split(".")[0],
                "species" : self.option("species")
            }
            self.predict.set_options(opts)
            self.predict.on('end', self.set_output)
            self.predict.run()
        else:
            files = os.listdir(os.path.join(self.work_dir, "fasta_dir_sort"))
            for file in files:
                self.predict = self.add_tool("tool_lab.srna_predict")
                opts= {
                    "fasta": os.path.join(os.path.join(self.work_dir,"fasta_dir_sort"),file),
                    "sample_name": file.rstrip(".fasta")
                }
                self.predict.set_options(opts)
                self.tools.append(self.predict)
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
        stat_file = self.output_dir + "/all.hgene.xls"
        if self.option("fasta").is_set:
            link_dir(self.predict.output_dir ,self.output_dir)
        else:
            for tool in self.tools:
                self.logger.info(tool)
                for file in os.listdir(tool.output_dir):
                    if os.path.exists(self.output_dir + "/" + file):
                        os.remove(self.output_dir + "/" + file)
                    os.link(tool.output_dir + "/" + file,self.output_dir + "/" + file)
                os.remove(self.output_dir + "/length_stat.txt")
        self.set_db()

    def set_db(self):
        api_main = self.api.api("tool_lab.srna_predict")
        if self.option('main_id'):
            mainid =  self.option('main_id')
        else:
            mainid = api_main.add_srna_main()
        self.logger.info(mainid)
        if self.option("fasta").is_set:
            api_main.add_srna_detail(main_id=mainid, file_path=self.output_dir)
        else:
            for tool in self.tools:
                api_main.add_srna_detail(main_id=mainid, file_path=tool.output_dir)
        api_main.add_srna_main(main_id=mainid)
        if os.path.exists(self.output_dir + "/length_stat.txt"):
            os.remove(self.output_dir + "/length_stat.txt")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "sRNA预测结果目录", 0],
            ["./*_assembly.tblout", "tblout", "sRNA预测详情表", 0],
            ["./*_assembly.txt", "txt", "sRNA预测结果表", 0],
        ])
        super(SrnaPredictWorkflow, self).end()
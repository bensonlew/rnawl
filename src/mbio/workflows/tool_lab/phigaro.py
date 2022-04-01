# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.tool_lab.common_function import fasta_dir_sort
import os
import subprocess
import gevent

class PhigaroWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PhigaroWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "file_style", "type": "string"},
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},
            {"name": "species", "type": "string", "default": "Bacteria"}, # 物种 Bacteria,Archaea,Viruses,Mitochondria
            {"name": "srna_predict", "type": "string","default": "false"},
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
        super(PhigaroWorkflow, self).run()

    def run_anno(self):
        """
        phigaro 前噬菌体预测
        :return:
        """
        if self.option("fasta").is_set:
            self.phigaro = self.add_tool("tool_lab.phigaro")
            opts = {
                "fasta": self.option("fasta"),
                "sample_name": self.option("fasta").prop['path'].split("/")[-1].split(".")[0],
            }
            self.phigaro.set_options(opts)
            self.phigaro.on('end', self.set_output)
            self.phigaro.run()
        else:
            files = os.listdir(os.path.join(self.work_dir, "fasta_dir_sort"))
            for file in files:
                self.phigaro = self.add_tool("tool_lab.phigaro")
                opts= {
                    "fasta": os.path.join(os.path.join(self.work_dir,"fasta_dir_sort"),file),
                    "sample_name": file.rstrip(".fasta"),
                }
                self.phigaro.set_options(opts)
                self.tools.append(self.phigaro)
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
            link_dir(self.phigaro.output_dir + '/' +self.option("fasta").prop['path'].split("/")[-1].split(".")[0],self.output_dir + '/' +self.option("fasta").prop['path'].split("/")[-1].split(".")[0])
        else:
            for tool in self.tools:
                self.logger.info(tool)
                for file in os.listdir(tool.output_dir):
                    if os.path.exists(self.output_dir + "/" + file):
                        os.removedirs(self.output_dir + "/" + file)
                    link_dir(tool.output_dir + "/" + file,self.output_dir + "/" + file)
        self.set_db()

    def set_db(self):
        api_main = self.api.api("tool_lab.phigaro")
        if self.option('main_id'):
            main_id =  self.option('main_id')
        else:
            main_id = api_main.add_phigaro_main()
        if self.option("fasta").is_set:
            api_main.add_phigaro_detail(main_id=main_id, file_path=self.output_dir)
        else:
            for tool in self.tools:
                if os.listdir(tool.output_dir):
                    api_main.add_phigaro_detail(main_id=main_id, file_path=tool.output_dir)
        api_main.add_phigaro_main(main_id)
        if os.path.exists(self.output_dir + "/stat"):
            os.remove(self.output_dir + "/stat")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "Phigaro前噬菌体预测结果目录", 0],
            ["./*", "dir", "Phigaro前噬菌体预测详细结果", 0],
            ["./*/phigaro_detail.xls", "xls", "前噬菌体预测详情表", 0],
            ["./*/phigaro.gff3", "gff3", "前噬菌体预测gff文件", 0],
            ["./*/phigaro_stat.xls", "xls", "前噬菌体预测统计表", 0],
        ])
        super(PhigaroWorkflow, self).end()
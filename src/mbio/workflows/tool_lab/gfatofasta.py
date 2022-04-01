# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.tool_lab.common_function import gfa_dir_sort
import os
import re, shutil
import subprocess
import gevent

class GfatofastaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GfatofastaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "file_style", "type": "string"},
            {"name": "gfa", "type": "infile", "format": "tool_lab.no_empty"},
            {"name": "gfa_dir", "type": "infile", "format": "tool_lab.input_dir"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.tools = []

    def check_options(self):
        if (not self.option("gfa").is_set) and (not self.option("gfa_dir").is_set):
            raise OptionError("请传入GFA文件或GFA文件夹！")

    def run(self):
        if self.option("gfa_dir").is_set:
            gfa_dir_sort(self.option("gfa_dir").prop['path'],self.work_dir + "/gfa_dir_sort")
        self.run_awk()
        super(GfatofastaWorkflow, self).run()

    def run_awk(self):
        """
        目的片段查找
        :return:
        """
        if self.option("gfa").is_set:
            self.gfatofasta = self.add_tool("tool_lab.gfatofasta")
            opts = {
                "gfa_file": self.option("gfa"),
                "sample_name": os.path.basename(self.option("gfa").prop["path"]).rstrip(".gfa")
            }
            self.gfatofasta.set_options(opts)
            self.gfatofasta.on('end', self.set_output)
            self.gfatofasta.run()
        else:
            files = os.listdir(os.path.join(self.work_dir, "gfa_dir_sort"))
            for file in files:
                self.gfatofasta = self.add_tool("tool_lab.gfatofasta")
                opts= {
                    "gfa_file": os.path.join(os.path.join(self.work_dir,"gfa_dir_sort"),file),
                    "sample_name": file.rstrip(".gfa")
                }
                self.gfatofasta.set_options(opts)
                self.tools.append(self.gfatofasta)
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
        if self.option("gfa").is_set:
            link_dir(self.gfatofasta.output_dir, self.output_dir)
        else:
            for tool in self.tools:
                self.logger.info(tool)
                for file in os.listdir(tool.output_dir):
                    if os.path.exists(self.output_dir + "/" + file):
                        os.remove(self.output_dir + "/" + file)
                    os.link(tool.output_dir + "/" + file,self.output_dir + "/" + file)
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_main = self.api.api("tool_lab.common")
        api_main.add_main_table("gfatofasta", main_id=self.option('main_id'))
        self.logger.info("导表完成")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "GFA转化FASTA结果目录", 0],
            ["./*.fasta", "fasta", "gfa转化得到的fasta文件", 0],
        ])
        super(GfatofastaWorkflow, self).end()
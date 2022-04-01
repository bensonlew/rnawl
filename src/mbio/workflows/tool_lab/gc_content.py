# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
import os
import subprocess

class GcContentWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GcContentWorkflow, self).__init__(wsheet_object)
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

    def check_options(self):
        if (not self.option("fasta").is_set) and (not self.option("fasta_dir").is_set):
            raise OptionError("请传入序列文件或者序列文件夹！", code="")

    def run(self):
        self.run_gc()
        super(GcContentWorkflow, self).run()

    def run_gc(self):
        """
        GC含量统计
        :return:
        """
        if self.option("fasta").is_set:
            self.file_name = self.option("fasta").prop['path'].split("/")[-1].split(".")[0]
            self.gc = self.add_tool("tool_lab.gc_content")
            opts = {
                "fasta": self.option("fasta"),
                "sample_name": self.file_name
            }
            self.gc.set_options(opts)
            self.gc.on('end', self.set_output)
            self.gc.run()
        else:
            self.gc = self.add_module("tool_lab.gc_content")
            opts = {
                "fasta_dir": self.option("fasta_dir")
            }
            self.gc.set_options(opts)
            self.gc.on('end', self.set_output)
            self.gc.run()


    def set_db(self):
        api_main = self.api.api("tool_lab.gc_content")
        api_main.add_gc_detail(main_id=self.option('main_id'), file_path=self.output_dir,)
        self.end()

    def set_output(self):
        link_dir(self.gc.output_dir, self.output_dir)
        self.set_db()
        #self.end()

    def get_list(self):
        if os.path.exists(self.option("fasta_dir").prop['path'] + '/' + "list.txt"):
            self.sample_dict = {}
            with open(self.option("fasta_dir").prop['path'] + '/' + "list.txt") as t:
                for i in t.readlines()[1:]:
                    if i.split("\t") not in self.sample_dict.keys():
                        self.sample_dict[i.split("\t")[0]] = [i.split("\t")[1].strip()]
                    else:
                        self.sample_dict[i.split("\t")[0]].append(i.split("\t")[1].strip())
        else:
            raise OptionError("list文件不存在！", code="")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "序列GC含量结果目录", 0],
            ["./*all.result.xls", "xls", "序列统计表", 0],
            ["./*sample.detail.result.xls", "xls", "序列详情统计表", 0]
        ])
        super(GcContentWorkflow, self).end()
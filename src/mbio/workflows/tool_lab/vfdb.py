# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.tool_lab.common_function import fasta_dir_sort
import os
import re, shutil
import subprocess
import gevent

class VfdbWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VfdbWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "seq_style", "type": "string","default":"nucl"},
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
        self.run_vfdb()
        super(VfdbWorkflow, self).run()

    def run_vfdb(self):
        """
        QUAST基因组评估
        :return:
        """
        if self.option("fasta").is_set:
            if os.path.exists(os.path.join(self.work_dir, "fasta_dir_sort")):
                shutil.rmtree(os.path.join(self.work_dir, "fasta_dir_sort"))
            os.mkdir(os.path.join(self.work_dir, "fasta_dir_sort"))
            name = self.option("fasta").prop["path"].split("/")[-1]
            if name.endswith(".gz"):
                os.system("gzip -d {}".format(self.option("fasta").prop["path"]))
            os.link(self.option("fasta").prop["path"].strip(".gz"),self.work_dir+"/"+"fasta_dir_sort"+"/"+name.split(".")[0] + ".fasta")
            self.vfdb = self.add_module("tool_lab.vfdb")
            opts = {
                "sample": name.split(".")[0],
                "query": self.work_dir+"/"+"fasta_dir_sort"+"/"+name.split(".")[0] + ".fasta",
                "query_type" : self.option("seq_style")
            }
            self.vfdb.set_options(opts)
            self.vfdb.on('end', self.set_output)
            self.vfdb.run()

        if self.option("fasta_dir").is_set:
            for file in os.listdir(self.work_dir + "/fasta_dir_sort"):
                self.vfdb = self.add_module("tool_lab.vfdb")
                opts = {
                        "sample": file.rstrip(".fasta"),
                        "query": self.work_dir + "/fasta_dir_sort/" + file,
                        "query_type" : self.option("seq_style")
                    }
                self.vfdb.set_options(opts)
                self.tools.append(self.vfdb)
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
            link_dir(self.vfdb.output_dir,self.output_dir)
        else:
            for tool in self.tools:
                for dir in os.listdir(tool.output_dir):
                    link_dir(tool.output_dir+"/"+dir, self.output_dir+"/"+dir)
        self.set_db()

    def set_db(self):
        api_main = self.api.api("tool_lab.vfdb")
        if self.option('main_id'):
            main_id =  self.option('main_id')
        else:
            main_id = api_main.add_vfdb_main()
        for dir in os.listdir(self.output_dir):
            api_main.add_vfdb_detail(main_id,file_path=self.output_dir + "/" + dir)
        api_main.add_vfdb_main(main_id)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "VFDB数据库毒力基因预测结果目录", 0],
            ["./*/", "dir", "VFDB数据库毒力基因预测结果目录", 0],
            ["./*/*vfdb_align.xls", "xls", "毒力基因预测表", 0],
            ["./*/*vfdb_anno.xls", "xls", "毒力基因注释表", 0],
            ["./*/*vfdb_level.xls", "xls", "毒力基因分级统计表", 0],
        ])
        super(VfdbWorkflow, self).end()
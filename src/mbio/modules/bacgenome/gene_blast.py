# -*- coding: utf-8 -*-
# __author__ = 'ysh'
# last_modify:20190409

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
import unittest


class GeneBlastModule(Module):
    def __init__(self, work_id):
        super(GeneBlastModule, self).__init__(work_id)
        options = [
            {"name": "sample_dir", "type": "infile", "format": "sequence.fasta_dir"},
            {"name": "search_sequence", "type": "infile", "format": "sequence.fasta"},
            {"name": "method", "type": "string"},
            {"name": "evalue", "type": "string", "default": "1e-5"},
            {"name": "max_target_num", "type": "string", "default": 10},
            {"name": "w_size", "type": "string"},
        ]
        self.add_option(options)
        self.tool_list = []
        self.stat_list = []
        self.sample_faa = {}

    def check_options(self):
        if not self.option("sample_dir"):
            raise OptionError("必须设置输入样本序列文件夹")
        if not self.option("method"):
            raise OptionError("必须输入比对方法")
        if not self.option("w_size"):
            raise OptionError("必须输入w_size")
        return True

    def run_gene_blast(self):
        if self.option("method") in ["blastn","blastx","tblastn"]:
            query_type = "nucl"
        else:  #blastp
            query_type = "prot"

        if self.option("method") in ["blastn","tblastn"]:
            reference_type = "nucl"
        elif self.option("method") in ["blastx","blastp"]:
            reference_type = "prot"

        samples = os.listdir(self.option("sample_dir").prop["path"])
        self.ref_name = []
        self.file_path_map_sample = {}
        for eachsample in samples:
            filepath = os.path.join(self.option("sample_dir").prop["path"],eachsample)
            self.file_path_map_sample[filepath] = eachsample
            self.sample_faa[eachsample] = filepath
            self.blast_tool = self.add_tool("align.blast")
            opts = {
                'query': self.option("search_sequence"),
                'query_type': query_type,
                'outfmt': 6,
                "blast" : self.option("method"),
                "reference": filepath,
                "reference_type":reference_type,
                "evalue": self.option("evalue"),
                "num_alignment": self.option("max_target_num"),
                "w_size" : self.option("w_size"),
                "database" : "customer_mode"
            }
            self.ref_name.append(eachsample)
            self.blast_tool.set_options(opts)
            self.tool_list.append(self.blast_tool)
        if len(self.tool_list) > 1:
            self.on_rely(self.tool_list, self.run_blast_stat)
            self.logger.info(self.tool_list)
        else:
            self.tool_list[0].on('end', self.run_blast_stat)
        for tool in self.tool_list:
            tool.run()

    def run_blast_stat(self):
        for index,tool in enumerate(self.tool_list):
            #sample = self.ref_name[index]
            sample = self.file_path_map_sample[tool.option('reference').prop['path']]
            eachtable = tool.option("outtable").prop["path"]
            self.stat_tool = self.add_tool("bacgenome.blast_detail")
            opts = {
                "stat_file": eachtable,
                #"faa_file": self.sample_faa[sample],
                "faa_file": tool.option('reference').prop['path'],
                "sample": sample
            }
            self.stat_tool.set_options(opts)
            self.stat_list.append(self.stat_tool)
        if len(self.stat_list) > 1:
            self.on_rely(self.tool_list, self.set_output)
        else:
            self.stat_list[0].on('end', self.set_output)
        for tool in self.stat_list:
            tool.run()

    def set_output(self):
        for each in self.stat_list:
            oldfiles = os.listdir(each.output_dir)
            for f in oldfiles:
                oldfile = each.output_dir + '/' + f
                newfile = os.path.join(self.output_dir,f)
                if os.path.exists(newfile):
                    os.remove(newfile)
                os.link(oldfile, newfile)
        self.end()

    def run(self):
        super(GeneBlastModule, self).run()
        self.run_gene_blast()

    def end(self):
        super(GeneBlastModule, self).end()



# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# last_modify:20181114

from biocluster.module import Module
import os,shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metabolome.common import link_file,link_dir


class SecAnnoModule(Module):
    def __init__(self, work_id):
        super(SecAnnoModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "lines", "type": "int", "default": 20000},
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "tax_table", "type": "infile", "format": "sequence.profile_table"},  # 分泌系统预测做一种注释的预测
            {"name": "tax_table_list", "type": "string"},  # ttss分析做三种不同注释的预测
            {"name": "anno_tool", "type": "string", "default": "signalp"},
            {"name": "signalp_out_format", "type": "string"},
            {"name": "signalp_type", "type": "string"},
            {"name": "sec_result_dir", "type": "outfile", "format": "annotation.mg_anno_dir"}
        ]
        self.add_option(options)
        self.split_fasta = self.add_tool("sequence.split_fasta")
        self.isolation_type = list()
        self.predict_tools = list()
        self.tax_tables = list()
        self.sec_anno_tool = ""
        self.sec_prof_tool = ""
        self.predict_result_path = os.path.join(self.work_dir, "predict_tmp")

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("must input option query")
        if self.option("lines") > 20000:
            raise OptionError("Max 20000 sequences are allowed")
        if self.option("signalp_type"):
            self.isolation_type = self.option("signalp_type").split(',')
            for predict_type in self.isolation_type:
                if predict_type not in ["gram-", "gram+", "euk"]:
                    raise OptionError("isolation type must be gram-,gram+ or euk")
        if self.option("tax_table").is_set:
            self.tax_tables.append(self.option("tax_table"))
        elif self.option("tax_table_list"):
            self.tax_tables = self.option("tax_table_list").split(",")
        else:
            raise OptionError("must input tax_table or tax_table_list")

    def run_split_fasta(self):
        self.split_fasta.set_options({
            "fasta": self.option("query"),
            "lines": self.option("lines")
        })
        self.split_fasta.on("end", self.predict)
        self.split_fasta.run()

    def predict(self):
        if self.option("anno_tool") == "signalp":
            self.run_signalp()
        else:
            self.end()

    def run_signalp(self):
        for f in os.listdir(self.split_fasta.output_dir):
            file_path = os.path.join(self.split_fasta.output_dir, f)
            for itype in self.isolation_type:
                signalp_tool = self.add_tool("predict.signalp")
                signalp_tool.set_options({
                    "query": file_path,
                    "type": itype,
                    "out_format": self.option("signalp_out_format"),
                    "mem": 10
                })
                self.predict_tools.append(signalp_tool)
        if len(self.predict_tools) > 1:
            self.on_rely(self.predict_tools, self.end)
        else:
            self.predict_tools[0].on('end', self.end)
        for tool in self.predict_tools:
            tool.run()

    def signalp_anno(self):
        input_path = self.predict_result_path
        1 if os.path.exists(input_path) else os.mkdir(input_path)
        for tool in self.predict_tools:
            for f in os.list_dir(tool.output_dir):
                file_path = os.path.join(tool.output_dir, f)
                new_path = os.path.join(input_path, f)
                link_file(file_path, new_path)
        self.sec_anno_tool.set_options({
            "predict_dir": input_path,
            "gene_nr_anno": self.option("tax_table")
        })
        self.sec_anno_tool.on('end', self.signalp_anno_stat)
        self.sec_anno_tool.run()

    def signalp_anno_stat(self):
        self.sec_prof_tool.set_options({
            'signalp_anno_table': self.sec_anno_tool.option("signalp_anno_result"),
            'reads_profile_table': self.option('reads_profile_table')
        })
        self.sec_prof_tool.on("end", self.set_output)
        self.sec_prof_tool.run()

    def set_output(self):
        self.option("sec_result_dir", self.output_dir)
        if self.option("anno_tool") == "signalp":
            link_dir(self.sec_prof_tool, self.output_dir)
        self.end()

    def run(self):
        super(SecAnnoModule, self).run()
        self.run_split_fasta()
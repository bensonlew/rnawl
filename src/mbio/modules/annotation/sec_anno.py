# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# last_modify:20181114

from biocluster.module import Module
import os, shutil, glob
from biocluster.core.exceptions import OptionError
from mbio.packages.metabolome.common import link_file,link_dir


class SecAnnoModule(Module):
    def __init__(self, work_id):
        super(SecAnnoModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "lines", "type": "int", "default": 200000},
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
        self.isolation_type = list()  # signalp預測模型gram+ gram- euk
        self.predict_tools = list()  # signalp預測需要三個不同的模型
        self.ttss_tool = ""  # ttss預測只用一個模型
        self.tax_tables = list()
        self.sec_stat_tool = self.add_tool("annotation.sec_anno_stat")  # signalp預測只用一個nr注釋表
        self.sec_stat_tools = list()  # ttss預測需要三個不同的nr注釋表
        self.predict_result_path = os.path.join(self.work_dir, "predict_tmp")

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("must input option query", code="21202401")
        if self.option("anno_tool") == "signalp" and self.option("lines") > 200000:
            raise OptionError("Max 200000 sequences are allowed in signalp", code="21202402")
        if self.option("signalp_type"):
            self.isolation_type = self.option("signalp_type").split(',')
            for predict_type in self.isolation_type:
                if predict_type not in ["gram-", "gram+", "euk"]:
                    raise OptionError("isolation type must be gram-,gram+ or euk", code="21202403")
        if self.option("tax_table").is_set:
            self.tax_tables.append(self.option("tax_table"))
        elif self.option("tax_table_list"):
            self.tax_tables = self.option("tax_table_list").split(",")
        else:
            raise OptionError("must input tax_table or tax_table_list", code="21202404")

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
            self.run_ttss()

    def run_signalp(self):
        for itype in self.isolation_type:
            signalp_tool = self.add_tool("predict.signalp")
            signalp_tool.set_options({
                "query_dir": self.split_fasta.output_dir,
                "type": itype,
                "out_format": self.option("signalp_out_format"),
                "mem": 80  # 增加到50 by ghd @20190305 by qingchen.zhang@20200415 增加到80G
            })
            self.predict_tools.append(signalp_tool)
        if len(self.predict_tools) > 1:
            self.on_rely(self.predict_tools, self.signalp_stat)
        else:
            self.predict_tools[0].on('end', self.signalp_stat)
        for tool in self.predict_tools:
            tool.run()

    def signalp_stat(self):
        input_path = self.predict_result_path
        1 if os.path.exists(input_path) else os.mkdir(input_path)
        for tool in self.predict_tools:
            tool_result = tool.option('result').path
            f = os.path.basename(tool_result)
            link_file(tool_result, os.path.join(input_path, f))
        self.sec_stat_tool.set_options({
            "predict_dir": input_path,
            "gene_nr_anno": self.option("tax_table"),
            "reads_profile_table": self.option("reads_profile_table")
        })
        self.sec_stat_tool.on('end', self.set_output)
        self.sec_stat_tool.run()

    def run_ttss(self):
        self.ttss_tool = self.add_tool("predict.ttss")
        self.ttss_tool.set_options({
            # "query_dir": self.split_fasta.output_dir,
            "query": self.option("query"),
            "mem": 45
        })
        self.ttss_tool.on('end', self.run_ttss_stat)
        self.ttss_tool.run()

    def run_ttss_stat(self):
        self.logger.info("DDDDDDDEBUG: NOW CHECK gene_nr_anno")
        for table in self.tax_tables:
            stat_tool = self.add_tool("annotation.sec_anno_stat")
            self.logger.info(table)
            stat_tool.set_options({
                "predict_dir": self.ttss_tool.output_dir,
                "gene_nr_anno": table,
                "nr_type": table.split('/')[-2],
                "reads_profile_table": self.option("reads_profile_table")
            })
            self.sec_stat_tools.append(stat_tool)
        self.logger.info("DDDDDDDDEBUG: debug end")
        if len(self.sec_stat_tools) > 1:
            self.on_rely(self.sec_stat_tools, self.set_output)
        else:
            self.sec_stat_tools[0].on('end', self.set_output)
        for tool in self.sec_stat_tools:
            tool.run()

    # def signalp_anno_stat(self):
    #     self.sec_prof_tool.set_options({
    #         'signalp_anno_table': self.sec_anno_tool.option("signalp_anno_result"),
    #         'reads_profile_table': self.option('reads_profile_table')
    #     })
    #     self.sec_prof_tool.on("end", self.set_output)
    #     self.sec_prof_tool.run()

    def set_output(self):
        self.option("sec_result_dir", self.output_dir)
        if self.option("anno_tool") == "signalp":
            link_dir(self.sec_stat_tool.output_dir, self.output_dir)
            link_dir(self.predict_result_path, self.output_dir)
        elif self.option("anno_tool") == "EffectiveT3":
            link_dir(self.ttss_tool.output_dir, self.output_dir)
            for tool in self.sec_stat_tools:
                # old_file1 = os.path.join(tool.output_dir, "ttss_summary.txt")
                # old_file2 = os.path.join(tool.output_dir, "ttss_fisher.txt")
                # tax_type = tool.option('gene_nr_anno').path.split('/')[-2]  # nr nr_lca nr_deunclassified
                # new_name1 = "ttss_%s_summary.txt" % tax_type
                # new_name2 = "ttss_%s_fisher.txt" % tax_type
                # link_file(old_file1, os.path.join(self.output_dir, new_name1))
                # link_file(old_file2, os.path.join(self.output_dir, new_name2))
                link_dir(tool.output_dir, self.output_dir)
        self.end()

    def run(self):
        super(SecAnnoModule, self).run()
        # self.run_split_fasta()
        if self.option("anno_tool") == "signalp":
            self.run_split_fasta()
        elif self.option("anno_tool") == "EffectiveT3":
            self.run_ttss()  # TTSS拆分序列后增加了运行内存，而实际运行速度没有明显增加

# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2018/12/25'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import Fasta
from mbio.packages.align.pocp import Pocp
from mbio.packages.metagbin.common_function import link_tax_dir

class BinTaxModule(Module):
    """
    binning物种注释
    """

    def __init__(self, work_id):
        super(BinTaxModule, self).__init__(work_id)
        option = [
            {"name": "query", "type": "infile", "format": "sequence.fasta,sequence.fasta_dir"},  # str/int/float/bool...
            {"name": "ref", "type": "infile", "format": "sequence.fasta,sequence.fasta_dir"},  # infile
            {"name": "method", "type": "string", "default": "blasr"}, # blasr/pocp/ani
            {"name": "sample", "type": "string"},
            {"name": "blasr_db", "type": "string"}, # 相对路径，去掉后缀
            {"name": "task_id", "type": "string"}
        ]
        self.add_option(option)
        self.query2ref = self.add_tool('align.blast')
        self.ref2query = self.add_tool('align.blast')
        self.blasr = self.add_tool("align.blasr")
        self.ani = self.add_tool("align.ani")
        self.stat = self.add_tool("metagbin.tax_stat")

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("query").is_set:
            raise OptionError("请输入分析基因集")
        if self.option("method") == "blasr":
            # if self.option("blasr_db") not in ["silva123", "silva128", "silva132", "greengene", "nt"]:
            #    raise OptionError("blasr数据库不在范围内")
            if not self.option("blasr_db"):
                raise OptionError("必须输入blasr数据库")
        elif self.option("method") in ["pocp", "ani"]:
            if not self.option("ref").is_set:
                raise OptionError("请输入参考基因集")
        else:
            raise OptionError("method分析方法不在范围内")
        return True

    def run(self):
        super(BinTaxModule, self).run()
        self.stat.on("end", self.set_output)
        if self.option("method") == "blasr":
            self.run_blasr()
        elif self.option("method") == "pocp":
            self.run_pocp()
        elif self.option("method") == "ani":
            self.run_ani()

    def run_blasr(self):
        # 需要获取到最长的一条序列
        new_fasta = os.path.join(self.work_dir, "query.fasta")
        obj = Fasta(self.option("query").path)
        obj.parse()
        obj.export_longest_seq(new_fasta)
        opts = {
            "query": new_fasta,
            "reference": self.option('blasr_db'),
            "outfmt": 4,
            "num_threads": 10,
        }
        if "nt" in self.option("blasr_db"):
            opts["memory"] = 40
        else:
            opts["memory"] = 25
        self.blasr.set_options(opts)
        self.blasr.on("end", self.tax_stat, "blasr")
        self.blasr.run()

    def run_pocp(self):
        self.on_rely([self.query2ref, self.ref2query], self.tax_stat, "pocp")
        # self.on_rely([self.query2ref, self.ref2query], self.pocp_stat)
        self.run_blastp("process1")
        self.run_blastp("process2")

    def generate_list(self, dir_path, outfile):
        with open(outfile, "w") as file:
            for name in os.listdir(dir_path):
                rel_path = os.path.join(dir_path, name)
                file.write(rel_path + "\n")

    def run_ani(self):
        opts = {}
        self.logger.info(self.option("query").format)
        if self.option("query").format == "sequence.fasta":
            opts["query"] = self.option("query")
        else:
            query_list = os.path.join(self.work_dir, "query_list")
            self.generate_list(self.option("query").path, query_list)
            opts["query_list"] = query_list
        if self.option("ref").format == "sequence.fasta":
            opts["ref"] = self.option("ref")
        else:
            ref_list = os.path.join(self.work_dir, "ref_list")
            self.generate_list(self.option("ref").path, ref_list)
            opts["ref_list"] = ref_list
        self.ani.set_options(opts)
        self.ani.on("end", self.tax_stat, "ani")
        self.ani.run()

    def run_blastp(self, process):
        opts = {
            'query_type': "prot",
            "database": "customer_mode",
            "outfmt": 5,
            "blast": "blastp",
            "reference_type": "prot"
        }
        if process == "process1":
            opts['query'] = self.option("query")
            opts['reference'] = self.option("ref")
            self.query2ref.set_options(opts)
            self.query2ref.run()
        elif process == "process2":
            opts['query'] = self.option("ref")
            opts['reference'] = self.option("query")
            self.ref2query.set_options(opts)
            self.ref2query.run()
        else:
            self.set_error("no such process named %s" % process)

    def tax_stat(self, event):
        self.logger.info("in stat")
        self.logger.info(event)
        ana = event["data"]
        opts = {"analysis": ana, "sample": self.option("sample"), "task_id": self.option("task_id")}
        if ana == "blasr":
            if os.path.getsize(self.blasr.option("outtable").prop["path"]) == 0:
                self.end()
            opts["blasr_table"] = self.blasr.option("outtable")
            opts["tax_file"] = self.option('blasr_db')
        elif ana == "pocp":
            opts.update({
                "query": self.option('query'),
                "ref": self.option("ref"),
                "query_blast": self.query2ref.option("outxml"),
                "ref_blast": self.ref2query.option("outxml")
            })
        elif ana == "ani":
            if os.path.getsize(self.ani.option("outfile").prop["path"]) == 0:
                self.end()
            self.logger.info("important: %s" % self.ani.option("outfile").path)
            # opts["query"] = self.option('query')
            # opts["ref"] = self.option('ref')
            opts["ani_table"] = self.ani.option("outfile")  # 加
        self.logger.info(opts)
        self.stat.set_options(opts)
        self.stat.run()

    def pocp_stat(self):
        result = os.path.join(self.work_dir, "pocp_result.xls")
        #self.logger.debug(self.query2ref.option("outtable").path)
        #self.logger.debug(self.ref2query.option("outtable").path)
        self.logger.debug(self.query2ref.option("outxml").path)
        self.logger.debug(self.ref2query.option("outxml").path)
        self.logger.debug(result)
        self.end()
        pocp_obj = Pocp()
        pocp_obj.parse_fa(self.option('query').path) # 允许输入为fasta,或者字符，处理到公共变量中
        pocp_obj.parse_fa(self.option('ref').path)
        # pocp_obj.parse_table(self.query2ref.option("outtable").path, "m6") # 允许输入为m5或m6格式的表，如果是m5不需要parse_fa
        # pocp_obj.parse_table(self.ref2query.option("outtable").path, "m6")
        pocp_obj.parse_table(self.query2ref.option("outxml").path, "m5")
        pocp_obj.parse_table(self.ref2query.option("outxml").path, "m5")
        self.logger.debug(pocp_obj.export(result))
        self.set_output()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        link_tax_dir(self.stat.output_dir, self.output_dir)
        self.logger.info("设置注释结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BinTaxModule, self).end()

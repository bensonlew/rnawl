# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 20180517

import os
import re
from Bio import SeqIO
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class SsrCompareWorkflow(Workflow):
    """
    交互分析：样本基因组SSR比较分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SsrCompareWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "query_fa1", "type": "infile", "format": "sequence.fasta"},  # 要比对的序列
            {"name": "query_fa2", "type": "infile", "format": "sequence.fasta"},  # 要比对的序列
            {"name": "ssr_result1", "type": "infile", 'format': 'bsa.vcf'},  # ssr结果
            {"name": "ssr_result2", "type": "infile", 'format': 'bsa.vcf'},  # ssr结果
            {"name": "dbname_nsq", "type": "string"},  # 比对库的路径
            {"name": "outfmt", "type": "string", "default": "6"},  # 输出格式 tab格式是6
            {"name": "num_threads", "type": "int", "default": 8},
            {"name": "evalue", "type": "float", "default": 1e-10},  # 增加参数evalue、num_alignments modified by zengjing 20180508
            {"name": "num_alignments", "type": "int", "default": 5},
            {"name": "seq_length", "type": "int", "default": 1000},  # 序列长度，用于进行筛选序列
            {"name": "dbtype", "type": "string", 'default': "nucl"},
            {"name": "parse_seqids", "type": "bool", "default": True},  # 是否在makeblastdb的时候增加参数-parse_seqids
            {"name": "is_same", "type": "bool"},  # 样本1和样本2之间相同或不同的SSR位点
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "project_type", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.make_blastdb = False
        self.times = 0
        self.query_dict = {}
        self.sample_id = ""

    def check_options(self):
        if not self.option("query_fa1").is_set:
            raise OptionError("请设置要比对的序列query_fa1", code="14501101")
        if not self.option("query_fa2").is_set:
            raise OptionError("请设置要比对的序列query_fa2", code="14501102")
        if not self.option("ssr_result1"):
            raise OptionError("请设置ssr_result1", code="14501103")
        if not self.option("ssr_result2"):
            raise OptionError("请设置ssr_result2", code="14501104")
        if not self.option("dbname_nsq"):
            raise OptionError("请设置比对库的路径", code="14501105")
        if self.option("is_same") not in [True, False]:
            raise OptionError("请设置样本间相同或者不同SSR", code="14501106")

    def get_long_fasta(self, query_fa):
        """
        筛选长度大于等于seq_length的序列
        """
        seq_length = self.option("seq_length")
        name = os.path.basename(query_fa).split(".")[0]
        new_query_fa = os.path.join(self.work_dir, name + ".fa")
        num = 0
        with open(new_query_fa, "w") as w:
            for seq_record in SeqIO.parse(query_fa, "fasta"):
                if len(seq_record.seq) >= seq_length:
                    num += 1
                    w.write('>{}\n{}\n'.format(seq_record.id, seq_record.seq))
        if name not in self.query_dict.keys():
            self.query_dict[name] = {"fa": new_query_fa, "ssr": self.option("ssr_result1").prop['path']}
        return num, name, new_query_fa

    def run_blastn(self):
        """
        进行blastn
        """
        self.logger.info("开始进行blastn")
        for name in self.query_dict.keys():
            query_fa = self.query_dict[name]["fa"]
            lines = int(self.query_dict[name]["num"] / 15)
            options = {
                "query_fa": query_fa,
                "lines": lines,
                "dbname_nsq": self.option("dbname_nsq"),
                "outfmt": self.option("outfmt"),
                "num_threads": self.option("num_threads"),
                "evalue": self.option("evalue"),
                "num_alignments": self.option("num_alignments"),
                "sample_id": name
            }
            self.blastn = self.add_module("wgs.blastn")
            self.blastn.set_options(options)
            self.blastn.on("end", self.run_blastn_out, name)
            self.blastn.run()

    def run_blastn_out(self, event):
        """
        对blastn的结果进行处理
        """
        self.logger.info("开始进行blastn_out")
        obj = event["bind_object"]
        name = event["data"]
        blast_out = os.path.join(obj.output_dir, name + ".blast.out.xls")
        query_fa = self.query_dict[name]["fa"]
        ssr_result = self.query_dict[name]["ssr"]
        options = {
            "blast_out": blast_out,
            "scafseq_out": query_fa,
            "ssr_result": ssr_result
        }
        self.blast_out = self.add_tool("wgs.blast_out")
        self.blast_out.set_options(options)
        self.blast_out.on("end", self.run_re_blast, name)
        self.blast_out.run()

    def run_re_blast(self, event):
        """
        将没有比对上数据库的序列提出来
        """
        obj = event["bind_object"]
        name = event["data"]
        result_out = os.path.join(obj.output_dir, name + ".result-2.out")
        query_fa = self.query_dict[name]["fa"]
        self.query_dict[name]["result1"] = os.path.join(obj.output_dir, name + ".result-1.out")
        self.query_dict[name]["result2"] = result_out
        options = {
            "result_out": result_out,
            "scafseq_out": query_fa
        }
        self.re_blast = self.add_tool("wgs.re_blast")
        self.re_blast.set_options(options)
        if not self.make_blastdb:
            self.make_blastdb = True
            self.re_blast.on("end", self.run_makeblastdb, name)
        else:
            self.re_blast.on("end", self.run_nomatch_blastn, "nomatchseq:" + name)
        self.re_blast.run()

    def run_makeblastdb(self, event):
        """
        建数据库
        """
        obj = event["bind_object"]
        name = event["data"]
        pop_fa = os.path.join(obj.output_dir, name + ".nomatch.seq")
        options = {
            "pop_fa": pop_fa,
            "dbtype": self.option("dbtype"),
            "parse_seqids": self.option("parse_seqids")
        }
        self.makeblastdb = self.add_tool("wgs.makeblastdb")
        self.makeblastdb.set_options(options)
        self.makeblastdb.on("end", self.run_nomatch_blastn, "nomatchdb:" + name)
        self.makeblastdb.run()

    def run_nomatch_blastn(self, event):
        """
        进行blastn:将样本1没比对上ref的序列blastn到样本2没比对上ref的序列将的数据库
        """
        obj = event["bind_object"]
        m = re.search("nomatchdb:(.+)", event["data"])
        n = re.search("nomatchseq:(.+)", event["data"])
        if m:
            name = m.group(1)
            self.nomatch_db = os.path.join(obj.output_dir, "dbname")
        if n:
            self.sample_id = n.group(1)
            self.nomatch_seq = os.path.join(obj.output_dir, self.sample_id + ".nomatch.seq")
        self.times += 1
        if self.times == 2:
            with open(self.nomatch_seq, "r")as f:
                lines = f.readlines()
                lines = int(len(lines) / 2 / 10)
            options = {
                "query_fa": self.nomatch_seq,
                "lines": lines,
                "dbname_nsq": self.nomatch_db,
                "outfmt": self.option("outfmt"),
                "num_threads": self.option("num_threads"),
                "evalue": self.option("evalue"),
                "num_alignments": self.option("num_alignments"),
                "sample_id": self.sample_id
            }
            self.blastn = self.add_module("wgs.blastn")
            self.blastn.set_options(options)
            self.blastn.on("end", self.run_nomatch_blastn_out, self.sample_id)
            self.blastn.run()

    def run_nomatch_blastn_out(self, event):
        """
        对没比对上的序列blastn的结果进行处理
        """
        obj = event["bind_object"]
        name = event["data"]
        blast_out = os.path.join(obj.output_dir, name + ".blast.out.xls")
        ssr_result = self.query_dict[name]["result2"]
        options = {
            "blast_out": blast_out,
            "scafseq_out": self.nomatch_seq,
            "ssr_result": ssr_result
        }
        self.blast_out = self.add_tool("wgs.blast_out")
        self.blast_out.set_options(options)
        self.blast_out.on("end", self.cat_files, name)
        self.blast_out.run()

    def cat_files(self, event):
        """
        将文件cat到一起
        """
        self.logger.info(self.query_dict)
        obj = event["bind_object"]
        name = event["data"]
        no_result1 = os.path.join(obj.output_dir, name + ".result-1.out")
        result1 = self.query_dict[name]["result1"]
        self.sample1_result = os.path.join(self.work_dir, name + ".ref.result.xls")
        os.system("cat {} {} > {}".format(result1, no_result1, self.sample1_result))
        for s in self.query_dict.keys():
            if s != name:
                self.sample2_result = os.path.join(self.work_dir, s + ".ref.result.xls")
                os.system("cat {} {} > {}".format(self.query_dict[s]["result1"], self.query_dict[s]["result2"], self.sample2_result))
        self.run_ssr_compare_stat()

    def run_ssr_compare_stat(self):
        """
        进行ssr统计
        """
        options = {
            "sample1_result": self.sample1_result,
            "sample2_result": self.sample2_result,
            "is_same": self.option("is_same")
        }
        self.ssr_compare_stat = self.add_tool("wgs.ssr_compare_stat")
        self.ssr_compare_stat.set_options(options)
        self.ssr_compare_stat.on("end", self.set_output)
        self.ssr_compare_stat.run()

    def set_output(self):
        for f in os.listdir(self.ssr_compare_stat.output_dir):
            new = os.path.join(self.output_dir, f)
            if os.path.exists(new):
                os.remove(new)
            os.link(os.path.join(self.ssr_compare_stat.output_dir, f), new)
        self.set_db()
        self.end()

    def set_db(self):
        """
        将结果保存到mongo
        """
        self.logger.info("将结果保存到mongo")
        api_ssr_compare = self.api.api("wgs.ssr_compare")
        if self.option("project_type"):
            api_ssr_compare._project_type = self.option("project_type")
        ssr_compare_id = self.option("main_id")
        ssr_stat = os.path.join(self.output_dir, "ssr_stat.xls")
        ssr_detail = os.path.join(self.output_dir, "ssr_detail.xls")
        api_ssr_compare.add_sg_ssr_compare_stat(ssr_compare_id, ssr_stat)
        # api_ssr_compare.add_sg_ssr_compare_detail(ssr_compare_id, ssr_detail)
        download_file = self._sheet.output.rstrip('/')
        api_ssr_compare.add_sg_ssr_compare_detail_new(ssr_compare_id, ssr_detail, download_file)

    def run(self):
        num1, name1, new_query_fa1 = self.get_long_fasta(self.option("query_fa1").prop["path"])
        self.query_dict[name1] = {"fa": new_query_fa1, "ssr": self.option("ssr_result1").prop['path'], "num": num1}
        num2, name2, new_query_fa2 = self.get_long_fasta(self.option("query_fa2").prop["path"])
        self.query_dict[name2] = {"fa": new_query_fa2, "ssr": self.option("ssr_result2").prop['path'], "num": num2}
        self.logger.info(self.query_dict)
        self.run_blastn()
        super(SsrCompareWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SsrCompareWorkflow, self).end()

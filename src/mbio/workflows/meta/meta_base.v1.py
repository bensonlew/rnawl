# -*- coding: utf-8 -*-
# __author__ = 'yuguo'

"""多样性基础分析"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class MetaBaseWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        """
        self._sheet = wsheet_object
        super(MetaBaseWorkflow, self).__init__(wsheet_object.id)
        options = [
            {'name': 'fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 输入fasta文件，序列名称格式为'>sampleID_seqID'.
            {'name': 'identity', 'type': 'float', 'default': 0.97},  # 相似性值，范围0-1.
            {'name': 'otu_table', 'type': 'outfile', 'format': 'meta.otu.otu_table'},  # 输出结果otu表
            {'name': 'otu_rep', 'type': 'outfile', 'format': 'sequence.fasta'},  # 输出结果otu代表序列
            {'name': 'otu_seqids', 'type': 'outfile', 'format': 'meta.otu.otu_seqids'},  # 输出结果otu中包含序列列表
            {'name': 'otu_biom', 'type': 'outfile', 'format': 'meta.otu.biom'},  # 输出结果biom格式otu表
            {'name': 'revcomp', 'type': 'bool'},  # 序列是否翻转
            {'name': 'confidence', 'type': 'float', 'default': 0.7},  # 置信度值
            {"name": "customer_mode", "type": "bool", "default": False},  # customer 自定义数据库
            {'name': 'database', 'type': 'string'},  # 数据库选择
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 参考fasta序列
            {'name': 'ref_taxon', 'type': 'infile', 'format': 'taxon.seq_taxon'},  # 参考taxon文件
            {'name': 'taxon_file', 'type': 'outfile', 'format': 'taxon.seq_taxon'},  # 输出序列的分类信息文件
            {"name": "estimate_indices", "type": "string", "default": "ace-chao-shannon-simpson"},  # 指数类型
            # {"name": "estimators", "type": "outfile", "format": "meta.alpha_diversity.estimators"},  # 输出结果
            {"name": "rarefy_indices", "type": "string", "default": "chao-shannon"},  # 指数类型
            {"name": "random_number", "type": "int", "default": 100},  # 随机取样数
            # {"name": "rarefaction", "type": "outfile", "format": "meta.alpha_diversity.rarefaction_dir"},  # 输出结果
            {'name': 'otu_taxon_biom', 'type': 'outfile', 'format': 'meta.otu.biom'},  # 输出的biom文件
            {'name': 'otu_taxon_table', 'type': 'outfile', 'format': 'meta.otu.otu_table'},  # 输出的biom文件
            {'name': 'otu_taxon_dir', 'type': 'outfile', 'format': 'meta.otu.tax_summary_abs_dir'}  # 输出的otu_taxon_dir文件夹
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.otu = self.add_tool("meta.otu.usearch_otu")
        self.tax = self.add_tool("taxon.qiime_assign")
        self.stat = self.add_tool("meta.otu.otu_taxon_stat")
        self.est = self.add_tool("meta.alpha_diversity.estimators")
        self.rarefy = self.add_tool("meta.alpha_diversity.rarefaction")

    def check_options(self):
        """
        检查参数设置
        """
        if not self.option("fasta").is_set:
            raise OptionError("必须设置输入fasta文件.")
        if self.option("identity") < 0 or self.option("identity") > 1:
            raise OptionError("identity值必须在0-1范围内.")
        if self.option("revcomp") not in [True, False]:
            raise OptionError("必须设置序列是否翻转")
        if self.option("customer_mode"):
            if not self.option("ref_fasta").is_set or not self.option("ref_taxon").is_set:
                raise OptionError("数据库自定义模式必须设置参考fasta序列和参考taxon文件")
        else:
            if self.option("database") not in ['silva119/16s_bacteria', 'silva119/16s_archaea', 'silva119/18s_eukaryota', 'unite6.0/its_fungi', 'fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS', 'fgr/nifH', 'fgr/pmoA', 'fgr/mmoX']:
                raise OptionError("数据库{}不被支持".fomat(self.option("database")))
        return True

    def run_otu(self):
        self.otu.set_options({
            "fasta": self.option("fasta"),
            "identity": self.option("identity")
            })

        self.on_rely(self.otu, self.run_taxon)
        self.on_rely(self.otu, self.run_estimate)
        self.on_rely(self.otu, self.run_rarefy)
        self.otu.on("end", self.set_output, "otu")
        self.otu.run()

    def run_taxon(self, relyobj):
        self.tax.set_options({
            "fasta": relyobj.rely[0].option("otu_rep"),
            "revcomp": self.option("revcomp"),
            "confidence": self.option("confidence"),
            "customer_mode": self.option("customer_mode"),
            "database": self.option("database")}
            )
        if self.option("customer_mode"):
            self.tax.set_options({
                "ref_fasta": self.option("ref_fasta"),
                "ref_taxon": self.option("ref_taxon")
            })
        self.on_rely(self.tax, self.run_stat)
        self.tax.on("end", self.set_output, "tax")
        self.tax.run()

    def run_stat(self, relyobj):
        self.stat.set_options({
            "otu_seqids": self.option("otu_seqids"),
            "taxon_file": self.option("taxon_file")
            })
        self.stat.on("end", self.set_output, "stat")
        self.stat.run()

    def run_estimate(self, relyobj):
        # self.est = self.add_tool("meta.alpha_diversity.estimators")
        self.est.set_options({
            "otutable": relyobj.rely[0].option("otu_table"),
            "indices": self.option("estimate_indices")
            })
        # self.est.on("end", self.set_output, "est")
        self.est.run()

    def run_rarefy(self, relyobj):
        # self.rarefy = self.add_tool("meta.alpha_diversity.rarefaction")
        self.rarefy.set_options({
            "otutable": relyobj.rely[0].option("otu_table"),
            "indices": self.option("rarefy_indices"),
            "random_number": self.option("random_number")
            })
        self.rarefy.run()

    def set_output(self, event):
        obj = event["bind_object"]
        if event['data'] is "otu":
            self.option("otu_table", obj.option("otu_table"))
            self.option("otu_rep", obj.option("otu_rep"))
            self.option("otu_seqids", obj.option("otu_seqids"))
            self.option("otu_biom", obj.option("otu_biom"))
        if event['data'] is "tax":
            self.option("taxon_file", obj.option("taxon_file"))
        if event['data'] is "stat":
            self.option("otu_taxon_biom", obj.option("otu_taxon_biom"))
            self.option("otu_taxon_table", obj.option("otu_taxon_table"))
            self.option("otu_taxon_dir", obj.option("otu_taxon_dir"))
        # if event['data'] is "est":
        #     self.option("estimators").set_path(obj.option("estimators").prop['path'])
        # if event['data'] is "rarefy":
        #     self.option("rarefaction").set_path(obj.option("rarefaction").prop['path'])

    def run(self):
        self.run_otu()
        self.on_rely([self.stat, self.est, self.rarefy], self.end)
        super(MetaBaseWorkflow, self).run()

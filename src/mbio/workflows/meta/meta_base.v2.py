# -*- coding: utf-8 -*-
# __author__ = 'yuguo'

"""多样性基础分析"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os


class MetaBaseWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        """
        self._sheet = wsheet_object
        super(MetaBaseWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq,sequence.fastq_dir'},  # 输入的fastq文件或fastq文件夹
            {'name': 'otu_fasta', 'type': 'outfile', 'format': 'sequence.fasta'},  # 输出的合并到一起的fasta，供后续的otu分析用
            {'name': 'identity', 'type': 'float', 'default': 0.97},  # 相似性值，范围0-1.
            {'name': 'otu_table', 'type': 'outfile', 'format': 'meta.otu.otu_table'},  # 输出结果otu表
            {'name': 'otu_rep', 'type': 'outfile', 'format': 'sequence.fasta'},  # 输出结果otu代表序列
            # {'name': 'otu_seqids', 'type': 'outfile', 'format': 'meta.otu.otu_seqids'},  # 输出结果otu中包含序列列表
            {'name': 'otu_biom', 'type': 'outfile', 'format': 'meta.otu.biom'},  # 输出结果biom格式otu表
            {'name': 'revcomp', 'type': 'bool'},  # 序列是否翻转
            {'name': 'confidence', 'type': 'float', 'default': 0.7},  # 置信度值
            # {"name": "customer_mode", "type": "bool", "default": False},  # customer 自定义数据库
            {'name': 'database', 'type': 'string'},  # 数据库选择
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 参考fasta序列
            {'name': 'ref_taxon', 'type': 'infile', 'format': 'taxon.seq_taxon'},  # 参考taxon文件
            {'name': 'taxon_file', 'type': 'outfile', 'format': 'taxon.seq_taxon'},  # 输出序列的分类信息文件
            {'name': 'otu_taxon_dir', 'type': 'outfile', 'format': 'meta.otu.tax_summary_dir'},  # 输出的otu_taxon_dir文件夹
            {"name": "estimate_indices", "type": "string", "default": "ace-chao-shannon-simpson-coverage"},
            {"name": "rarefy_indices", "type": "string", "default": "sobs-shannon"},  # 指数类型
            {"name": "rarefy_freq", "type": "int", "default": 100},
            {"name": "alpha_level", "type": "string", "default": "otu"},  # level水平
            {"name": "beta_analysis", "type": "string",
                "default": "pca,hcluster"},
            {"name": "beta_level", "type": "string", "default": "otu"},
            {"name": "dis_method", "type": "string", "default": "bray_curtis"},
            # {"name": "phy_newick", "type": "infile", "format": "meta.beta_diversity.newick_tree"},
            {"name": "permutations", "type": "int", "default": 999},
            {"name": "linkage", "type": "string", "default": "average"},
            {"name": "envtable", "type": "infile", "format": "meta.env_table"}
            # {"name": "group", "type": "infile", "format": "meta.otu.group_table"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.qc = self.add_module("meta.qc.miseq_qc")
        self.otu = self.add_tool("meta.otu.usearch_otu")
        self.phylo = self.add_tool("phylo.phylo_tree")
        self.tax = self.add_tool("taxon.qiime_assign")
        self.stat = self.add_tool("meta.otu.otu_taxon_stat")
        # self.est = self.add_tool("meta.alpha_diversity.estimators")
        # self.rarefy = self.add_tool("meta.alpha_diversity.rarefaction")
        self.alpha = self.add_module("meta.alpha_diversity.alpha_diversity")
        self.beta = self.add_module("meta.beta_diversity.beta_diversity")
        self.step.add_steps("qcstat", "otucluster", "taxassign", "alphadiv", "betadiv")

    def check_options(self):
        """
        检查参数设置
        """
        # if not self.option("fasta").is_set:
        #     raise OptionError("必须设置输入fasta文件.")
        if self.option("identity") < 0 or self.option("identity") > 1:
            raise OptionError("identity值必须在0-1范围内.")
        if self.option("revcomp") not in [True, False]:
            raise OptionError("必须设置参数revcomp")
        if self.option('database') == "custom_mode":
            if not self.option("ref_fasta").is_set or not self.option("ref_taxon").is_set:
                raise OptionError("数据库自定义模式必须设置ref_fasta和ref_taxon参数")
        else:
            if self.option("database") not in ['silva119/16s_bacteria', 'silva119/16s_archaea',
                                               'silva119/18s_eukaryota', 'unite6.0/its_fungi', 'fgr/amoA', 'fgr/nosZ',
                                               'fgr/nirK', 'fgr/nirS', 'fgr/nifH', 'fgr/pmoA', 'fgr/mmoX']:
                raise OptionError("数据库{}不被支持".format(self.option("database")))
        return True

    def run_qc(self):
        # if self.option('multi_fastq').is_set:
        #     self.qc.set_options({
        #         'multi_fastq': self.option('multi_fastq')
        #     })
        # if self.option('single_fastq').is_set:
        #     self.qc.set_options({
        #         'single_fastq': self.option('single_fastq')
        #     })

        self.qc.set_options({
                "in_fastq": self.option("in_fastq")
            })
        # if self.option("fastq").format is 'sequence.fastq_dir':
        #     self.qc.set_options({
        #         "filename_sample": self.option("filename_sample")
        #         })
        # if self.option("fastq").format is 'sequence.fastq':
        #     self.qc.set_options({
        #         "seqname_sample": self.option("seqname_sample")
        #         })
        self.qc.on("start", self.set_step, {'start': self.step.qcstat})
        self.qc.on("end", self.set_output, "qc")
        self.qc.run()

    def run_otu(self):
        self.otu.set_options({
            "fasta": self.option("otu_fasta"),
            "identity": self.option("identity")
            })
        self.otu.on("end", self.set_output, "otu")
        self.otu.on("start", self.set_step, {'start': self.step.otucluster})
        self.otu.run()

    def run_phylotree(self):
        self.phylo.set_options({
            "fasta_file": self.otu.output_dir+"/otu_reps.fasta"
            })
        self.phylo.run()

    def run_taxon(self, relyobj):
        self.tax.set_options({
            "fasta": relyobj.rely[0].option("otu_rep"),
            "revcomp": self.option("revcomp"),
            "confidence": self.option("confidence"),
            "database": self.option("database")}
            )
        if self.option("database") == "custom_mode":
            self.tax.set_options({
                "ref_fasta": self.option("ref_fasta"),
                "ref_taxon": self.option("ref_taxon")
            })
        self.tax.on("end", self.set_output, "tax")
        self.tax.on("start", self.set_step, {'end': self.step.otucluster, 'start': self.step.taxassign})
        self.tax.run()

    def run_stat(self):
        self.stat.set_options({
            "in_otu_table": self.option("otu_table"),
            "taxon_file": self.option("taxon_file")
            })
        self.stat.on("end", self.set_output, "stat")
        self.stat.on("end", self.set_step, {'end': self.step.taxassign})
        self.stat.run()

    def run_alpha(self):
        self.alpha.set_options({
            'otu_table': self.option('otu_taxon_dir'),
            "level": self.option('alpha_level'),
            'estimate_indices': self.option('estimate_indices'),
            'rarefy_indices': self.option('rarefy_indices'),
            'rarefy_freq': self.option('rarefy_freq')
            })
        self.alpha.on("end", self.set_output, "alpha")
        self.alpha.on("start", self.set_step, {'start': self.step.alphadiv})
        self.alpha.on("end", self.set_step, {'end': self.step.alphadiv})
        self.alpha.run()

    def run_beta(self):
        self.beta.set_options({
            'analysis': self.option('beta_analysis'),
            'dis_method': self.option('dis_method'),
            'otutable': self.option('otu_taxon_dir'),
            "level": self.option('beta_level'),
            'permutations': self.option('permutations'),
            'phy_newick': self.phylo.option('phylo_tre').prop['path']
            })
        # if self.option('phy_newick').is_set:
        #     self.beta.set_options({
        #         'phy_newick': self.option('phy_newick')
        #         })
        # if self.option('envtable').is_set:
        #     self.beta.set_options({
        #         'envtable': self.option('envtable')
        #         })
        # if self.option('group').is_set:
        #     self.beta.set_options({
        #         'group': self.option('group')
        #         })
        self.beta.on("end", self.set_output, "beta")
        self.beta.on("start", self.set_step, {'start': self.step.betadiv})
        self.beta.on("end", self.set_step, {'end': self.step.betadiv})
        self.beta.run()

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def set_output(self, event):
        obj = event["bind_object"]
        if event['data'] is "qc":
            self.option("otu_fasta", obj.option("otu_fasta"))
            os.system('cp -r '+obj.output_dir+' '+self.output_dir+"/QC_stat")
        if event['data'] is "otu":
            self.option("otu_table", obj.option("otu_table"))
            self.option("otu_rep", obj.option("otu_rep"))
            self.option("otu_biom", obj.option("otu_biom"))
            os.system('cp -r '+obj.output_dir+' '+self.output_dir+"/Otu")
        if event['data'] is "tax":
            self.option("taxon_file", obj.option("taxon_file"))
            os.system('cp -r '+obj.output_dir+' '+self.output_dir+"/Tax_assign")
        if event['data'] is "stat":
            # self.option("otu_taxon_biom", obj.option("otu_taxon_biom"))
            # self.option("otu_taxon_table", obj.option("otu_taxon_table"))
            self.option("otu_taxon_dir", obj.option("otu_taxon_dir"))
            os.system('cp -r '+obj.output_dir+' '+self.output_dir+"/OtuTaxon_summary")
        if event['data'] is "alpha":
            os.system('cp -r '+obj.output_dir+' '+self.output_dir+"/Alpha_diversity")
        if event['data'] is "beta":
            os.system('cp -r '+obj.output_dir+' '+self.output_dir+"/Beta_diversity")

    def set_db(self):
        # 设置报告文件到数据库
        # 设置QC报告文件
        api_samples = self.api.sample
        sample_info_path = self.qc.output_dir+"/samples_info/samples_info.txt"
        if not os.path.isfile(sample_info_path):
            raise Exception("找不到报告文件:{}".format(sample_info_path))
        api_samples.add_samples_info(sample_info_path)
        base_info_path = ""
        with open(self.qc.output_dir+"/samples_info/samples_info.txt") as f:
            f.readline()
            for line in f:
                s = line.split('\t')[0]
                base_info_path = self.qc.output_dir+"/base_info/{}.fastq.fastxstat.txt".format(s)
                if not os.path.isfile(base_info_path):
                    raise Exception("找不到报告文件:{}".format(base_info_path))
                api_samples.add_base_info(s, base_info_path)
        for step in (20, 50, 100, 200):
            reads_len_info_path = self.qc.output_dir+"/reads_len_info/step_{}.reads_len_info.txt".format(str(step))
            if not os.path.isfile(reads_len_info_path):
                raise Exception("找不到报告文件:{}".format(base_info_path))
            api_samples.add_reads_len_info(step, reads_len_info_path)
        # 设置OTU table文件
        api_otu = self.api.meta
        otu_path = self.output_dir+"/OtuTaxon_summary/otu_taxon.xls"
        rep_path = self.output_dir+"/Otu/otu_reps.fasta"
        if not os.path.isfile(otu_path):
            raise Exception("找不到报告文件:{}".format(otu_path))
        otu_id = api_otu.add_otu_table(otu_path, major=True, rep_path=rep_path)
        # 设置进化树文件
        api_tree = self.api.newicktree
        tree_path = self.phylo.option('phylo_tre').prop['path']
        if not os.path.isfile(tree_path):
            raise Exception("找不到报告文件:{}".format(tree_path))
        api_tree.add_tree_file(tree_path, major=True, level=9, table_id=otu_id, table_type='otu', tree_type='phylo')
        # 设置alpha多样性文件
        api_est = self.api.estimator
        est_path = self.output_dir+"/Alpha_diversity/estimators.xls"
        if not os.path.isfile(est_path):
            raise Exception("找不到报告文件:{}".format(est_path))
        api_est.add_est_table(est_path, major=True, level=9, otu_id=otu_id)
        api_rare = self.api.rarefaction
        rare_path = self.work_dir+"/AlphaDiversity/Rarefaction/output/"
        api_rare.add_rare_table(rare_path, level=9, otu_id=otu_id)
        # 设置beta多样性文件
        api_dist = self.api.distance
        dist_path = self.beta.option('dis_matrix').prop['path']
        if not os.path.isfile(dist_path):
            raise Exception("找不到报告文件:{}".format(dist_path))
        dist_id = api_dist.add_dist_table(dist_path, level=9, otu_id=otu_id, major=True)
        # 设置hcluster树文件
        api_hcluster = self.api.newicktree
        hcluster_path = self.beta.output_dir+"/Hcluster/hcluster.tre"
        if not os.path.isfile(hcluster_path):
            raise Exception("找不到报告文件:{}".format(hcluster_path))
        api_hcluster.add_tree_file(hcluster_path, major=True, table_id=dist_id, table_type='dist', tree_type='cluster')
        self.end()

    def run(self):
        self.run_qc()
        self.on_rely(self.qc, self.run_otu)
        self.on_rely(self.otu, self.run_taxon)
        self.on_rely(self.otu, self.run_phylotree)
        self.on_rely(self.tax, self.run_stat)
        self.on_rely(self.stat, self.run_alpha)
        self.on_rely([self.stat, self.phylo], self.run_beta)
        self.on_rely([self.alpha, self.beta], self.set_db)
        super(MetaBaseWorkflow, self).run()

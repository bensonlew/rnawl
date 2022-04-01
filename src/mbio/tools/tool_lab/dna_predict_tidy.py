# -*- coding: utf-8 -*-
# 20210429
import os
import re
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import Counter


class DnaPredictTidyAgent(Agent):
    """
    prodigal小工具
    """

    def __init__(self, parent):
        super(DnaPredictTidyAgent, self).__init__(parent)
        options = [
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  # 基因预测结果
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "gene_prefix", "type": "string", "default": "ORF"},  # 基因前缀，自定义前缀,默认gene
            {"name": "gene", "type": "outfile", "format": "gene_structure.gff3"},  # 基因预测结果,排序整理后的，并加上累积长度起始位置
            {"name": "gene_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "gene_faa", "type": "outfile", "format": "sequence.fasta"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("gene_gff").is_set:
            raise OptionError("必须设置参数gene_gff", code="31403101")
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome", code="31403104")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(DnaPredictTidyAgent, self).end()


class DnaPredictTidyTool(Tool):
    def __init__(self, config):
        super(DnaPredictTidyTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.sample_name = (self.genome_fasta).split('/')[-1].rsplit('.')[0]
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/bioinfo/seq/EMBOSS-6.6.0/lib"
        self.set_environ(LD_LIBRARY_PATH=self.LD_LIBRARY_PATH)
        self.transeq = "/bioinfo/seq/EMBOSS-6.6.0/emboss/transeq"

    def run_drop_pseudogene_by_trna_and_sort(self):
        gene_table_path = self.option("gene_gff").prop["path"]
        gene_table = pd.read_table(gene_table_path, sep='\t', header=0)
        gene = gene_table.ix[:, ["Sequence id", "Start", "End"]]
        if len(gene_table) > 0:
            gene["location"] = gene["Sequence id"].str.split(" ", expand=True)[0]
            gene_sort = gene["location"].str.split("caffold", expand=True)
        summary = pd.concat([gene.ix[:, ["Sequence id", "s", "sort"]]])
        summary = summary.sort_values(by=["sort", "s"])
        gene_id_list = []
        gene_suffix_len = len(str(len(summary)))
        if gene_suffix_len < 5:
            gene_suffix_len = 4
        for n in range(len(summary)):
            k = n + 1
            gene_id = self.option("gene_prefix") + "0" * (gene_suffix_len - len(str(k))) + str(k)
            gene_id_list.append(gene_id)
        summary["Gene ID"] = gene_id_list
        del summary["s"]
        del summary["sort"]
        # "提取排序好的gene_id, 形成各预测的gff文件"
        gene = summary.merge(gene_table, on='Sequence id', how='inner')
        del gene["Gene id"]
        gene_file_path = self.work_dir + "/predict.gff"
        gene = gene.set_index("Gene ID")
        if len(gene) > 0:
            location = gene["Sequence id"].str.split("_", expand=True)[0]
            loc_nu = (Counter(location))
            for i in range(len(gene)):
                seq_id = gene["Sequence id"][i]
                gene_suffix_len = len(str(loc_nu[seq_id.split("_ORF")[0]]))
                nu = seq_id.split("_ORF")[1].lstrip("0")
                tmp = "_ORF" + "0" * (gene_suffix_len - len(nu)) + nu
                gene["Sequence id"][i] = re.subn("_ORF[0-9]*", tmp, gene["Sequence id"][i])[0]
        gene.to_csv(gene_file_path, sep="\t")

    def run_get_fa_bygff(self, genome, gff):
        cmd = '{} {}get_fa_bygff.pl {} {} {} {}'.format(self.perl_path, self.config.PACKAGE_DIR + "/tool_lab/", genome,
                                                        gff, "0", self.output_dir)
        name = os.path.basename(gff).split(".")[0]
        command = self.add_command(str(name).lower(), cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取fnn和整理GFF文件运行完成")
        else:
            self.set_error("提取fnn和整理GFF文件运行出错!", code="31403101")

    def run_transeq(self):
        nul_seq = self.output_dir + "/" + self.sample_name + "_predict.fnn"
        port_seq = self.work_dir + "/" + self.sample_name + "_predict.faa"
        cmd = '{} -trim -table 11 -sequence {} -outseq {}'.format(self.transeq, nul_seq, port_seq)
        command = self.add_command("transeq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("翻译蛋白序列运行完成")
            self.run_faa2m()
        else:
            self.set_error("翻译蛋白序列运行出错!", code="31403102")

    def run_faa2m(self):
        cmd = '{} {}faa2m.pl {} {}'.format(self.perl_path, self.perl_script,
                                           self.work_dir + "/" + self.sample_name + "_predict.faa",
                                           self.output_dir + "/" + self.sample_name + "_predict.faa")
        command = self.add_command("faa2m", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("矫正起始密码子运行完成")
        else:
            self.set_error("矫正起始密码子运行出错!", code="31403103")

    def set_output(self):
        self.option('gene', self.output_dir + "/" + self.sample_name + "_predict.gff")
        if os.path.getsize(self.output_dir + "/" + self.sample_name + "_predict.fnn"):
            self.option('gene_fnn', self.output_dir + "/" + self.sample_name + "_predict.fnn")
            self.option('gene_faa', self.output_dir + "/" + self.sample_name + "_predict.faa")

    def run(self):
        super(DnaPredictTidyTool, self).run()
        self.run_drop_pseudogene_by_trna_and_sort()
        self.run_get_fa_bygff(self.option("input_genome").prop["path"], self.work_dir + "/predict.gff")
        if os.path.getsize(self.output_dir + "/" + self.sample_name + "_predict.fnn"):
            self.run_transeq()
        self.set_output()
        self.end()

# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# version 1.0
# last_modify: 2019.09.20

import os
import re
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import Counter


class DnaPredictTidyStatAgent(Agent):
    """
    对基因预测、tRNA预测和rRNA预测的结果进行整理：
    1. 利用tRNA结果去除假基因（ORF预测基因中与tRNA预测结果重叠的部分）；
    2. 将预测结果按位置信息进行统一排序；
    3. 统一基因id：默认为gene0001开始，可设置前缀替换“gene”
    """

    def __init__(self, parent):
        super(DnaPredictTidyStatAgent, self).__init__(parent)
        options = [
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  # 基因预测结果
            {"name": "trna_gff", "type": "infile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "rrna_gff", "type": "infile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "gene_prefix", "type": "string", "default": "gene"},  # 基因前缀，自定义前缀,默认gene
            {"name": "plasmid_prefix", "type": "string", "default": ""},
            # 质粒前缀,字典形式含一个或多个质粒，"{‘plassmidA’: ‘plasa', ‘plassmidB’: ‘plasb'}"
            {"name": "gene", "type": "outfile", "format": "gene_structure.gff3"},  # 基因预测结果,排序整理后的，并加上累积长度起始位置
            {"name": "trna", "type": "outfile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "rrna", "type": "outfile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "gene_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "gene_faa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "trna_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "rrna_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "genome_name", "type": "string"}, #用于设置文件名称
            {"name": "genome_type", "type": "string"}, #基因组的类型，是完成图还是扫描图，这个字段主要是用于计算N50进行判断
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("gene_gff").is_set:
            raise OptionError("必须设置参数gene_gff")
        # if not self.option("trna_gff").is_set:
        #     raise OptionError("必须设置参数trna_gff")
        # if not self.option("rrna_gff").is_set:
        #     raise OptionError("必须设置参数rrna_gff")
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(DnaPredictTidyStatAgent, self).end()


class DnaPredictTidyStatTool(Tool):
    def __init__(self, config):
        super(DnaPredictTidyStatTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.sample_name = self.option('genome_name')
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.transeq = "/bioinfo/seq/EMBOSS-6.6.0/emboss/transeq"
        self.scripts = self.config.PACKAGE_DIR + "/bac_comp_genome/"
        self.LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/bioinfo/seq/EMBOSS-6.6.0/lib:" + self.config.SOFTWARE_DIR + "/library/lib:" + self.config.SOFTWARE_DIR + "/library/lib64"
        self.set_environ(LD_LIBRARY_PATH=self.LD_LIBRARY_PATH)

    def run_drop_pseudogene_by_trna_and_sort(self):
        """
        根据tRNA结果校正基因序列文件和排序
        :return:
        """
        gene_table_path = self.option("gene_gff").prop["path"]
        gene_table = pd.read_table(gene_table_path, sep='\t', header=0)
        gene = gene_table.ix[:, ["Sequence id", "Start", "End"]]
        if len(gene_table) > 0:
            gene["location"] = gene["Sequence id"].str.split("_", expand=True)[0]
            gene_sort = gene["location"].str.split("caffold", expand=True)
        trna_table_path = self.option("trna_gff").prop["path"]
        trna_table = pd.read_table(trna_table_path, sep='\t', header=0)
        trna = trna_table.ix[:, ["Sequence id", "Start", "End"]]
        trna["s"] = trna.ix[:, ["Start", "End"]].apply(lambda x: x.min(), axis=1)
        trna["e"] = trna.ix[:, ["Start", "End"]].apply(lambda x: x.max(), axis=1)
        if len(trna_table) > 0:
            trna["location"] = trna["Sequence id"].str.split("_", expand=True)[0]
            trna_sort = trna["location"].str.split("caffold", expand=True)
            if len(gene_table) > 0:
                gene_drop = []
                for n in range(len(trna)):
                    s = trna.iat[n, 3]
                    e = trna.iat[n, 4]
                    loca = trna.iat[n, 5]
                    tmp = gene[gene["location"] == loca]
                    if len(tmp) > 0:
                        drop = tmp[((s < tmp["Start"]) & (tmp["Start"] < e)) | ((s < tmp["End"]) & (tmp["End"] < e))]
                        if list(drop["Sequence id"]) > 0:
                            gene_drop = gene_drop + list(drop["Sequence id"])
                gene.index = gene["Sequence id"]
                if len(gene_drop) > 0:
                    gene.index = gene["Sequence id"]
                    gene = gene.drop(gene_drop, axis=0)
                gene["s"] = gene.ix[:, ["Start", "End"]].apply(lambda x: x.min(), axis=1)
                gene_sort = gene["location"].str.split("caffold", expand=True)
        rrna_table_path = self.option("rrna_gff").prop["path"]
        rrna_table = pd.read_table(rrna_table_path, sep='\t', header=0)
        rrna = rrna_table.ix[:, ["Sequence id", "Start", "End"]]
        rrna["s"] = rrna.ix[:, ["Start", "End"]].apply(lambda x: x.min(), axis=1)
        if len(rrna_table) > 0:
            rrna["location"] = rrna["Sequence id"].str.split("_", expand=True)[0]
            rrna_sort = rrna["location"].str.split("caffold", expand=True)
        if len(list(gene["location"])) > 1:
            if re.match("^[A-Za-z]caffold[0-9]+$", list(gene["location"])[1]):
                if len(trna_table) > 0:
                    trna["sort"] = trna_sort[1].astype(int)
                if len(gene_table) > 0:
                    gene["sort"] = gene_sort[1].astype(int)
                if len(rrna_table) > 0:
                    rrna["sort"] = rrna_sort[1].astype(int)
        else:
            if len(trna_table) > 0:
                trna["sort"] = trna_sort[0]
            if len(gene_table) > 0:
                gene["sort"] = gene_sort[0]
            if len(rrna_table) > 0:
                rrna["sort"] = rrna_sort[0]
        summary = pd.concat([gene.ix[:, ["Sequence id", "s", "sort"]], trna.ix[:, ["Sequence id", "s", "sort"]],
                             rrna.ix[:, ["Sequence id", "s", "sort"]]])
        summary = summary.sort_values(by=["sort", "s"])
        gene_id_list = []
        if self.option("plasmid_prefix") != "":
            pre = eval(self.option("plasmid_prefix"))  # 如果字符串引号为双引号，则通不过该不
            for k in sorted(pre.keys()):
                if len(summary) >0:
                    plas_one = summary[summary["sort"] == k]
                    gene_suffix_len = len(str(len(plas_one)))
                    if gene_suffix_len < 5:
                        gene_suffix_len = 4
                    for n in range(len(plas_one)):
                        p = n + 1
                        gene_id = pre[k] + "0" * (gene_suffix_len - len(str(p))) + str(p)
                        gene_id_list.append(gene_id)
        else:
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
        rrna = summary.merge(rrna_table, on='Sequence id', how='inner')
        del rrna["Sequence Name"]
        del rrna["Type"]
        rrna_file_path = self.work_dir + "/rRNA.gff"
        rrna = rrna.set_index("Gene ID")
        if len(rrna) > 0:
            location = rrna["Sequence id"].str.split("_", expand=True)[0]
            loc_nu = (Counter(location))
            for i in range(len(rrna)):
                seq_id = rrna["Sequence id"][i]
                gene_suffix_len = len(str(loc_nu[seq_id.split("_rRNA")[0]]))
                nu = seq_id.split("_rRNA")[1].lstrip("0")
                tmp = "_rRNA" + "0" * (gene_suffix_len - len(nu)) + nu
                rrna["Sequence id"][i] = re.subn("_rRNA[0-9]*", tmp, rrna["Sequence id"][i])[0]
        rrna.to_csv(rrna_file_path, sep="\t")
        trna = summary.merge(trna_table, on='Sequence id', how='inner')
        del trna["Sequence Name"]
        trna_file_path = self.work_dir + "/tRNA.gff"
        trna = trna.set_index("Gene ID")
        if len(trna) > 0:
            location = trna["Sequence id"].str.split("_", expand=True)[0]
            loc_nu = (Counter(location))
            for i in range(len(trna)):
                seq_id = trna["Sequence id"][i]
                gene_suffix_len = len(str(loc_nu[seq_id.split("_tRNA")[0]]))
                nu = seq_id.split("_tRNA")[1].lstrip("0")
                tmp = "_tRNA" + "0" * (gene_suffix_len - len(nu)) + nu
                trna["Sequence id"][i] = re.subn("_tRNA[0-9]*", tmp, trna["Sequence id"][i])[0]
        trna.to_csv(trna_file_path, sep="\t")

    def run_get_fa_bygff(self, genome, gff, type):
        """
        根据gff文件从中提取出序列和整理标准化gff文件
        :param genome:
        :param gff:
        :return:
        """
        cmd = '{} {}get_fa_bygff.pl {} {} {} {}'.format(self.perl_path, self.perl_script, genome,
                                                        gff, "0", self.output_dir)
        name = os.path.basename(gff).split(".gff")[0]
        command = self.add_command(str(name).lower(), cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取fnn和整理GFF文件运行完成")
        else:
            self.set_error("提取fnn和整理GFF文件运行出错!")

    def run_transeq(self):
        """
        对提取的序列进行翻译
        :return:
        """
        nul_seq = self.output_dir + "/" + self.sample_name + "_predict.ffn"
        port_seq = self.work_dir + "/" + self.sample_name + "_predict.faa"
        cmd = '{} -trim -table 11 -sequence {} -outseq {}'.format(self.transeq, nul_seq, port_seq)
        command = self.add_command("transeq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("翻译蛋白序列运行完成")
            self.run_faa2m()
        else:
            self.set_error("翻译蛋白序列运行出错!")

    def run_faa2m(self):
        """
        校正起始密码子
        :return:
        """
        cmd = '{} {}faa2m.pl {} {}'.format(self.perl_path, self.perl_script,
                                           self.work_dir + "/" + self.sample_name + "_predict.faa",
                                           self.output_dir + "/" + self.sample_name + "_predict.faa")
        command = self.add_command("faa2m", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("矫正起始密码子运行完成")
        else:
            self.set_error("矫正起始密码子运行出错!")

    def run_base_stat(self):
        """
        统计基本信息 共有八个参数
        :return:
        """
        outfile_path = self.output_dir + '/' +  self.sample_name + "_stat.xls"
        self.logger.info("正在统计基础信息")
        cmd = '{} {}base_n50.pl {} {} {} {} {} {} {} {}'.format(
            self.perl_path, self.scripts, self.option('input_genome').prop['path'], self.sample_name,
            self.option('genome_name'), self.option('genome_type'), self.option('gene_gff').prop['path'],
            self.option('rrna_gff').prop['path'], self.option('trna_gff').prop['path'], outfile_path
        )
        self.logger.info(cmd)
        command = self.add_command("stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("统计基础信息运行完成")
        else:
            self.set_error("统计基础信息运行出错!")

    def set_output(self):
        self.logger.info("开始设置结果目录")
        if os.path.exists(self.output_dir + "/" + self.sample_name + "_CDS.gff"):
            self.option('gene', self.output_dir + "/" + self.sample_name + "_CDS.gff")
        elif os.path.exists(self.output_dir + "/" + self.sample_name + "_tRNA.gff"):
            self.option('trna', self.output_dir + "/" + self.sample_name + "_tRNA.gff")
        elif os.path.exists(self.output_dir + "/" + self.sample_name + "_rRNA.gff"):
            self.option('rrna', self.output_dir + "/" + self.sample_name + "_rRNA.gff")

        if os.path.exists(self.output_dir + "/" + self.sample_name + "_CDS.ffn") and os.path.getsize(self.output_dir + "/" + self.sample_name + "_CDS.ffn"):
            self.option('gene_fnn', self.output_dir + "/" + self.sample_name + "_CDS.ffn")
            self.option('gene_faa', self.output_dir + "/" + self.sample_name + "_CDS.faa")
        if os.path.exists(self.output_dir + "/" + self.sample_name + "_tRNA.ffn") and os.path.getsize(self.output_dir + "/" + self.sample_name + "_tRNA.ffn"):
            self.option('trna_fnn', self.output_dir + "/" + self.sample_name + "_tRNA.ffn")
        if os.path.exists(self.output_dir + "/" + self.sample_name + "_rRNA.ffn") and os.path.getsize(self.output_dir + "/" + self.sample_name + "_rRNA.ffn"):
            self.option('rrna_fnn', self.output_dir + "/" + self.sample_name + "_rRNA.ffn")
        self.logger.info("设置结果目录完成")

    def run(self):
        super(DnaPredictTidyStatTool, self).run()
        if self.option("trna_gff").is_set and self.option("rrna_gff").is_set:
            self.run_drop_pseudogene_by_trna_and_sort()
            for i in ["CDS.gff", "tRNA.gff", "rRNA.gff"]:
                if re.search(r'CDS', i):
                    type = 'CDS'
                elif re.search(r'tRNA', i):
                    type = 'tRNA'
                else:
                    type = 'rRNA'
                self.run_get_fa_bygff(self.option("input_genome").prop["path"], self.work_dir +'/'+ i, type)
        else:
            if self.option("trna_gff").is_set and (not self.option("rrna_gff").is_set):
                gene_gff = self.option("gene_gff").prop['path']
                trna_gff = self.option("trna_gff").prop["path"]
                self.run_get_fa_bygff(self.option("input_genome").prop["path"], gene_gff, 'CDS')
                self.run_get_fa_bygff(self.option("input_genome").prop["path"], trna_gff, 'tRNA')
            elif self.option("rrna_gff").is_set and (not self.option("trna_gff").is_set):
                gene_gff = self.option("gene_gff").prop['path']
                rrna_gff = self.option("rrna_gff").prop["path"]
                self.run_get_fa_bygff(self.option("input_genome").prop["path"], gene_gff, 'CDS')
                self.run_get_fa_bygff(self.option("input_genome").prop["path"], rrna_gff, 'rRNA')
            else:
                gene_gff = self.option("gene_gff").prop['path']
                self.run_get_fa_bygff(self.option("input_genome").prop["path"], gene_gff, 'CDS')
        self.run_drop_pseudogene_by_trna_and_sort()
        for i in ["/predict.gff", "/tRNA.gff", "/rRNA.gff"]:
            self.run_get_fa_bygff(self.option("input_genome").prop["path"], self.work_dir + i)
        if os.path.getsize(self.output_dir + "/" + self.sample_name + "_predict.fna"):
            self.run_transeq()
        self.set_output()
        self.end()

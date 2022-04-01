# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# version 1.0
# last_modify: 20190915

import os
import re
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import Counter


class DnaPredictTidyAgent(Agent):
    """
    对基因预测、tRNA预测和rRNA预测的结果进行整理：
    1. 利用tRNA结果去除假基因（ORF预测基因中与tRNA预测结果重叠的部分）；
    2. 将预测结果按位置信息进行统一排序；
    3. 统计序列的基本信息结果这里主要是统计GC含量、genome的大小、N50、cds数目、rRNA数目和tRNA数目;
    """

    def __init__(self, parent):
        super(DnaPredictTidyAgent, self).__init__(parent)
        options = [
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  # 基因预测结果
            {"name": "trna_gff", "type": "infile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "rrna_gff", "type": "infile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "gene_prefix", "type": "string"},  # 基因前缀，自定义前缀,默认选择物种名称的前三位+随机三位数字
            {"name": "plasmid_prefix", "type": "string",'default': ""},
            # 质粒前缀,字典形式含一个或多个质粒，"{‘plassmidA’: ‘plasa', ‘plassmidB’: ‘plasb'}"
            {"name": "genome_type", "type": "string"}, #基因组的类型，是完成图还是扫描图，这个字段主要是用于计算N50进行判断
            {"name": "genome_name", "type": "string"}, #用于设置文件名称
            {"name": "sample", "type": "string"}, #用于设置文件名称
            {"name": "gene", "type": "outfile", "format": "gene_structure.gff3"},  # 基因预测结果,排序整理后的，并加上累积长度起始位置
            {"name": "trna", "type": "outfile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "rrna", "type": "outfile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "gene_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "gene_faa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "trna_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "rrna_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "is_gff", "type": "string", "default": "seq"},##类型为seq或者gff+seq
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
        if not self.option("genome_type"):
            raise OptionError("必须设置参数基因组类型")
        if not self.option("genome_name"):
            raise OptionError("必须设置参数基因组名称")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(DnaPredictTidyAgent, self).end()


class DnaPredictTidyTool(Tool):
    def __init__(self, config):
        super(DnaPredictTidyTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        if self.option("genome_name"):
            self.sample_name = self.option("genome_name")
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.transeq = "/bioinfo/seq/EMBOSS-6.6.0/emboss/transeq"
        self.scripts = self.config.PACKAGE_DIR + "/bac_comp_genome/"
        self.LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/bioinfo/seq/EMBOSS-6.6.0/lib:" + self.config.SOFTWARE_DIR + "/library/lib:" + self.config.SOFTWARE_DIR + "/library/lib64"
        self.set_environ(LD_LIBRARY_PATH=self.LD_LIBRARY_PATH)

    def run_drop_pseudogene_by_trna_and_sort(self):
        """
        根据tRNA预测结果去除掉假基因结果
        :return:
        """
        if self.option("gene_gff").is_set:
            gene_gff = self.option('gene_gff').prop['path']
        else:
            gene_gff = self.work_dir + "/CDS.gff"
            os.system('touch %s'%(gene_gff))
        if self.option('rrna_gff').is_set:
            rrna_gff = self.option('rrna_gff').prop['path']
        else:
            rrna_gff = self.work_dir + "/rRNA.gff"
            os.system('touch %s'%(rrna_gff))
        if self.option('trna_gff').is_set:
            trna_gff = self.option('trna_gff').prop['path']
        else:
            trna_gff = self.work_dir + "/tRNA.gff"
            os.system('touch %s'%(trna_gff))

        gene_table_path = gene_gff
        gene_table = pd.read_table(gene_table_path, sep='\t', header=0)
        gene = gene_table.ix[:, ["Sequence id", "Start", "End"]]
        if len(gene_table) > 0:
            gene["location"] = gene["Sequence id"].str.split("\_ORF", expand=True)[0]
            gene_sort = gene["location"].str.split("caffold", expand=True)
            self.logger.info("+++++++++gene_sort{}\n".format(gene_sort.head()))
        trna_table_path = trna_gff
        trna_table = pd.read_table(trna_table_path, sep='\t', header=0)
        trna = trna_table.ix[:, ["Sequence id", "Start", "End"]]
        trna["s"] = trna.ix[:, ["Start", "End"]].apply(lambda x: x.min(), axis=1)
        trna["e"] = trna.ix[:, ["Start", "End"]].apply(lambda x: x.max(), axis=1)
        if len(trna_table) > 0:
            trna["location"] = trna["Sequence id"].str.split("\_tRNA", expand=True)[0]
            trna_sort = trna["location"].str.split("caffold", expand=True)
            self.logger.info("+++++++++trna_sort{}\n".format(trna_sort.head()))
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
        rrna_table_path = rrna_gff
        rrna_table = pd.read_table(rrna_table_path, sep='\t', header=0)
        rrna = rrna_table.ix[:, ["Sequence id", "Start", "End"]]
        rrna["s"] = rrna.ix[:, ["Start", "End"]].apply(lambda x: x.min(), axis=1)
        if len(rrna_table) > 0:
            rrna["location"] = rrna["Sequence id"].str.split("\_rRNA", expand=True)[0]
            rrna_sort = rrna["location"].str.split("caffold", expand=True)
            self.logger.info("+++++++++rrna_sort{}\n".format(rrna_sort.head()))
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
            pre = eval(self.option("plasmid_prefix"))  # 如果字符串引号为双引号，则通不过这一步
            self.logger.info("+++++++++summary{}\n".format(summary.head()))
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
        gene_file_path = self.work_dir + "/CDS.gff"
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
        从gff文件中获取fasta序列
        :param genome:
        :param gff:
        :return:
        """
        self.logger.info("正在从gff文件中获取fasta序列")
        cmd = '{} {}get_fa_bygff2.pl {} {} {} {} {} {} {}'.format(self.perl_path, self.scripts, genome,
                                                        gff, "0", self.output_dir, self.option('genome_name'), type, self.option('genome_type'))
        name = os.path.basename(gff).split(".gff")[0]
        self.logger.info(cmd)
        command = self.add_command(str(name).lower(), cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取fnn和整理GFF文件运行完成")
        else:
            self.set_error("提取fnn和整理GFF文件运行出错!")

    def run_transeq(self):
        """
        根据获得的cds的fasta序列，完成翻译序列的任务
        :return:
        """
        self.logger.info("正在根据获得的fasta序列翻译为蛋白序列")
        nul_seq = self.output_dir + "/" + self.sample_name + "_CDS.fna"
        port_seq = self.work_dir + "/" + self.sample_name + "_CDS.faa"
        cmd = '{} -trim -table 11 -sequence {} -outseq {}'.format(self.transeq, nul_seq, port_seq)
        self.logger.info(cmd)
        command = self.add_command("transeq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("翻译蛋白序列运行完成")
            self.run_faa2m()
        else:
            self.set_error("翻译蛋白序列运行出错!")

    def run_faa2m(self):
        """
        矫正起始密码子
        :return:
        """
        self.logger.info("正在矫正起始密码子")
        cmd = '{} {}faa2m.pl {} {}'.format(self.perl_path, self.perl_script,
                                           self.work_dir + "/" + self.sample_name + "_CDS.faa",
                                           self.output_dir + "/" + self.sample_name + "_CDS.faa")
        self.logger.info(cmd)
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
        if self.option("gene_gff").is_set:
            gene_gff = self.option('gene_gff').prop['path']
        else:
            gene_gff = self.work_dir + "/CDS.gff"
            os.system('touch %s'%(gene_gff))
        if self.option('rrna_gff').is_set:
            rrna_gff = self.option('rrna_gff').prop['path']
        else:
            rrna_gff = self.work_dir + "/rRNA.gff"
            os.system('touch %s'%(rrna_gff))
        if self.option('trna_gff').is_set:
            trna_gff = self.option('trna_gff').prop['path']
        else:
            trna_gff = self.work_dir + "/tRNA.gff"
            os.system('touch %s'%(trna_gff))
        cmd = '{} {}base_n50.pl {} {} {} {} {} {} {} {}'.format(
            self.perl_path, self.scripts, self.option('input_genome').prop['path'], self.sample_name,
            self.option('genome_name'), self.option('genome_type'), gene_gff,
            rrna_gff, trna_gff, outfile_path
        )
        self.logger.info(cmd)
        command = self.add_command("stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("统计基础信息运行完成")
        else:
            self.set_error("统计基础信息运行出错!")


    def set_output(self):
        """
        设置结果文件
        :return:
        """
        self.logger.info("开始设置结果目录")
        if os.path.exists(self.output_dir + "/" + self.sample_name + "_CDS.gff"):
            self.option('gene', self.output_dir + "/" + self.sample_name + "_CDS.gff")
        elif os.path.exists(self.output_dir + "/" + self.sample_name + "_tRNA.gff"):
            self.option('trna', self.output_dir + "/" + self.sample_name + "_tRNA.gff")
        elif os.path.exists(self.output_dir + "/" + self.sample_name + "_rRNA.gff"):
            self.option('rrna', self.output_dir + "/" + self.sample_name + "_rRNA.gff")

        if os.path.exists(self.output_dir + "/" + self.sample_name + "_CDS.fna") and os.path.getsize(self.output_dir + "/" + self.sample_name + "_CDS.fna"):
            self.option('gene_fnn', self.output_dir + "/" + self.sample_name + "_CDS.fna")
            self.option('gene_faa', self.output_dir + "/" + self.sample_name + "_CDS.faa")
        if os.path.exists(self.output_dir + "/" + self.sample_name + "_tRNA.fna") and os.path.getsize(self.output_dir + "/" + self.sample_name + "_tRNA.fna"):
            self.option('trna_fnn', self.output_dir + "/" + self.sample_name + "_tRNA.fna")
        if os.path.exists(self.output_dir + "/" + self.sample_name + "_rRNA.fna") and os.path.getsize(self.output_dir + "/" + self.sample_name + "_rRNA.fna"):
            self.option('rrna_fnn', self.output_dir + "/" + self.sample_name + "_rRNA.fna")
        self.logger.info("设置结果目录完成")

    def run(self):
        """
        运行，为了兼容tRNA和rrna的预测结果为空
        :return:
        """
        self.logger.info("开始运行啦")
        super(DnaPredictTidyTool, self).run()
        if self.option("trna_gff").is_set and self.option("rrna_gff").is_set:
            if self.option("is_gff") not in ["gff+seq"]:## gff+seq的情况不做去重和统一排序的功能，默认是排好序和唯一的
                self.run_drop_pseudogene_by_trna_and_sort()
            else:
                os.system("cp {} {}".format(self.option("trna_gff").prop['path'], self.work_dir + "/tRNA.gff"))
                os.system("cp {} {}".format(self.option("gene_gff").prop['path'], self.work_dir + "/CDS.gff"))
                os.system("cp {} {}".format(self.option("rrna_gff").prop['path'], self.work_dir + "/rRNA.gff"))
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
        if os.path.getsize(self.output_dir + "/" + self.sample_name + "_CDS.fna"):
            self.run_transeq()
        self.run_base_stat()
        self.set_output()
        self.end()

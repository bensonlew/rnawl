# -*- coding: utf-8 -*-
# __author__ = 'chenyanyan, qindanhua'

import os
# import glob
# import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
# from mbio.files.sequence.file_sample import FileSampleFile
import json
# from mbio.packages.gene_structure.snp_anno import snp_anno
import gevent.subprocess as subprocess
import sys
from Bio import SeqIO
import glob2
import gevent
import pandas as pd
from biocluster.config import Config
import shutil

class SamRnaModule(Module):
    """
    star:序列比对
    picard:处理比对结果sam文件
    gatk：snp calling软件
    version 1.0
    author: chenyanyan
    last_modify: 2017.01.03 by qindanhua
    """
    def __init__(self, work_id):
        super(SamRnaModule, self).__init__(work_id)
        options = [
            {"name": "ref_genome", "type": "string", 'default': 'customer_mode'},  # 参考基因组类型
            {"name": "ref_genome_custom", "type": "infile", "format": "ref_rna_v2.common"},  # 自定义参考基因组文件
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf,gene_structure.gff3"},  # 基因组gtf文件
            {"name": "bamlist", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "input_bam", "type": "infile", "format": "align.bwa.bam"},  # bam格式文件,排序过的
            {"name": "fq_list", "type": "string"},  # 文件里面根据单双端的不同有3列或者2列
            {"name": "des", "type": "string"},
            {"name": "des_type", "type": "string"},
        ]
        self.add_option(options)
        self.samples = []
        self.fastq_files = []
        self.step.add_steps('picard', 'gatk')  # 添加步骤
        self.count = 0
        self.ref_name = "customer_mode"
        self.annovars = []      # add by qindanhua
        self.ref_link = ''
        self.end_times = 0
        self.ref_gft = ''
        self.end_times = 0
        self.samtools_tools = []
        self.samtools_path = Config().SOFTWARE_DIR + "/bioinfo/align/samtools-1.6/samtools-1.6/samtools"
        self.split = self.add_tool('ref_rna_v2.split_ref_fasta')

    def check_options(self):
        """
        检查参数
        """
        if self.option("ref_genome") == "customer_mode" and not self.option("ref_genome_custom").is_set:
            raise OptionError("请传入自定义参考序列!", code = "23701301")

    def split_ref_run(self):
        self.split.set_options({
            "fasta": self.option("ref_genome_custom").path,
        })
        self.split.run()

    # def make_bamlist(self, fq_list):
    #     df1 = pd.read_table(self.option("bamlist").prop['path'], header=None)
    #     dir = os.path.dirname(df1.iloc[0,0])
    #     df = pd.read_table(fq_list, header=None)
    #     sample_list = df.iloc[:, 0].tolist()
    #     with open(self.work_dir + "/bamlist", "w") as fw:
    #         for i in sample_list:
    #             new_dir = os.path.join(dir, i + ".bam")
    #             fw.write(new_dir + "\n")

    def make_bamlist(self):
        old=self.option("bamlist").prop['path']
        new=self.work_dir + "/bamlist"
        if os.path.exists(new):
            os.remove(new)
        os.link(old,new)

    def samtool_run(self):
        # 这里需要先建立索引, 这个文件会产生在当前路径,之前测试的手是因为文件路径已经存在fai文件
        if not os.path.exists(self.option("ref_genome_custom").prop["path"] + ".fai"):
            cmd = "{} faidx {}".format(self.samtools_path, self.option("ref_genome_custom").prop["path"])
            self.logger.info("开始进行samtools建索引！")
            subprocess.check_call(cmd, shell=True)
        # self.make_bamlist(self.option("fq_list"))
        self.make_bamlist()
        files_list = glob2.glob(self.split.output_dir + "/split_file/*fasta*")
        for m in files_list:
            self.samtool = self.add_tool('ref_rna_v2.snp')
            self.samtool.set_options({
                "trinity_fa": self.option("ref_genome_custom").path,
                "bamlist": self.work_dir + "/bamlist",
                "bed_file": m,
            })
            self.samtools_tools.append(self.samtool)
        for j in range(len(self.samtools_tools)):
            self.samtools_tools[j].on('end', self.set_output, 'snp_indel')
        if self.samtools_tools:
            if len(self.samtools_tools) > 1:
                self.on_rely(self.samtools_tools, self.get_vcf_list)
            elif len(self.samtools_tools) == 1:
                self.samtools_tools[0].on('end', self.get_vcf_list)
        else:
            self.set_error("self.samtools_tools列表为空！", code = "23701302")
        for tool in self.samtools_tools:
            gevent.sleep(1)
            tool.run()
        self.logger.info("self.samtool is running!")


    def bcftools_vcf_run(self):
        self.bcftools_vcf = self.add_tool("ref_rna_v2.bcftool_vcf")
        # self.get_vcf_list()
        self.bcftools_vcf.set_options({
            "vcf_list": os.path.join(self.output_dir, "vcf.list"),
            # "isoform_unigene": self.option('isoform_unigene')
        })
        self.bcftools_vcf.on('end', self.set_output, 'vcf_call')
        self.bcftools_vcf.on("end", self.snp_anno)
        self.bcftools_vcf.run()

    def get_vcf_list(self):
        gevent.sleep(2)
        path = os.path.join(self.work_dir, 'snp_indel')
        files = os.listdir(path)
        vcf = os.path.join(self.output_dir, "vcf.list")
        with open(vcf, 'w') as w:
            for m in files:
                if "call" in m:
                    w.write("{}\n".format(os.path.join(path, m)))
            else:
                pass
        self.bcftools_vcf_run()

    def snp_anno(self):
        annovar = self.add_tool('ref_rna_v2.annovar')
        # 由于我们使用的tophat2和hisat2没有对sam添加文件头,所以我们在这里需要替换下最后生成的vcf的#CHROM那一行
        with open(self.bcftools_vcf.output_dir + "/pop.variant.vcf") as f, open(self.bcftools_vcf.output_dir + "/variant.vcf", "w") as w:
            for line in f:
                if not line.startswith("#CHROM"):
                    w.write(line)
                    continue
                new_line_1 = line.strip().split("\t")[0:9]
                new_line_2 = [x.split("/")[-1].split(".bam")[0] for x in line.strip().split("\t")[9:]]
                w.write("\t".join(new_line_1 + new_line_2) + "\n")
        options = {
            "ref_genome": self.option("ref_genome"),
            "input_file": self.bcftools_vcf.output_dir + "/variant.vcf",
            "combine_vcf": False,
            "des_type": self.option("des_type"),
            "des": self.option("des"),
        }
        if self.option("ref_genome") == "customer_mode":
            ref_genome_name = os.path.basename(self.option("ref_genome_custom").prop["path"])
            if ref_genome_name.endswith(".fasta") or ref_genome_name.endswith(".fa"):
                options["ref_fasta"] = self.work_dir + "/" + os.path.basename(self.option("ref_genome_custom").prop["path"])
            else:
                options["ref_fasta"] = self.work_dir + "/" + os.path.basename(self.option("ref_genome_custom").prop["path"]) + '.fa'

            if os.path.exists(options["ref_fasta"]):
                os.remove(options["ref_fasta"])
            if os.path.exists(options["ref_fasta"] + '.fai'):
                os.remove(options["ref_fasta"] + '.fai')

            os.link(self.option("ref_genome_custom").prop["path"], options["ref_fasta"])
            os.link(self.option("ref_genome_custom").prop["path"] + '.fai', options["ref_fasta"] + '.fai')
            options["ref_gtf"] = self.option("ref_gtf")
        annovar.set_options(options)
        self.annovars.append(annovar)
        self.logger.info("set output done")
        annovar.on("end", self.predeal, "annovar")
        # annovar.on('end', self.end)
        annovar.run()

    def predeal(self, event):
        obj = event["bind_object"]
        predeal = self.add_tool('ref_rna_v2.predeal_snpresults')
        options = {
            "method": "gatk",
            "snp_detail": obj.output_dir + "/snp_anno.xls",
            "snp_anno": obj.work_dir + "/snp_annotation.xls"
        }
        predeal.set_options(options)
        self.logger.info("predeal done")
        predeal.on("end", self.set_output, "predeal")
        predeal.run()

    def set_output(self, event):
        self.logger.info("set output started!!!")
        obj = event["bind_object"]
        if event['data'] == 'predeal':
            self.logger.info("llllllllllllooking for event data")
            for file in os.listdir(obj.output_dir):
                if file.endswith(".xls") or file.endswith("info"):
                    old = os.path.join(obj.output_dir, file)
                    new = os.path.join(self.output_dir, file)
                    if os.path.exists(new):
                        os.remove(new)
                    os.link(old, new)
            self.logger.info("set output done")
            self.end()
            os.remove(self.output_dir + "/vcf.list")
        if event['data'] == 'snp_indel':
            self.linkdir(obj.output_dir, 'snp_indel')
        if event['data'] == "vcf_call":
            self.linkdir(obj.output_dir, 'vcf_call')

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.link(oldfiles[i], newdir)

    def run(self):
        super(SamRnaModule, self).run()
        self.split.on('end', self.samtool_run)
        self.split_ref_run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r".", "", "结果输出目录"],
            [r"./filtered_vcf/", "文件夹", "过滤后的vcf格式的SNP位点文件结果输出目录"]
        ])
        super(SamRnaModule, self).end()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
# import glob
# import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
# from mbio.files.sequence.file_sample import FileSampleFile
import json
# from mbio.packages.gene_structure.snp_anno import snp_anno
# import subprocess
import sys
from Bio import SeqIO
import glob2
import gevent
import pandas as pd
from biocluster.config import Config


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
            {"name": "ref_genome", "type": "string"},  # 参考基因组类型
            {"name": "ref_genome_custom", "type": "infile", "format": "sequence.fasta"},  # 自定义参考基因组文件
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf,gene_structure.gff3"},  # 基因组gtf文件
            {"name": "bamlist", "type": "infile", "format": "prok_rna.common"},
            {"name": "input_bam", "type": "infile", "format": "align.bwa.bam"},  # bam格式文件,排序过的
            {"name": "fq_list", "type": "string"},  # 文件里面根据单双端的不同有3列或者2列
            {"name": "id2name", "type": "string"}
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

    def check_options(self):
        """
        检查参数
        """
        if self.option("ref_genome") == "customer_mode" and not self.option("ref_genome_custom").is_set:
            raise OptionError("请传入自定义参考序列", code = "25001001")

    def split_ref_run(self):
        split_file = "{}/split_file".format(self.work_dir)
        if os.path.exists(split_file) or os.path.isdir(split_file):
            os.system('rm -r %s' % split_file)
            os.mkdir(split_file)
        else:
            os.mkdir(split_file)

        def write_buffer(filename, data):
            with open(filename, 'wb') as target:
                for d in data:
                    target.write('{}\t{}\t{}\n'.format(d.id, 0, len(d.seq)))
        buffers = list()
        i = 0
        file_size = os.path.getsize(self.option("ref_genome_custom").prop['path'])/float(1024*1024*1024)
        file_size = round(file_size, 2)
        # 按照大约50个tool投递
        line_limit = 800000*file_size
        for seq_record in SeqIO.parse(self.option("ref_genome_custom").prop['path'], 'fasta'):
            i += 1
            filename = split_file + '/fasta_{}'.format(i)
            if len(seq_record.seq) >= line_limit:
                write_buffer(filename, [seq_record])
                continue
            buffers.append(seq_record)
            if sum(len(i.seq) for i in buffers) >= line_limit:
                write_buffer(filename, buffers)
                del buffers[:]
        filename = split_file + '/fasta_{}'.format(i+1)
        if buffers:
            write_buffer(filename, buffers)
            del buffers[:]

    def make_bamlist(self, fq_list):
        df1 = pd.read_table(self.option("bamlist").prop['path'], header=None)
        dir = os.path.dirname(df1.iloc[0,0])
        df = pd.read_table(fq_list, header=None)
        sample_list = df.iloc[:, 0].tolist()
        with open(self.work_dir + "/bamlist", "w") as fw:
            for i in sample_list:
                new_dir = os.path.join(dir, str(i) + ".bam")
                fw.write(new_dir + "\n")

    def samtool_run(self):
        self.split_ref_run()
        # 这里需要先建立索引, 这个文件会产生在当前路径,之前测试的手是因为文件路径已经存在fai文件
        # 在module当中禁止调用subprocess，会启动双进程
        # if not os.path.exists(self.option("ref_genome_custom").prop["path"] + ".fai"):
        #     cmd = "{} faidx {}".format(self.samtools_path, self.option("ref_genome_custom").prop["path"])
        #     self.logger.info("开始进行samtools建索引！")
        #     subprocess.check_call(cmd, shell=True)
        self.make_bamlist(self.option("fq_list"))
        files_list = glob2.glob(self.work_dir + "/split_file/*fasta*")
        for m in files_list:
            self.samtool = self.add_tool('prok_rna.snp')
            self.samtool.set_options({
                "trinity_fa": self.option("ref_genome_custom"),
                "bamlist": self.work_dir + "/bamlist",
                "bed_file": m,
            })
            self.samtools_tools.append(self.samtool)
        for j in range(len(self.samtools_tools)):
            self.samtools_tools[j].on('end', self.set_output, 'snp_indel')
        if self.samtools_tools:
            if len(self.samtools_tools) > 1:
                self.on_rely(self.samtools_tools, self.bcftools_vcf_run)
            elif len(self.samtools_tools) == 1:
                self.samtools_tools[0].on('end', self.bcftools_vcf_run)
        else:
            self.set_error("self.samtools_tools列表为空", code = "25001002")
        for tool in self.samtools_tools:
            gevent.sleep(1)
            tool.run()
        self.logger.info("self.samtool is running!")


    def bcftools_vcf_run(self):
        self.bcftools_vcf = self.add_tool("prok_rna.bcftool_vcf")
        self.get_vcf_list()
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

    def snp_anno(self):
        annovar = self.add_tool('prok_rna.annovar')
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
            "id2name": self.option("id2name")
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
        annovar.on("end", self.set_output, 'annovar')
        annovar.on('end', self.end)
        annovar.run()

    def set_output(self, event):
        self.logger.info("set output started!!!")
        obj = event["bind_object"]
        if event['data'] == 'annovar':
            output_name = self.output_dir + "/" + "snp_anno.xls"
            self.logger.info(output_name)
            if os.path.exists(output_name):
                os.remove(output_name)
            if os.path.exists(self.work_dir + "/snp_annotation.xls"):
                os.remove(self.work_dir + "/snp_annotation.xls")
            os.link(obj.output_dir + "/snp_anno.xls", output_name)
            os.link(obj.work_dir + "/snp_annotation.xls", self.work_dir + "/snp_annotation.xls")
            # os.link(obj.work_dir + "/snp_transition_tranversion_statistics.xls", self.output_dir + "/snp_transition_tranversion_statistics.xls")
            # os.link(obj.work_dir + "/snp_freq_statistics.xls", self.output_dir + "/snp_freq_statistics.xls")
            # os.link(obj.work_dir + "/snp_depth_statistics.xls", self.output_dir + "/snp_depth_statistics.xls")
            # os.link(obj.work_dir + "/snp_position_distribution.xls", self.output_dir + "/snp_position_distribution.xls")
            # os.link(obj.work_dir + "/indel_position_distribution.xls", self.output_dir + "/indel_position_distribution.xls")
            self.logger.info("set output done")
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
        self.samtool_run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r".", "", "结果输出目录"],
            [r"./filtered_vcf/", "文件夹", "过滤后的vcf格式的SNP位点文件结果输出目录"]
        ])
        super(SamRnaModule, self).end()

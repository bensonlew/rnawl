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


class SnpRnaModule(Module):
    """
    star:序列比对
    picard:处理比对结果sam文件
    gatk：snp calling软件
    version 1.0
    author: chenyanyan
    last_modify: 2017.01.03 by qindanhua
    """
    def __init__(self, work_id):
        super(SnpRnaModule, self).__init__(work_id)
        options = [
            {"name": "ref_genome", "type": "string"},  # 参考基因组类型
            {"name": "ref_genome_custom", "type": "infile", "format": "sequence.fasta"},  # 自定义参考基因组文件
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf,gene_structure.gff3"},  # 基因组gtf文件
            # {"name": "in_bam", "type": "infile", "format": "align.bwa.bam_dir"},
            # {"name": "readFilesIN", "type": "infile", "format": "sequence.fastq"},  # 用于比对的单端序列文件
            # {"name": "readFilesIN1", "type": "infile", "format": "sequence.fastq, sequence.fasta"},  # 双端序列←
            # {"name": "readFilesIN2", "type": "infile", "format": "sequence.fastq, sequence.fasta"},  # 双端序列右
            {"name": "in_bam", "type": "infile", "format": "align.bwa.bam_dir"},  # bam格式文件
            {"name": "seq_method", "type": "string"},  # 比对方式
            {"name": "input_bam", "type": "infile", "format": "align.bwa.bam"}  # bam格式文件,排序过的
        ]
        self.add_option(options)
        self.samples = []
        self.fastq_files = []
        self.mapping_tools = []  # star的tools
        self.picards = []
        self.gatks = []
        self.step.add_steps('picard', 'gatk')  # 添加步骤
        self.count = 0
        self.ref_name = "customer_mode"
        self.annovars = []      # add by qindanhua
        self.ref_link = ''
        self.end_times = 0
        self.ref_gft = ''

    def check_options(self):
        """
        检查参数
        """
        if self.option("ref_genome") == "customer_mode" and not self.option("ref_genome_custom").is_set:
            raise OptionError("请传入自定义参考序列!")

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def picard_run(self):
        for i in sorted(os.listdir(self.option("in_bam").prop["path"])):
            self.samples.append(i[:-4])
            f_path = os.path.join(self.option("in_bam").prop["path"], i)
            self.logger.info(f_path)  # 打印出f_path的信息，是上一步输出文件的路径
            picard = self.add_tool('gene_structure.picard_rna')
            picard.set_options({
                "in_bam": f_path
            })
            self.picards.append(picard)
            picard.on("end", self.gatk_run, i[:-4])  # event["data"]传为样本名称
        if len(self.picards) == 1:
            self.picards[0].on("end", self.finish_update, "picard")
        else:
            self.on_rely(self.picards, self.finish_update, 'picard')
        for picard in self.picards:
            picard.run()

    def gatk_run(self, event):
        obj = event["bind_object"]
        picard_output = os.listdir(obj.output_dir)
        for i in picard_output:
            if i.endswith(".bam"):
                f_path = os.path.join(obj.output_dir, i)
                self.logger.info(f_path)
        gatk = self.add_tool('gene_structure.gatk')
        self.gatks.append(gatk)
        self.logger.info("因为参考数据库中未含有gatk所需dict，gatk使用自定义模式进行")
        #
        # if self.option("ref_genome") != "customer_mode":
        #     with open(self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.json", "r") as f:
        #         dict = json.loads(f.read())
        #         ref_fasta = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/" + dict[self.option("ref_genome")]["dna_fa"]
        gatk.set_options({
            "ref_fa": self.option("ref_genome_custom"),
            "input_bam": f_path,
            "ref_genome": "customer_mode"
        })
        # else:
        #     self.ref_name = self.option("ref_genome")
        #     gatk.set_options({
        #         "ref_genome": self.ref_name,
        #         "input_bam": f_path
        #     })
        gatk.on("end", self.finish_update, 'gatk')
        gatk.on("end", self.snp_anno, event['data'])
        self.logger.info("gatk is running!")
        gatk.run()

    # add by qindanhua  20170103
    def snp_anno(self, event):
        obj = event["bind_object"]
        self.logger.info(event['data'])
        gatk_output = os.listdir(obj.output_dir)
        vcf_path = ""
        for i in gatk_output:
            if i.endswith(".vcf"):
                vcf_path = os.path.join(obj.output_dir, i)
        self.logger.info(vcf_path)
        if not os.path.exists(self.work_dir + "/output_vcf"):
            os.mkdir(self.work_dir + "/output_vcf")
        self.logger.info(self.work_dir + "/output_vcf/" + event['data'] + ".vcf")
        output_vcf = self.work_dir + "/output_vcf/" + event['data'] + ".vcf"
        if not os.path.exists(output_vcf):
            os.link(vcf_path, output_vcf)
        self.end_times += 1
        if self.end_times == len(self.samples):
            annovar = self.add_tool('gene_structure.annovar')
            options = {
                "ref_genome": self.ref_name,
                "input_file": self.work_dir + "/output_vcf",
                "combine_vcf": True
                # "ref_gtf": self.option("ref_gtf").prop["path"]
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
            annovar.on("end", self.set_output, "annovar_" + event['data'])
            annovar.run()

    def set_output(self, event):
        self.logger.info("set output started!!!")
        obj = event["bind_object"]
        if event['data'][:7] == 'annovar':
            # output_name = self.output_dir + "/" + "{}.snp_anno.xls".format(event['data'][8:])
            output_name = self.output_dir + "/" + "snp_anno.xls"
            self.logger.info("llllllllllllooking for event data")
            self.logger.info(event['data'][8:])
            self.logger.info(output_name)
            if os.path.exists(output_name):
                os.remove(output_name)
            # snp_freq_stat(obj.work_dir + "snp.vcf", obj.output_dir + "/snp_anno.xls", output_name)
            os.link(obj.output_dir + "/snp_anno.xls", output_name)
            # self.end_times += 1
            # if self.end_times == len(self.samples):
            self.logger.info("set output done")
            self.end()
        # self.logger.info("set output done")

    def linkdir(self, dirpath, dirname, output_dir):
        files = os.listdir(dirpath)
        newdir = os.path.join(output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in files]
        newfiles = [os.path.join(newdir, i) for i in files]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(files)):
            if os.path.exists(newfiles[i]):
                os.remove(newfiles[i])
            os.link(oldfiles[i], newfiles[i])

    def run(self):
        super(SnpRnaModule, self).run()
        self.picard_run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r".", "", "结果输出目录"],
            [r"./filtered_vcf/", "文件夹", "过滤后的vcf格式的SNP位点文件结果输出目录"]
        ])
        super(SnpRnaModule, self).end()

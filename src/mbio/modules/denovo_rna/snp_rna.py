#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import json
import unittest
from biocluster.config import Config

import gevent.subprocess as subprocess



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
            {"name": "ref_genome_custom", "type": "infile", "format": "ref_rna_v2.common"},  # 自定义参考基因组文件
            {"name": "ref_gtf", "type": "infile", "format": "ref_rna_v2.common"},  # 基因组gtf文件
            {"name": "in_bam", "type": "infile", "format": "align.bwa.bam_dir"},  # bam格式文件
            {"name": "input_bam", "type": "infile", "format": "align.bwa.bam"},  # bam格式文件,排序过的,
            {"name": "des", "type": "string"},
            {"name": "des_type", "type": "string"},
        ]
        self.add_option(options)
        self.samples = []
        self.fastq_files = []
        self.mapping_tools = []  # star的tools
        self.picards = []
        self.gatks = []
        self.step.add_steps('picard', 'gatk',"vcffilter")  # 添加步骤
        self.count = 0
        self.ref_name = "customer_mode"
        self.annovars = []      # add by qindanhua
        self.ref_link = ''
        self.end_times = 0
        self.ref_gft = ''
        self.ref_fasta = ""
        self.ref_index1 = ""
        self.ref_index2 = ""
        self.samtools_path = Config().SOFTWARE_DIR + "/bioinfo/ref_rna_v2/miniconda2/bin/"
        self.picard_path = Config().SOFTWARE_DIR + "/bioinfo/gene-structure/"
        self.java_path = Config().SOFTWARE_DIR + "/program/sun_jdk1.8.0/bin/java"

    def check_options(self):
        """
        检查参数
        """
        if self.option("ref_genome") == "customer_mode" and not self.option("ref_genome_custom").is_set:
            raise OptionError("请传入自定义参考序列!", code = "23701401")

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run(self):
        super(SnpRnaModule, self).run()
        self.genome_split()
        self.picard_run()

    def genome_split(self):
        """
        将所有样本的bam总文件，进行分割，然后进行并行call snp
        :return:
        """
        ref_genome_name = os.path.basename(self.option("ref_genome_custom").prop["path"])
        if ref_genome_name.endswith(".fasta") or ref_genome_name.endswith(".fa"):
            self.ref_fasta = self.work_dir + "/" + os.path.basename(self.option("ref_genome_custom").prop["path"])
        else:
            self.ref_fasta = self.work_dir + "/" + os.path.basename(self.option("ref_genome_custom").prop["path"]) + '.fa'

        if os.path.exists(self.ref_fasta):
            os.remove(self.ref_fasta)
        if os.path.exists(self.ref_fasta + '.fai'):
            os.remove(self.ref_fasta + '.fai')
        if os.path.exists(self.work_dir+"/"+".".join(os.path.basename(self.option("ref_genome_custom").prop["path"]).split(".")[:-1])+ ".dict" ):
            os.remove(self.work_dir+"/"+".".join(os.path.basename(self.option("ref_genome_custom").prop["path"]).split(".")[:-1])+ ".dict" )
        os.link(self.option("ref_genome_custom").prop["path"], self.ref_fasta)
        self.ref_index1 = self.ref_fasta + ".fai"
        self.ref_index2 = self.work_dir+"/"+".".join(os.path.basename(self.option("ref_genome_custom").prop["path"]).split(".")[:-1])+ ".dict"
        if not os.path.exists(self.option("ref_genome_custom").prop["path"] + ".fai"):
            cmd = "{} faidx {}".format(self.samtools_path, self.option("ref_genome_custom").prop["path"])
            self.logger.info("开始进行samtools建索引！")
            subprocess.check_call(cmd, shell=True)
            os.link(self.option("ref_genome_custom").prop["path"] + '.fai', self.ref_index1)
        else:
            os.link(self.option("ref_genome_custom").prop["path"] + '.fai', self.ref_index1)
        if not os.path.exists(os.path.dirname(self.option("ref_genome_custom").prop["path"]) + "/" + ".".join(os.path.basename(self.option("ref_genome_custom").prop["path"]).split(".")[:-1]) + ".dict"):
            dirname = os.path.split(self.option("ref_genome_custom").prop["path"])[-1]
            dictname = os.path.splitext(dirname)[0] + ".dict"
            self.dict(self.ref_fasta, os.path.join(os.path.split(self.option("ref_genome_custom").prop["path"])[0],dictname))
            os.link(os.path.join(os.path.split(self.option("ref_genome_custom").prop["path"])[0], dictname), self.ref_index2)
        else:
            os.link(os.path.dirname(self.option("ref_genome_custom").prop["path"]) + "/" + ".".join(os.path.basename(self.option("ref_genome_custom").prop["path"]).split(".")[:-1]) + ".dict", self.ref_index2)

    def dict(self, ref_fasta, dict_name):
        """
        使用picard对参考基因组构建字典
        """

        cmd = "{} -jar {}picard.jar CreateSequenceDictionary R={} O={}" \
            .format(self.java_path,self.picard_path, ref_fasta, dict_name)
        if os.path.exists(dict_name):
            os.remove(dict_name)
        print cmd
        self.logger.info("开始用picard对参考基因组构建字典")
        subprocess.check_call(cmd, shell=True)
        self.logger.info("参考基因组构建dict done!")

    def picard_run(self):
        for i in os.listdir(self.option("in_bam").prop["path"]):
            self.samples.append(i[:-4])
            f_path = os.path.join(self.option("in_bam").prop["path"], i)
            self.logger.info(f_path)  # 打印出f_path的信息，是上一步输出文件的路径
            picard = self.add_tool('ref_rna_v3.picard_rna')
            picard.set_options({
                "in_bam": f_path,
                "in_reference": self.option("ref_genome_custom")
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
        gatk = self.add_tool('ref_rna_v2.gatk')
        self.gatks.append(gatk)
        self.logger.info("因为参考数据库中未含有gatk所需dict，gatk使用自定义模式进行")
        gatk.set_options({
            "ref_fa": self.option("ref_genome_custom"),
            "input_bam": f_path,
            "ref_genome": "customer_mode"
        })
        gatk.on("end", self.finish_update, 'gatk')
        gatk.on("end", self.vcf_filter_run, event['data'])
        self.logger.info("gatk is running!")
        gatk.run()

    def vcf_filter_run(self,event):
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
            self.vcffilter = self.add_tool('ref_rna_v2.vcf_filter_gatk')
            options = {
                "ref_genome": self.ref_name,
                "input_file": self.work_dir + "/output_vcf",
                "combine_vcf": True,
                "ref_fasta": self.option("ref_genome_custom"),
            }
            self.vcffilter.set_options(options)
            self.vcffilter.on("end", self.finish_update, 'vcffilter')
            self.vcffilter.on("end",self.snp_anno,"predeal")
            self.logger.info("vcf_filter is running!")
            self.vcffilter.run()

    # add by qindanhua  20170103
    def snp_anno(self, event):
        obj = event["bind_object"]
        self.logger.info(event['data'])
        vcf_output = os.path.join(obj.output_dir,"final.vcf")
        # vcf_path = ""
        # for i in gatk_output:
        #     if i.endswith(".vcf"):
        #         vcf_path = os.path.join(obj.output_dir, i)
        # self.logger.info(vcf_path)
        # if not os.path.exists(self.work_dir + "/output_vcf"):
        #     os.mkdir(self.work_dir + "/output_vcf")
        # self.logger.info(self.work_dir + "/output_vcf/" + event['data'] + ".vcf")
        # output_vcf = self.work_dir + "/output_vcf/" + event['data'] + ".vcf"
        # if not os.path.exists(output_vcf):
        #     os.link(vcf_path, output_vcf)
        # self.end_times += 1
        # if self.end_times == len(self.samples):
        annovar = self.add_tool('ref_rna_v3.annovar')
        options = {
                "ref_genome": self.ref_name,
                "input_file": vcf_output,
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
        annovar.on("end", self.predeal, "annovar")
        self.annovars.append(annovar)
        self.logger.info("annovar done")
        annovar.run()

    def predeal(self, event):
        obj = event["bind_object"]
        predeal = self.add_tool('ref_rna_v3.predeal_snpresults')
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

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r".", "", "结果输出目录"],
            [r"./filtered_vcf/", "文件夹", "过滤后的vcf格式的SNP位点文件结果输出目录"]
        ])
        super(SnpRnaModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "snp_test_gatk" + str(random.randint(1, 10000))+"-xxx",
            "type": "module",
            "name": "denovo_rna.snp_rna",
            "instant": False,
            "options": dict(
                ref_genome="customer_mode",
                ref_genome_custom="/mnt/ilustre/users/sanger-dev/workspace/20190416/Denovorna_tsg_33857/DenovoAssemble2/output/Trinity.filter.unigene.fasta",
                #ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/gtf/Mus_musculus.GRCm38.89.gtf",
                in_bam="/mnt/ilustre/users/sanger-dev/workspace/20190416/Denovorna_tsg_33857/Snp/bam_folder",
                #des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/biomart/Mus_musculus.GRCm38.biomart_gene.txt",
                #des_type="type1",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()


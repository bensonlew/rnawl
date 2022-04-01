# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last_modify:20180404

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import unittest
import gevent.subprocess as subprocess
from biocluster.config import Config


class CallSnpIndelSentieonModule(Module):
    """
    用于对bam文件进行call突变位点，该模块中包含了三种方式的call 突变位点，分别是gatks，samtools，freebayes，其中进行call的
    时候是要将所有的bam文件合并在一起，然后将合并后的基因组分成25份进行并行计算，计算完成后再合并到一起
    last modified by fwy @ 20190620 改写为ref_rna_v2专用
    """
    def __init__(self, work_id):
        super(CallSnpIndelSentieonModule, self).__init__(work_id)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "bam_list", "type": "string"},  # 如果是gatk 样本+建库类型+bam文件，否则样本bam列表
            #{"name": "ref_dict", "type": "string"},  # 参考基因组的配置文件，里面可以解析出有多少条染色体
            {"name": "ref_gtf", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "call_type", "type": "string", "default": "sentieon"},  # call snp的方式
            {"name": "des", "type": "string"},
            {"name": "des_type", "type": "string"},
            {"name": "analysis_format", "type": "string", "default": "bam"},
            {"name": "align_method", "type": "string", "default": "hisat"},
        ]
        self.add_option(options)
        self.file_list = []
        self.samtools_call_tools = []
        self.freebayes_call_tools = []
        self.haplotype_tools = []
        self.gvcf_typing_tools = []
        #self.bcftools_vcf = self.add_tool("wgs.bcftool_vcf")
        self.gvcf_typing_v2 = None
        self.ref_fasta=""
        self.ref_index1=""
        self.ref_index2=""
        self.samtools_path = Config().SOFTWARE_DIR + "/bioinfo/ref_rna_v2/miniconda2/bin/"
        self.picard_path = Config().SOFTWARE_DIR + "/bioinfo/gene-structure/"
        self.java_path=Config().SOFTWARE_DIR+"/program/sun_jdk1.8.0/bin/java"
        self.annovars = []

    def check_options(self):
        if not self.option("ref_fasta").is_set:
            raise OptionError("缺少ref_fasta参数", code="25600305")
        if not self.option("bam_list"):
            raise OptionError("缺少bam_list参数", code="25600306")
        #if not self.option("ref_dict"):
            #raise OptionError("缺少ref_dict参数", code="24500403")
        if self.option("call_type"):
            if self.option("call_type") not in ["samtools", "gatk", "freebayes", "sentieon"]:
                raise OptionError("call_type方式不合法，必须是samtools, gatk, freebayes, sentieon", code="25600307")
        else:
            raise OptionError("缺少call_type参数", code="25600308")
        ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        ref_file_name = os.path.basename(self.option("ref_fasta").prop['path']).split(".fa")[0]
        #if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
                #or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
            #raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="24500406")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def genome_split(self):
        """
        将所有样本的bam总文件，进行分割，然后进行并行call snp
        :return:
        """
        ref_genome_name = os.path.basename(self.option("ref_fasta").prop["path"])
        if ref_genome_name.endswith(".fasta") or ref_genome_name.endswith(".fa"):
            self.ref_fasta = self.work_dir + "/" + os.path.basename(self.option("ref_fasta").prop["path"])
        else:
            self.ref_fasta = self.work_dir + "/" + os.path.basename(self.option("ref_fasta").prop["path"]) + '.fa'

        if os.path.exists(self.ref_fasta):
            os.remove(self.ref_fasta)
        if os.path.exists(self.ref_fasta + '.fai'):
            os.remove(self.ref_fasta + '.fai')
        if os.path.exists(self.work_dir+"/"+".".join(os.path.basename(self.option("ref_fasta").prop["path"]).split(".")[:-1])+ ".dict" ):
            os.remove(self.work_dir+"/"+".".join(os.path.basename(self.option("ref_fasta").prop["path"]).split(".")[:-1])+ ".dict" )
        os.link(self.option("ref_fasta").prop["path"], self.ref_fasta)
        self.ref_index1 = self.ref_fasta + ".fai"
        self.ref_index2 = self.work_dir+"/"+".".join(os.path.basename(self.option("ref_fasta").prop["path"]).split(".")[:-1])+ ".dict"
        if not os.path.exists(self.option("ref_fasta").prop["path"] + ".fai"):
            cmd = "{} faidx {}".format(self.samtools_path, self.option("ref_fasta").prop["path"])
            self.logger.info("开始进行samtools建索引！")
            subprocess.check_call(cmd, shell=True)
            os.link(self.option("ref_fasta").prop["path"] + '.fai', self.ref_index1)
        else:
            os.link(self.option("ref_fasta").prop["path"] + '.fai', self.ref_index1)
        # if not os.path.exists(os.path.dirname(self.option("ref_fasta").prop["path"]) + "/" + ".".join(os.path.basename(self.option("ref_fasta").prop["path"]).split(".")[:-1]) + ".dict"):
        #     dirname = os.path.split(self.option("ref_fasta").prop["path"])[-1]
        #     dictname = os.path.splitext(dirname)[0] + ".dict"
        #     self.dict(self.ref_fasta, os.path.join(os.path.split(self.option("ref_fasta").prop["path"])[0],dictname))
        #     os.link(os.path.join(os.path.split(self.option("ref_fasta").prop["path"])[0], dictname), self.ref_index2)
        # else:
        #     os.link(os.path.dirname(self.option("ref_fasta").prop["path"]) + "/" + ".".join(os.path.basename(self.option("ref_fasta").prop["path"]).split(".")[:-1]) + ".dict", self.ref_index2)
        # n = 0
        # split_file = "{}/split_file".format(self.work_dir)
        # if os.path.exists(split_file) or os.path.isdir(split_file):
        #     os.system('rm -r %s' % split_file)
        #     os.mkdir(split_file)
        # else:
        #     os.mkdir(split_file)
        # with open(self.ref_index2, "r") as r:
        #     data = r.readlines()[1:]
        #     for line in data:
        #         temp = line.strip().split("\t")
        #         chr_name = temp[1].strip().split(":")[1]
        #         length = temp[2].strip().split(':')[1]
        #         n += 1
        #         file_name = n % 25
        #         if self.option("call_type") == "gatk":
        #             with open(split_file + "/{}.intervals".format(file_name), 'a+') as w:
        #                 w.write(chr_name + "\n")
        #         else:
        #             with open(split_file + "/{}.bed".format(file_name), 'a+') as w:
        #                 w.write(chr_name + "\t" + "0\t" + length + "\n")
        # self.get_bed_intervals_file(split_file)

    def get_bed_intervals_file(self, split_file):
        """
        获取分割后的bed文件或者intervals文件的个数
        :return:
        """
        files = os.listdir(split_file)
        if len(files) == 0:
            self.set_error("基因组分割后的列表为空，不能进行后续分析下，请检查ref.dict文件！", code="25600307")
        else:
            for m in files:
                self.file_list.append(os.path.join(split_file, m))
            self.logger.info("file_list为{}".format(self.file_list))

    def samtools_call_run(self):
        """
        按照work_dir中的split_file有多少个文件就运行多少个tool
        :return:
        """
        self.set_bam_list()
        n = 0
        for m in self.file_list:
            samtools_call = self.add_tool("wgs.samtools_call")
            samtools_call.set_options({
                "ref_fasta": self.option("ref_fasta").prop["path"],
                "bed_file": m,
                "bam_list": "{}/bam.list".format(self.work_dir)
            })
            self.samtools_call_tools.append(samtools_call)
            n += 1
        for j in range(len(self.samtools_call_tools)):
            self.samtools_call_tools[j].on('end', self.set_output, 'snp_indel')
        if self.samtools_call_tools:
            if len(self.samtools_call_tools) > 1:
                self.on_rely(self.samtools_call_tools, self.bcftools_vcf_run)
            elif len(self.samtools_call_tools) == 1:
                self.samtools_call_tools[0].on('end', self.bcftools_vcf_run)
            for tool in self.samtools_call_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("samtools_call_tools列表为空！", code="25600308")

    def freebayes_call_run(self):
        self.set_bam_list()
        n = 0
        for m in self.file_list:
            freebayes_call = self.add_tool("wgs.freebayes_call")
            freebayes_call.set_options({
                "ref_fasta": self.option("ref_fasta").prop["path"],
                "bed_file": m,
                "bam_list": "{}/bam.list".format(self.work_dir)
            })
            self.freebayes_call_tools.append(freebayes_call)
            n += 1
        for j in range(len(self.freebayes_call_tools)):
            self.freebayes_call_tools[j].on('end', self.set_output, 'snp_indel')
        if self.freebayes_call_tools:
            if len(self.freebayes_call_tools) > 1:
                self.on_rely(self.freebayes_call_tools, self.bcftools_vcf_run)
            elif len(self.freebayes_call_tools) == 1:
                self.freebayes_call_tools[0].on('end', self.bcftools_vcf_run)
            for tool in self.freebayes_call_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("freebayes_call_tools列表为空！", code="25600309")

    def haplotype_run(self):
        sample_info = self.get_sample_bam()
        n = 0
        for keys in sample_info.keys():
            haplotype = self.add_tool("wgs.haplotype")
            haplotype.set_options({
                "ref_fasta": self.option("ref_fasta").prop["path"],
                "sample_bam": sample_info[keys],
                "sample_id": keys
            })
            self.haplotype_tools.append(haplotype)
            n += 1
        for j in range(len(self.haplotype_tools)):
            self.haplotype_tools[j].on('end', self.set_output, 'haplotype')
        if self.haplotype_tools:
            if len(self.haplotype_tools) > 1:
                self.on_rely(self.haplotype_tools, self.gvcf_typing_run)
            elif len(self.haplotype_tools) == 1:
                self.haplotype_tools[0].on('end', self.gvcf_typing_run)
            for tool in self.haplotype_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("haplotype_tools列表为空！", code="25600310")

    def gvcf_typing_run(self):
        n = 0
        for m in self.file_list:
            gvcf_typing = self.add_tool("wgs.gvcf_typing")
            gvcf_typing.set_options({
                "ref_fasta": self.option("ref_fasta").prop["path"],
                "gvcf_path": self.output_dir + "/haplotype",
                "intervals_file": m
            })
            self.gvcf_typing_tools.append(gvcf_typing)
            n += 1
        for j in range(len(self.gvcf_typing_tools)):
            self.gvcf_typing_tools[j].on('end', self.set_output, 'snp_indel')
        if self.gvcf_typing_tools:
            if len(self.gvcf_typing_tools) > 1:
                self.on_rely(self.gvcf_typing_tools, self.bcftools_vcf_run)
            elif len(self.gvcf_typing_tools) == 1:
                self.gvcf_typing_tools[0].on('end', self.bcftools_vcf_run)
            for tool in self.gvcf_typing_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("gvcf_typing_tools列表为空！", code="25600311")

    def haplotype_v2_run(self):
        sample_info = self.get_sample_bam()
        n = 0
        for keys in sample_info.keys():
            haplotype = self.add_module("ref_rna_v3.pre4sentieon")
            haplotype.set_options({
                "ref_fasta": self.option("ref_fasta").prop["path"],
                "in_bam": sample_info[keys],
                "call_type": "sentieon",
                "analysis_format":self.option("analysis_format"),
                "align_method":self.option("align_method")
            })
            self.haplotype_tools.append(haplotype)
            n += 1
        for j in range(len(self.haplotype_tools)):
            self.haplotype_tools[j].on('end', self.set_output, 'haplotype')
        if self.haplotype_tools:
            if len(self.haplotype_tools) > 1:
                self.on_rely(self.haplotype_tools, self.gvcf_typing_v2_run)
            elif len(self.haplotype_tools) == 1:
                self.haplotype_tools[0].on('end', self.gvcf_typing_v2_run)
        else:
            self.set_error("haplotype_tools列表为空！", code="25600312")
        for tool in self.haplotype_tools:
            gevent.sleep(1)
            tool.run()

    def gvcf_typing_v2_run(self):
        """
        DE1_10  /mnt/ilustre/users/caixia.tian/newmdt/Project/MJ20180723036-gaojun/03.ref/08.haplotype/DE1_10.g.vcf
        :return:
        """
        self.make_gvcf_list()
        self.gvcf_typing_v2 = self.add_tool("ref_rna_v2.gvcf_typing_v2")
        # self.gvcf_typing_v2 = self.add_tool("ref_rna_v3.snp.gvcf_typing_v2")
        self.gvcf_typing_v2.set_options({
            "vcf_list": self.output_dir + "/haplotype/gvcf.list",
            "fa_file": self.option("ref_fasta").prop["path"]
        })
        self.gvcf_typing_v2.on('end', self.set_output, 'vcf_call')
        self.gvcf_typing_v2.on('end', self.vcf_filter_run)
        self.gvcf_typing_v2.run()

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

    def vcf_filter_run(self):
        self.vcffilter = self.add_tool('ref_rna_v3.snp.vcf_filter_gatk')
        options = {
            "ref_genome": "customer_mode",
            "input_file": self.gvcf_typing_v2.output_dir + "/pop.variant.vcf",
            "combine_vcf": False,
            "ref_fasta": self.ref_fasta,
        }
        self.vcffilter.set_options(options)
        sample_number = len(open(self.option("bam_list"), 'r').readlines())
        if sample_number >= 48:
            self.vcffilter.on("end", self.split_vcf)
        else:
            self.vcffilter.on("end", self.snp_anno)
        self.logger.info("vcf_filter is running!")
        self.vcffilter.run()

    def split_vcf(self):
        self.split_vcf = self.add_tool('ref_rna_v3.large.split_vcf')
        options = {
            "fasta": self.ref_fasta,
            "vcf": self.vcffilter.output_dir + "/final.vcf"
        }
        self.split_vcf.set_options(options)
        self.split_vcf.on("end", self.run_anno)
        self.split_vcf.run()

    def run_anno(self):
        for file in os.listdir(self.split_vcf.output_dir):
            annovar = self.add_tool('ref_rna_v3.annovar')
            options = {
                "ref_genome": "customer_mode",
                "ref_fasta": self.ref_fasta,
                "input_file": os.path.join(self.split_vcf.output_dir, file),
                "combine_vcf": False,
                "des_type": self.option("des_type"),
                "des": self.option("des"),
                "ref_gtf": self.option("ref_gtf"),
            }
            annovar.set_options(options)
            self.annovars.append(annovar)
            # annovar.run()
        if len(self.annovars) > 1:
            self.on_rely(self.annovars, self.merge_anno)
        elif len(self.annovars) == 1:
            self.annovars[0].on('end', self.merge_anno)
        for annovar in self.annovars:
            annovar.run()


    def merge_anno(self):
        self.cats = list()
        snp_anno_list = os.path.join(self.work_dir, "snp_anno_list.txt")
        with open(snp_anno_list, "w") as f1:
            for tool in self.annovars:
                f1.write(os.path.join(tool.work_dir, "snp_anno.xls") + "\n")
        merge_anno = self.add_tool('ref_rna_v3.large.cat')
        options = {
            "input_list": snp_anno_list,
            "out_file": os.path.join(self.work_dir, "snp_anno.xls"),
        }
        merge_anno.set_options(options)
        merge_anno.run()

        snp_annotation_list = os.path.join(self.work_dir, "snp_annotation_list.txt")
        with open(snp_annotation_list, "w") as f2:
            for tool in self.annovars:
                f2.write(os.path.join(tool.work_dir, "snp_annotation.xls") + "\n")
        merge_annotation = self.add_tool('ref_rna_v3.large.cat')
        options = {
            "input_list": snp_annotation_list,
            "out_file": os.path.join(self.work_dir, "snp_annotation.xls"),
        }
        merge_annotation.set_options(options)
        merge_annotation.run()
        self.on_rely([merge_anno, merge_annotation], self.predeal1)

    def snp_anno(self):
        annovar = self.add_tool('ref_rna_v3.annovar')
        options = {
            "ref_genome": "customer_mode",
            "ref_fasta":self.ref_fasta,
            "input_file": self.vcffilter.output_dir + "/final.vcf",
            "combine_vcf": False,
            "des_type": self.option("des_type"),
            "des": self.option("des"),
            "ref_gtf":self.option("ref_gtf"),
        }
        annovar.set_options(options)
        self.annovars.append(annovar)
        self.logger.info("set output done")
        annovar.on("end", self.predeal, "annovar")
        # annovar.on('end', self.end)
        annovar.run()

    def predeal1(self):
        predeal = self.add_tool('ref_rna_v3.predeal_snpresults2')
        options = {
            "method": "gatk",
            "snp_detail": os.path.join(self.work_dir, "snp_anno.xls"),
            "snp_anno": os.path.join(self.work_dir, "snp_annotation.xls")
        }
        predeal.set_options(options)
        self.logger.info("predeal done")
        predeal.on("end", self.set_output, "predeal")
        predeal.run()

    def predeal(self, event):
        obj = event["bind_object"]
        predeal = self.add_tool('ref_rna_v3.predeal_snpresults2')
        options = {
            "method": "gatk",
            "snp_detail": obj.output_dir + "/snp_anno.xls",
            "snp_anno": obj.work_dir + "/snp_annotation.xls"
        }
        predeal.set_options(options)
        self.logger.info("predeal done")
        predeal.on("end", self.set_output, "predeal")
        predeal.run()


    def make_gvcf_list(self):
        """
        DE1_10  /mnt/ilustre/users/caixia.tian/newmdt/Project/MJ20180723036-gaojun/03.ref/08.haplotype/DE1_10.g.vcf
        :return:
        """
        with open(self.output_dir + "/haplotype/gvcf.list", 'w') as w:
            for m in os.listdir(self.output_dir + "/haplotype"):
                n = re.match('(.*)\.g\.vcf$', m)
                if n:
                    w.write('{}\t{}\n'.format(n.group(1), os.path.join(self.output_dir + "/haplotype", m)))

    def get_sample_bam(self):
        """
        生成样本与bam文件对应的字典文件，这里只是gatk call snp需要， bam文件来源于04.bam-sort或者05.bam-mkdup文件夹
        :return:
        """
        sample_info = {}
        with open(self.option("bam_list"), "r") as r:
            data = r.readlines()
            for line in data:
                temp = line.strip().split("\t")
                if len(temp) == 3:
                    if temp[1].upper() not in ['WGS', 'RAD', 'GBS', 'WES']:
                        self.logger.warn("{}建库类型不合法!".format(temp[1]))
                    sample_info[temp[0]] = temp[2]
                elif len(temp) == 2:
                    sample_info[temp[0]] = temp[1]
                else:
                    sample_info[line.split("/")[-1].split(".bam")[0]] =line.strip()
                    #self.set_error("%sbam_list格式不正确！", variables=(self.option("bam_list")), code="24500418")
        self.logger.info("样本与bam文件对应关系：{}".format(sample_info))
        return sample_info

    def get_vcf_list(self):
        gevent.sleep(2)
        path = os.path.join(self.output_dir, 'snp_indel')
        files = os.listdir(path)
        vcf = os.path.join(self.output_dir, "vcf.list")
        with open(vcf, 'w') as w:
            for m in files:
                if self.option("call_type") == "samtools":
                    if re.match(r".*\.bed\.vcf\.gz$", m):
                        w.write("{}\n".format(os.path.join(path, m)))
                elif self.option("call_type") == "gatk":
                    if re.match(r".*\.noid\.vcf$", m):
                        w.write("{}\n".format(os.path.join(path, m)))
                elif self.option("call_type") == "freebayes":
                    if re.match(r".*\.bed\.vcf$", m):
                        w.write("{}\n".format(os.path.join(path, m)))
                else:
                    pass

    def bcftools_vcf_run(self):
        self.get_vcf_list()
        self.bcftools_vcf.set_options({
            "vcf_list": os.path.join(self.output_dir, "vcf.list")
        })
        self.bcftools_vcf.on('end', self.set_output, 'vcf_call')
        self.bcftools_vcf.on("end", self.end)
        self.bcftools_vcf.run()

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
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
                # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'snp_indel':
            self.linkdir(obj.output_dir, 'snp_indel')
        elif event['data'] == 'haplotype':
            self.linkdir(obj.output_dir, 'haplotype')
        elif event['data'] == 'gvcf_typing':
            self.linkdir(obj.output_dir, 'gvcf_typing')
        elif event['data'] == "vcf_call":
            self.linkdir(obj.output_dir, 'vcf_call')
        elif event['data'] == "predeal":
            self.linkdir(obj.output_dir, 'predeal')
            self.end()
        else:
            pass

    def run(self):
        super(CallSnpIndelSentieonModule, self).run()
        self.logger.info("开始切割")
        self.genome_split()
        self.logger.info("切割完成")
        if self.option("call_type") == "samtools":
            self.samtools_call_run()
        elif self.option("call_type") == "gatk":
            self.haplotype_run()
        elif self.option("call_type") == "freebayes":
            self.freebayes_call_run()
        elif self.option("call_type") == "sentieon":
            self.haplotype_v2_run()

    def set_bam_list(self):
        """
        samtools和freebayes与gatk的bamlist格式不一致，这里进行处理下
        """
        self.logger.info("开始设置bamlist！")
        with open(self.option("bam_list"), "r") as r, open("{}/bam.list".format(self.work_dir), "w") as w:
            data = r.readlines()
            for line in data:
                temp = line.strip().split("\t")
                w.write(temp[1] + "\n")
        self.logger.info("设置bamlist成功！")

    def end(self):
        super(CallSnpIndelSentieonModule, self).end()



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
            "id": "sentieonsmallnewtest" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "ref_rna_v3.call_snp_indel_sentieon",
            "instant": False,
            "options": dict(
                bam_list="/mnt/ilustre/users/isanger/workspace/20210315/Refrna_n34u_49qekvdhgok1ed33c1m5hf/RnaseqMapping/output/bamlist",
                call_type="sentieon",
                ref_fasta="/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/plants/Brassica_napus/genoscope/dna/Brassica_napus_v4.1.chromosomes.fa",
                des="/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/plants/Brassica_napus/genoscope/biomart/Brassica_napus.biomart",
                des_type="type3",
                ref_gtf="/mnt/ilustre/users/isanger/workspace/20210315/Refrna_n34u_49qekvdhgok1ed33c1m5hf/FileCheck/Brassica_napus.gtf",
                analysis_format="bam",
                align_method="hisat"


                #test1
                # bam_list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test1/Refrna_i-sanger_233268/output/bamlist",
                # call_type="sentieon",
                # ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Magnaporthe_oryzae/MG8_ensembl_45/dna/Magnaporthe_oryzae.MG8.dna.toplevel.fa",
                # des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Magnaporthe_oryzae/MG8_ensembl_45/biomart/biomart1.txt",
                # des_type="type1",
                # ref_gtf="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test1/Refrna_i-sanger_233268/output/ref_and_new.gtf",
                # # analysis_format="cram"

                #test2

                # bam_list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test2/Refrna_majorbio_233735/output/bamlist",
                # call_type="sentieon",
                # ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/dna/Homo_sapiens.GRCh38.dna.toplevel.fa",
                # des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/biomart/biomart1.txt",
                # des_type="type1",
                # ref_gtf="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test2/Refrna_majorbio_233735/output/ref_and_new.gtf",
                # # analysis_format="cram"

                #test3
                # bam_list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test3/majorbio_233634/output/bamlist",
                # call_type="sentieon",
                # ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Castanea_mollissima/v1.0_v1.0/dna/HBY_2_genome_assembly.fasta",
                # des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Castanea_mollissima/v1.0_v1.0/biomart/biomart",
                # des_type="type2",
                # ref_gtf="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test3/majorbio_233634/output/ref_and_new.gtf",
                # # analysis_format="cram"

                # # #test4
                # bam_list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test4/i-sanger_233520/output/bamlist",
                # call_type="sentieon",
                # ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/dna/Homo_sapiens.GRCh38.dna.toplevel.fa",
                # des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/biomart/biomart1.txt",
                # des_type="type1",
                # ref_gtf="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test4/i-sanger_233520/output/ref_and_new.gtf",
                # # analysis_format="cram"

                # # # #test5
                # bam_list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test5/Refrna_majorbio_232290/output/bamlist",
                # call_type="sentieon",
                # ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Cordyceps_militaris/GCF_000225605.1_CmilitarisCM01_v01/dna/GCF_000225605.1_CmilitarisCM01_v01_genomic.fna",
                # des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Cordyceps_militaris/GCF_000225605.1_CmilitarisCM01_v01/biomart/biomart",
                # des_type="type2",
                # ref_gtf="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test5/Refrna_majorbio_232290/output/ref_and_new.gtf",
                # # analysis_format="cram"

                # # test5
                # bam_list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test5/Refrna_majorbio_232290/output/bamlist",
                # call_type="sentieon",
                # ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Cordyceps_militaris/GCF_000225605.1_CmilitarisCM01_v01/dna/GCF_000225605.1_CmilitarisCM01_v01_genomic.fna",
                # des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Cordyceps_militaris/GCF_000225605.1_CmilitarisCM01_v01/biomart/biomart",
                # des_type="type2",
                # ref_gtf="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test5/Refrna_majorbio_232290/output/ref_and_new.gtf",
                # analysis_format="cram"

                # # # test6 majorbio_231937
                # bam_list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test6/Refrna_majorbio_231937/output/bamlist",
                # call_type="sentieon",
                # ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Litsea_cubeba/Lcu_v1.0_Lcu_v1.0/dna/Lcu_new_hic.fasta",
                # des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Litsea_cubeba/Lcu_v1.0_Lcu_v1.0/biomart/biomart",
                # des_type="type2",
                # ref_gtf="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test6/Refrna_majorbio_231937/output/ref_and_new.gtf",
                # analysis_format="cram"

                # # test7 majorbio_230798
                # bam_list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test7/Refrna_majorbio_230798/output/bamlist",
                # call_type="sentieon",
                # ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Triticum_turgidum/Svevo.v1_ensembl_45/dna/Triticum_turgidum.Svevo.v1.dna_sm.toplevel.fa",
                # des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Triticum_turgidum/Svevo.v1_ensembl_45/biomart/biomart3.txt",
                # des_type="type3",
                # ref_gtf="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test7/Refrna_majorbio_230798/output/ref_and_new.gtf",
                # analysis_format="cram"

                # test8 majorbio_229067
                # bam_list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test8/Refrna_majorbio_229067/output/bamlist",
                # call_type="sentieon",
                # ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Triticum_aestivum/iwgsc_refseq/dna/iwgsc_refseqv1.0.genome.fasta",
                # des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Triticum_aestivum/iwgsc_refseq/biomart/iwgsc_refseqv1.0.biomart",
                # des_type="type3",
                # ref_gtf="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test8/Refrna_majorbio_229067/output/ref_and_new.gtf",
                # analysis_format="cram"

            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last_modify:20180404

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import unittest
from biocluster.config import Config
import gevent.subprocess as subprocess



class SentieontestModule(Module):
    """
    用于对bam文件进行call突变位点，该模块中包含了三种方式的call 突变位点，分别是gatks，samtools，freebayes，其中进行call的
    时候是要将所有的bam文件合并在一起，然后将合并后的基因组分成25份进行并行计算，计算完成后再合并到一起
    last modified by hongdong @ 20190215 add sentieon method
    """
    def __init__(self, work_id):
        super(SentieontestModule, self).__init__(work_id)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "ref_rna_v2.common"},
            #{"name": "bam_list", "type": "string"},  # 如果是gatk 样本+建库类型+bam文件，否则样本bam列表
            {"name": "ref_dict", "type": "string"},  # 参考基因组的配置文件，里面可以解析出有多少条染色体
            {"name": "in_bam", "type": "infile", "format": "align.bwa.bam_dir"},  # bam格式文件
            {"name": "call_type", "type": "string", "default": "samtools"},  # call snp的方式
            {"name": "des", "type": "string"},
            {"name": "des_type", "type": "string"},
            {"name": "ref_gtf", "type": "infile", "format": "ref_rna_v2.common"}
        ]
        self.add_option(options)
        self.file_list = []
        self.samtools_call_tools = []
        self.freebayes_call_tools = []
        self.haplotype_tools = []
        self.gvcf_typing_tools = []
        self.samples = []
        self.picards = []
        self.ref_fasta = ""
        self.ref_index1 = ""
        self.ref_index2 = ""
        self.annovars = []
        self.samtools_path = Config().SOFTWARE_DIR + "/miniconda2/bin/"
        self.picard_path = Config().SOFTWARE_DIR + "/bioinfo/gene-structure/"
        self.java_path = Config().SOFTWARE_DIR + "/program/sun_jdk1.8.0/bin/java"
        self.bcftools_vcf = self.add_tool("wgs.bcftool_vcf")
        self.gvcf_typing_v2 = None

    def picard_run(self):
        for i in os.listdir(self.option("in_bam").prop["path"]):
            self.samples.append(i[:-4])
            f_path = os.path.join(self.option("in_bam").prop["path"], i)
            self.logger.info(f_path)  # 打印出f_path的信息，是上一步输出文件的路径
            picard = self.add_tool('denovo_rna_v2.samtools_presnp')
            picard.set_options({
                "in_bam": f_path
            })
            self.picards.append(picard)
            #picard.on("end", self.gatk_run, i[:-4])  # event["data"]传为样本名称
        for j in range(len(self.picards)):
            self.picards[j].on("end", self.set_output, 'picards')
        if len(self.picards) == 1:
            self.picards[0].on("end", self.haplotype_run2)
        else:
            self.on_rely(self.picards, self.haplotype_run2)
        for picard in self.picards:
            picard.run()

    def check_options(self):
        if not self.option("ref_fasta").is_set:
            raise OptionError("缺少ref_fasta参数", code="24500401")
        # if not self.option("bam_list"):
        #     raise OptionError("缺少bam_list参数", code="24500402")
        #if not self.option("ref_dict"):
            #raise OptionError("缺少ref_dict参数", code="24500403")
        if self.option("call_type"):
            if self.option("call_type") not in ["samtools", "gatk", "freebayes", "sentieon"]:
                raise OptionError("call_type方式不合法，必须是samtools, gatk, freebayes, sentieon", code="24500404")
        else:
            raise OptionError("缺少call_type参数", code="24500405")
        ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        ref_file_name = os.path.basename(self.option("ref_fasta").prop['path']).split(".")[0]
        #if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
                #or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
            #raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="24500406")
        #return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def picard_run(self):
        for i in os.listdir(self.option("in_bam").prop["path"]):
            self.samples.append(i[:-4])
            f_path = os.path.join(self.option("in_bam").prop["path"], i)
            self.logger.info(f_path)  # 打印出f_path的信息，是上一步输出文件的路径
            picard = self.add_tool('ref_rna_v2.picard_rna')
            picard.set_options({
                "in_bam": f_path
            })
            self.picards.append(picard)
            #picard.on("end", self.gatk_run, i[:-4])  # event["data"]传为样本名称
        for j in range(len(self.picards)):
            self.picards[j].on("end", self.set_output, 'picards')
        if len(self.picards) == 1:
            self.picards[0].on("end", self.haplotype_run2)
        else:
            self.on_rely(self.picards, self.haplotype_run2)
        for picard in self.picards:
            picard.run()


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
            cmd = "{}samtools faidx {}".format(self.samtools_path, self.option("ref_fasta").prop["path"])
            self.logger.info("开始进行samtools建索引！")
            subprocess.check_call(cmd, shell=True)
            os.link(self.option("ref_fasta").prop["path"] + '.fai', self.ref_index1)
        else:
            os.link(self.option("ref_fasta").prop["path"] + '.fai', self.ref_index1)
        if not os.path.exists(os.path.dirname(self.option("ref_fasta").prop["path"]) + "/" + ".".join(os.path.basename(self.option("ref_fasta").prop["path"]).split(".")[:-1]) + ".dict"):
            dirname = os.path.split(self.option("ref_fasta").prop["path"])[-1]
            dictname = os.path.splitext(dirname)[0] + ".dict"
            self.dict(self.option("ref_fasta").prop["path"], os.path.join(os.path.split(self.option("ref_fasta").prop["path"])[0],dictname))
            os.link(os.path.join(os.path.split(self.option("ref_fasta").prop["path"])[0], dictname), self.ref_index2)
        else:
            os.link(os.path.dirname(self.option("ref_fasta").prop["path"]) + "/" + ".".join(os.path.basename(self.option("ref_fasta").prop["path"]).split(".")[:-1]) + ".dict", self.ref_index2)
        n = 0
        split_file = "{}/split_file".format(self.work_dir)
        if os.path.exists(split_file) or os.path.isdir(split_file):
            os.system('rm -r %s' % split_file)
            os.mkdir(split_file)
        else:
            os.mkdir(split_file)
        with open(self.ref_index2, "r") as r:
            data = r.readlines()[1:]
            for line in data:
                temp = line.strip().split("\t")
                chr_name = temp[1].strip().split(":")[1]
                length = temp[2].strip().split(':')[1]
                n += 1
                file_name = n % 25
                if self.option("call_type") == "gatk":
                    with open(split_file + "/{}.intervals".format(file_name), 'a+') as w:
                        w.write(chr_name + "\n")
                else:
                    with open(split_file + "/{}.bed".format(file_name), 'a+') as w:
                        w.write(chr_name + "\t" + "0\t" + length + "\n")
        self.get_bed_intervals_file(split_file)

    def get_bed_intervals_file(self, split_file):
        """
        获取分割后的bed文件或者intervals文件的个数
        :return:
        """
        files = os.listdir(split_file)
        if len(files) == 0:
            self.set_error("基因组分割后的列表为空，不能进行后续分析下，请检查ref.dict文件！", code="24500413")
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
            self.set_error("samtools_call_tools列表为空！", code="24500414")

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
            self.set_error("freebayes_call_tools列表为空！", code="24500415")

    def haplotype_run2(self):
        #sample_info = self.get_sample_bam()
        n = 0
        with open(os.path.join(self.work_dir,"bam_list"),"w") as blst:
          for i in os.listdir(os.path.join(self.output_dir,"picards")):
            if i.endswith("bam"):
                bam_path=os.path.join(os.path.join(self.output_dir,"picards"),i)
                blst.write(i.split(".")[0]+"\t"+bam_path+"\n")
        self.bam_list=os.path.join(self.work_dir,"bam_list")
        with open(self.bam_list) as f:
            lines = f.readlines()
            for line in lines:
                sample_name = line.strip().split("\t")[0]
                sample_path = line.strip().split("\t")[1]
                haplotype = self.add_tool("ref_rna_v2.sentieontest")
                options = ({
                    #"ref_genome":"customer_mode",
                    "bam_file": sample_path,
                    "fa_file": self.option("ref_fasta").prop["path"],
                    "name": sample_name
                })
                haplotype.set_options(options)
                self.haplotype_tools.append(haplotype)
            n += 1
            for j in range(len(self.haplotype_tools)):
                self.haplotype_tools[j].on('end', self.set_output, 'haplotype')
            if self.haplotype_tools:
                if len(self.haplotype_tools) > 1:
                    self.on_rely(self.haplotype_tools, self.gvcf_typing_v2_run)
                elif len(self.haplotype_tools) == 1:
                    self.haplotype_tools[0].on('end', self.gvcf_typing_v2_run)
                for tool in self.haplotype_tools:
                    gevent.sleep(1)
                    tool.run()
            else:
                self.set_error("haplotype_tools列表为空！", code="24500416")

    def gvcf_typing_run(self):
        n = 0
        for m in self.file_list:
            gvcf_typing = self.add_tool("ref_rna_v2.gvcf_typing")
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
            self.set_error("gvcf_typing_tools列表为空！", code="24500417")

    def haplotype_v2_run(self):
        sample_info = self.get_sample_bam()
        n = 0
        for keys in sample_info.keys():
            haplotype = self.add_tool("wgs_v2.haplotyper_v2")
            haplotype.set_options({
                "fa_file": self.option("ref_fasta").prop["path"],
                "bam_file": sample_info[keys],
                "name": keys + '.g.vcf'
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
            raise Exception("haplotype_tools列表为空！")
        for tool in self.haplotype_tools:
            gevent.sleep(1)
            tool.run()

    def gvcf_typing_v2_run(self):
        """
        DE1_10  /mnt/ilustre/users/caixia.tian/newmdt/Project/MJ20180723036-gaojun/03.ref/08.haplotype/DE1_10.g.vcf
        :return:
        """
        self.make_gvcf_list()
        self.gvcf_typing_v2 = self.add_tool("wgs_v2.gvcf_typing_v2")
        self.gvcf_typing_v2.set_options({
            "vcf_list": self.output_dir + "/haplotype/gvcf.list",
            "fa_file": self.option("ref_fasta").prop["path"]
        })
        self.gvcf_typing_v2.on('end', self.set_output, 'vcf_call')
        self.gvcf_typing_v2.on('end', self.vcf_filter_run)
        self.gvcf_typing_v2.run()

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
                    sample_name=line.split("/")[-1].split(".bam")[0]
                    sample_info[sample_name]=line.strip()
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
        self.bcftools_vcf.on("end", self.vcf_filter_run)
        self.bcftools_vcf.run()

    def vcf_filter_run(self):
        self.vcffilter = self.add_tool('ref_rna_v2.vcf_filter_gatk')
        options = {
            "ref_genome": "customer_mode",
            "input_file": self.gvcf_typing_v2.output_dir + "/pop.variant.vcf",
            "combine_vcf": False,
            "ref_fasta": self.ref_fasta,
        }
        self.vcffilter.set_options(options)
        #self.vcffilter.on("end", self.snp_anno)
        self.vcffilter.on("end", self.set_output, 'predeal')
        self.logger.info("vcf_filter is running!")
        self.vcffilter.run()

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
        elif event['data'] == "picards":
            self.linkdir(obj.output_dir, 'picards')
        elif event['data'] == "predeal":
            self.linkdir(obj.output_dir, 'predeal')
            self.end()
        else:
            pass

    def run(self):
        super(SentieontestModule, self).run()
        self.logger.info("开始切割")
        self.genome_split()
        self.logger.info("切割完成")
        if self.option("call_type") == "samtools":
            self.samtools_call_run()
        elif self.option("call_type") == "gatk":
            self.picard_run()
        elif self.option("call_type") == "freebayes":
            self.freebayes_call_run()
        elif self.option("call_type") == "sentieon":
            self.picard_run()

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
        super(SentieontestModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "RnaseqMapping_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "denovo_rna.sentieontest",
            "instant": False,
            "options": dict(
                ref_fasta="/mnt/ilustre/users/sanger-dev/workspace/20190416/Denovorna_tsg_33857/DenovoAssemble2/output/Trinity.filter.unigene.fasta",
                #bam_list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/gatk_new_test/bamlist_new",
                #ref_dict="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38_Ensembl_96/dna/Mus_musculus.GRCm38.dna.toplevel.old.dict",
                call_type="sentieon",
                #in_bam="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/gatk_diff_methods/bam_realign",
                in_bam="/mnt/ilustre/users/sanger-dev/workspace/20190416/Denovorna_tsg_33857/Snp/bam_folder",
                #des = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38_Ensembl_96/biomart/biomart.txt",
                #des_type = "type1",
                #ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38_Ensembl_96/gtf/Mus_musculus.GRCm38.96.gtf",

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
        unittest.main()
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


class CallSnpIndelGatkModule(Module):
    """
    用于对bam文件进行call突变位点，该模块中包含了三种方式的call 突变位点，分别是gatks，samtools，freebayes，其中进行call的
    时候是要将所有的bam文件合并在一起，然后将合并后的基因组分成25份进行并行计算，计算完成后再合并到一起
    last modified by hongdong @ 20190215 add sentieon method
    """
    def __init__(self, work_id):
        super(CallSnpIndelGatkModule, self).__init__(work_id)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "denovo_rna_v2.trinity_fasta"},
            {"name": "fq_type", "type": "string", "default": "PE"},  # fq类型，PE、SE
            {"name": "fq_list", "type": "string"},  # 文件里面根据单双端的不同有3列或者2列
            #{"name": "bam_list", "type": "string"},  # 如果是gatk 样本+建库类型+bam文件，否则样本bam列表
            #{"name": "ref_dict", "type": "string"},  # 参考基因组的配置文件，里面可以解析出有多少条染色体
            {"name": "in_bam", "type": "string", 'default':None},
            {"name": "call_type", "type": "string", "default": "gatk"},  # call snp的方式
            {"name": "anno", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "cds_bed", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "allt2g", "type": "infile", "format": "denovo_rna_v2.common"},
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
        self.bwa_tools = []
        self.samtools_path = Config().SOFTWARE_DIR + "/miniconda2/bin/"
        self.picard_path = Config().SOFTWARE_DIR + "/bioinfo/gene-structure/"
        self.java_path = Config().SOFTWARE_DIR + "/program/sun_jdk1.8.0/bin/java"
        self.bcftools_vcf = self.add_tool("wgs.bcftool_vcf")
        self.gvcf_typing_v2 = None
        if self.option('in_bam'):
            self.in_bam=self.option('in_bam')

    def check_options(self):
        if not self.option("ref_fasta").is_set:
            raise OptionError("缺少ref_fasta参数", code="200101")
        # if not self.option("bam_list"):
        #     raise OptionError("缺少bam_list参数", code="24500402")
        # if not self.option("ref_dict"):
        #     raise OptionError("缺少ref_dict参数", code="24500403")
        if self.option("call_type"):
            if self.option("call_type") not in ["samtools", "gatk", "freebayes", "sentieon"]:
                raise OptionError("call_type方式不合法，必须是samtools, gatk, freebayes, sentieon", code="200102")
        else:
            raise OptionError("缺少call_type参数", code="200103")
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

    def bwa_run(self, fq_type, fq_list):
        self.in_bam = os.path.join(self.work_dir,"bwa")
        n = 0
        with open(fq_list) as f:
            for line in f:
                self.bwa = self.add_tool("denovo_rna_v2.bwa")
                n += 1
                self.bwa.set_options({
                    "fq_list": line.strip(),
                    "trinity_fa": self.option('ref_fasta'),
                    "fq_type": fq_type,
                    "num": str(n)
                })
                self.bwa_tools.append(self.bwa)
        for j in range(len(self.bwa_tools)):
            self.bwa_tools[j].on('end', self.set_output, 'bwa')
        if self.bwa_tools:
            if len(self.bwa_tools) > 1:
                self.on_rely(self.bwa_tools, self.picard_run)
            elif len(self.bwa_tools) == 1:
                self.bwa_tools[0].on('end', self.picard_run)
        else:
            self.set_error("self.bwa_tools列表为空！", code = "200101")
        for tool in self.bwa_tools:
            gevent.sleep(1)
            tool.run()

    def picard_run(self):
        for i in os.listdir(self.in_bam):
            self.samples.append(i[:-4])
            f_path = os.path.join(self.in_bam, i)
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
            self.picards[0].on("end", self.haplotype_run)
        else:
            self.on_rely(self.picards, self.haplotype_run)
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

    def get_bed_intervals_file(self, split_file):
        """
        获取分割后的bed文件或者intervals文件的个数
        :return:
        """
        files = os.listdir(split_file)
        if len(files) == 0:
            self.set_error("基因组分割后的列表为空，不能进行后续分析下，请检查ref.dict文件！", code="200102")
        else:
            for m in files:
                self.file_list.append(os.path.join(split_file, m))
            self.logger.info("file_list为{}".format(self.file_list))

    def haplotype_run(self):
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
                haplotype = self.add_tool("ref_rna_v3.gatk")
                options = ({
                    "ref_genome":"customer_mode",
                    "sample_bam": sample_path,
                    "ref_fa": self.option("ref_fasta").prop["path"],
                    "sample_id": sample_name
                })
                haplotype.set_options(options)
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
                self.set_error("haplotype_tools列表为空！", code="200103")

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
            self.set_error("gvcf_typing_tools列表为空！", code="200104")

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
            "input_file": self.bcftools_vcf.output_dir + "/pop.variant.vcf",
            "combine_vcf": False,
            "ref_fasta": self.ref_fasta,
        }
        self.vcffilter.set_options(options)
        self.vcffilter.on("end", self.set_output, "vcf_filter")
        if self.inter_inbam:
            self.vcffilter.on("end", self.snpfinal_run)
        else:
            self.vcffilter.on("end",self.end)
        self.logger.info("vcf_filter is running!")
        self.vcffilter.run()

    def snpfinal_run(self):
        with open(os.path.join(self.work_dir,"bam_list")) as len_bam:
            sample_num = len(len_bam.readlines())
        self.snpfinal = self.add_tool('denovo_rna_v3.snpfinal_new2')
        self.snpfinal.set_options({
            "bamlist": sample_num,
            'method':'gatk',
            "filted_snp_vcf": self.vcffilter.output_dir + "/pop.snp.filter.recode.vcf",
            "filted_indel_vcf": self.vcffilter.output_dir + "/pop.indel.filter.recode.vcf",
            'cds_bed': self.option("cds_bed"),
            'allt2g' : self.option('allt2g'),
            'anno' : self.option('anno')
            #"qual": self.option('qual'),
           # "dp": self.option('dp'),
        })
        # 这个地方需要在最后跑的那个tool那里加一个on('end', self.end)，这样module才会终止，并且这个
        # self.end和这最后一个tool的set_output是一起跑的，不会先end导致set_output失败
        self.snpfinal.on('end', self.set_output, 'snpfinal')
        self.snpfinal.run()


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

    def linkdirs(self, dirpath, dirname):
        allfiles = os.listdir(os.path.join(dirpath, "ready_sort_bam"))
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, "ready_sort_bam", i) for i in allfiles]
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
        elif event['data'] == "bwa":
            self.linkdirs(obj.output_dir, 'bwa')
        elif event['data'] == "snpfinal":
            self.linkdir(obj.output_dir, 'snpfinal')
            self.end()
        else:
            pass

    def run(self):
        super(CallSnpIndelGatkModule, self).run()
        self.logger.info("开始切割")
        self.genome_split()
        self.logger.info("切割完成")
        if self.option("in_bam"):
            self.in_bam=self.option("in_bam")
            self.inter_inbam = self.option('in_bam')
        else:
            self.in_bam=None
            self.inter_inbam = None
        if self.option("call_type") == "gatk":
            if self.in_bam:
                self.picard_run()
            else:
                self.bwa_run(self.option("fq_type"), self.option("fq_list"))
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
        super(CallSnpIndelGatkModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Snptest_" + str(random.randint(1, 10000))+"yyyy",
            "type": "module",
            "name": "denovo_rna_v3.call_snp_indel_gatk",
            "instant": False,
            "options": dict(
                ref_fasta="/mnt/ilustre/users/sanger-dev/workspace/20190416/Denovorna_tsg_33857/DenovoAssemble2/output/Trinity.filter.unigene.fasta",
                fq_list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/denovo_snp/final_snp/test_files/test_bam/fq/fq_list",
                fq_type="PE",
                cds_bed="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/denovo_snp/final_snp/cds/all_predicted.bed",
                allt2g="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/denovo_snp/final_snp/cds/all_tran2gen.txt",
                anno="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/denovo_snp/final_snp/Snpfinal_new/gene_anno_detail.xls",
                in_bam="/mnt/ilustre/users/sanger-dev/workspace/20190814/Single_denovo_snp8413yyyy/Snp/bwa"

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
        unittest.main()
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


class CallSomaticSentieonModule(Module):
    """
    用于对bam文件进行call突变位点，该模块中包含了三种方式的call 突变位点，分别是gatks，samtools，freebayes，其中进行call的
    时候是要将所有的bam文件合并在一起，然后将合并后的基因组分成25份进行并行计算，计算完成后再合并到一起
    last modified by fwy @ 20190620 改写为ref_rna_v2专用
    """
    def __init__(self, work_id):
        super(CallSomaticSentieonModule, self).__init__(work_id)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "bam_list", "type": "string"},  # 如果是gatk 样本+建库类型+bam文件，否则样本bam列表
            {"name": "tn_pair", "type": "infile", "format": "ref_rna_v2.common"},
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
        self.tnhaplotypers =[]
        self.gvcf_typing_tools = []
        #self.bcftools_vcf = self.add_tool("wgs.bcftool_vcf")
        self.gvcf_typing_v2 = None
        self.ref_fasta=""
        self.ref_index1=""
        self.ref_index2=""
        self.samtools_path = Config().SOFTWARE_DIR + "/bioinfo/ref_rna_v2/miniconda2/bin/"
        self.picard_path = Config().SOFTWARE_DIR + "/bioinfo/gene-structure/"
        self.java_path=Config().SOFTWARE_DIR+"/program/sun_jdk1.8.0/bin/java"
        self.vcffilter = self.add_tool('medical_transcriptome.snp.vcf_filter_gatk')
        self.predeal = self.add_tool('medical_transcriptome.snp.predeal_snpresults')
        self.somatic_predict = self.add_module('medical_transcriptome.somatic.somatic_predict_new')
        self.annovar = self.add_tool('medical_transcriptome.somatic.annovar')
        self.annovars = []
        self.pair_relation = ""

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

    def sentieon_somatic_run(self):
        sample_info = self.get_sample_bam()
        n = 0
        for keys in sample_info.keys():
            haplotype = self.add_module("medical_transcriptome.somatic.pre4sentieon")
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
            self.haplotype_tools[j].on('end', self.set_output, 'somatic_pre')
        if self.haplotype_tools:
            if len(self.haplotype_tools) > 1:
                self.on_rely(self.haplotype_tools, self.run_TNhaplotyper)
            elif len(self.haplotype_tools) == 1:
                self.haplotype_tools[0].on('end', self.run_TNhaplotyper)
        else:
            self.set_error("TNhaplotyper_pre_tools列表为空！")
        for tool in self.haplotype_tools:
            gevent.sleep(1)
            tool.run()



    def run_TNhaplotyper(self):
        self.pair_relation = self.get_pair_relation()
        for pair in self.pair_relation:
            normal,tumor = pair
            normal_path = os.path.join(self.output_dir,"somatic_pre","{}.bam".format(normal))
            tumor_path = os.path.join(self.output_dir, "somatic_pre", "{}.bam".format(tumor))
            tnhaplotyper = self.add_tool("medical_transcriptome.somatic.sentieon_somatic")
            tnhaplotyper.set_options({
                "tumor_file": tumor_path,
                "tumor_name": tumor,
                "normal_file": normal_path,
                "normal_name":normal,
                "fa_file": self.option("ref_fasta").prop["path"]
            })
            self.tnhaplotypers.append(tnhaplotyper)
        for j in range(len(self.tnhaplotypers)):
            self.tnhaplotypers[j].on('end', self.set_output, 'tnhaplotyper')
        if self.tnhaplotypers:
            if len(self.tnhaplotypers) > 1:
                self.on_rely(self.tnhaplotypers, self.snp_anno)
            elif len(self.tnhaplotypers) == 1:
                self.tnhaplotypers[0].on('end', self.snp_anno)
        else:
            self.set_error("TNhaplotyper_tools列表为空！")
        for tool in self.tnhaplotypers:
            gevent.sleep(1)
            tool.run()


    def get_pair_relation(self):
        pair_relation = []
        with open(self.option("tn_pair").prop["path"],"r") as f:
            f.readline()
            for line in f.readlines():
                line = line.strip().split("\t")
                normal = line[0]
                tumor = line[1]
                pair_relation.append((normal,tumor))
        return pair_relation


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
        options = {
            "ref_genome": "customer_mode",
            "input_file": self.gvcf_typing_v2.output_dir + "/pop.variant.vcf",
            "combine_vcf": False,
            "ref_fasta": self.option("ref_fasta").prop["path"],
        }
        self.vcffilter.set_options(options)
        self.vcffilter.on("end", self.snp_anno)
        self.logger.info("vcf_filter is running!")
        self.vcffilter.run()


    def snp_anno(self):
        # self.annovar = self.add_tool('medical_transcriptome.somatic.annovar')
        options = {
            "ref_genome": "customer_mode",
            "ref_fasta":self.option("ref_fasta").prop["path"],
            "input_file": os.path.join(self.output_dir ,"tnhaplotyper"),
            "combine_vcf": True,
            "des_type": self.option("des_type"),
            "des": self.option("des"),
            "ref_gtf":self.option("ref_gtf"),
        }
        self.annovar.set_options(options)
        # self.annovars.append(self.annovar)
        self.logger.info("set output done")
        self.annovar.on("end", self.run_predeal, "annovar")
        self.annovar.on("end", self.predict)
        # annovar.on('end', self.end)
        self.annovar.run()

    def predict(self):
        options = {
            "ref_fasta":self.option("ref_fasta").prop["path"],
            "combine_vcf": os.path.join(self.annovar.work_dir ,"filter_Combine_Variants"),
        }
        self.somatic_predict.set_options(options)
        self.logger.info("set output done")
        # self.somatic_predict.on("end", self.predeal, "annovar")
        # annovar.on('end', self.end)
        self.somatic_predict.run()

    def run_predeal(self, event):
        obj = event["bind_object"]
        options = {
            "method": "gatk",
            "snp_detail": obj.output_dir + "/snp_anno.xls",
            "snp_anno": obj.work_dir + "/snp_annotation.xls"
        }
        self.predeal.set_options(options)
        self.logger.info("predeal done")
        self.predeal.on("end", self.set_output, "predeal")
        self.predeal.run()


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
        elif event['data'] == 'somatic_pre':
            self.linkdir(obj.output_dir, 'somatic_pre')
        elif event['data'] == 'tnhaplotyper':
            self.linkdir(obj.output_dir, 'tnhaplotyper')
        elif event['data'] == "vcf_call":
            self.linkdir(obj.output_dir, 'vcf_call')
        elif event['data'] == "predeal":
            self.linkdir(obj.output_dir, 'predeal')
        elif event['data'] == "somatic_predict":
            self.linkdir(obj.output_dir, 'somatic_predict')
        elif event['data'] == "end":
            self.end()
        else:
            pass

    def set_output_final(self):
        self.linkdir(self.predeal.output_dir, 'predeal')
        self.linkdir(self.somatic_predict.output_dir, 'somatic_predict')
        self.end()


    def run(self):
        super(CallSomaticSentieonModule, self).run()
        self.logger.info("开始切割")
        # self.genome_split()
        self.logger.info("切割完成")
        self.on_rely([self.somatic_predict, self.predeal], self.set_output_final)
        self.sentieon_somatic_run()

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
        super(CallSomaticSentieonModule, self).end()



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
            "name": "medical_transcriptome.somatic.call_somatic_sentieon",
            "instant": False,
            "options": dict(
                bam_list="/mnt/ilustre/users/sanger-dev/workspace/20201117/Snp_tsg_38314_8828_9278/bamlist_new",
                call_type="sentieon",
                ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa",
                des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/biomart/biomart.txt",
                des_type="type1",
                ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/gtf/Homo_sapiens.GRCh38.96.gtf",
                analysis_format="bam",
                align_method="hisat",
                tn_pair = "/mnt/ilustre/users/sanger-dev/workspace/20201117/Snp_tsg_38314_8828_9278/tn",

                # bam_list="/mnt/ilustre/users/sanger-dev/workspace/20201105/MedicalTranscriptome_v6dpivfr84k4ooq2lmvh47dpfs/RnaseqMapping/output/bamlist",
                # call_type="sentieon",
                # ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa",
                # des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/biomart/biomart.txt",
                # des_type="type1",
                # ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/gtf/Homo_sapiens.GRCh38.96.gtf",
                # analysis_format="bam",
                # align_method="hisat",
                # tn_pair="/mnt/ilustre/users/sanger-dev/workspace/20201105/MedicalTranscriptome_v6dpivfr84k4ooq2lmvh47dpfs/RnaseqMapping/output/tn",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
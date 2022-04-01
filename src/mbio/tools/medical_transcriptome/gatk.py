# -*- coding: utf-8 -*-
# __author__ = 'chenyanyan'
# last modify 2016.09.28

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import glob
from biocluster.core.exceptions import OptionError
import subprocess
import re 
import shutil
import json


class GatkAgent(Agent):
    """
    SNP工具,gatk的每一步里面的R都是fa文件和它的字典一起的，如若不是，gatk将无法正确运行
    """
    def __init__(self, parent):
        super(GatkAgent, self).__init__(parent)
        options = [
            # {"name":"ref_genome_custom", "type": "infile", "format": "sequence.fasta"},
            {"name": "ref_genome", "type": "string"},  # 本地参考基因组的名字，从已有列表里面选择
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # 参考基因组文件,需要和picard所建立的dict和samtools所建立的fai文件一起
            {"name": "input_bam", "type": "infile", "format": "align.bwa.bam"},  # bam文件类型，输入
            {"name": "is_indexed", "type": "bool", "default": True},
            # {"name": "known_vcf", "type": "string", "default": ""},
        ]
        self.add_option(options)
        self.step.add_steps('gatk')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.gatk.start()
        self.step.update()

    def step_end(self):
        self.step.gatk.finish()
        self.step.update()
        
    def check_options(self):
        #if not self.option("ref_genome") in self._ref_genome_lst:
        #    raise OptionError("请选择参考基因组！")
        if self.option("ref_genome") == "customer_mode" and not self.option("ref_fa").is_set:
            raise OptionError("自定义参考基因组文件未提供！", code = "32200901")
        if not self.option("input_bam").is_set:
            raise OptionError("用于分析的bam文件未提供！", code = "32200902")
            
    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 20
        self._memory = '80G'
        
    def end(self):
      
        super(GatkAgent, self).end()    


class GatkTool(Tool):
    """
    GATK3.6
    """
    def __init__(self, config):
        super(GatkTool, self).__init__(config)
        self.gatk_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/"
        self.gatk4_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/gatk-4.beta.5/"
        self.java_path = "/program/sun_jdk1.8.0/bin/"
        self.samtools_path = "/bioinfo/align/samtools-1.3.1/samtools"
        self.picard_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/"
        self.known_vcf = "/mnt/ilustre/users/sanger-dev/sg-users/qindanhua/denovo_rna/test_file/human_ref/known.vcf"

    def dict(self, ref_fasta, dict_name):
        """
        使用picard对参考基因组构建字典
        """

        cmd = "program/sun_jdk1.8.0/bin/java -jar {}picard.jar CreateSequenceDictionary R={} O={}"\
            .format(self.picard_path, ref_fasta, dict_name)
        if os.path.exists(dict_name):
            os.remove(dict_name)
        print cmd
        self.logger.info("开始用picard对参考基因组构建字典")
        command = self.add_command("dict", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("参考基因组构建dict done!")
        else:
            self.set_error("构建dict过程error！", code = "32200903")


    def samtools_faidx(self, ref_fasta):

        """
        使用samtools对参考基因组建索引
        """
        cmd = "{} faidx {}".format(self.samtools_path, ref_fasta)
        self.logger.info("开始进行samtools建索引！")
        command = self.add_command("samtools_faidx", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("samtools建索引完成！")
        else:
            self.logger.info("重新运行samtools建索引过程")
            if command.return_code == None:
                command.rerun()
                if command.return_code == 0:
                    self.logger.info("samtools建索引完成！")
                else:
                    self.set_error("samtools建索引出错！", code = "32200904")


    # add by qindanhua 20170221
    def realignment(self, ref_fasta):
        cmd = "{}java -jar {}GenomeAnalysisTK.jar -T RealignerTargetCreator -R {} -I {} -o {} -known {}"\
            .format(self.java_path, self.gatk_path, ref_fasta, "split.bam", "realn.intervals", self.known_vcf)
        self.logger.info("开始进行indel比对")
        command = self.add_command("realignment", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("indel比对完成！")
        else:
            self.set_error("indel比对出错！", code = "32200905")


    def indel_realn(self, ref_fasta):
        cmd = "{}java -jar {}GenomeAnalysisTK.jar -T IndelRealigner -R {} -I {} -o {} -targetIntervals {}"\
            .format(self.java_path, self.gatk_path, ref_fasta, "split.bam", "indel_realn.bam", "realn.intervals")
        self.logger.info("开始进行indel_realn")
        command = self.add_command("indel_realn", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("indel_realn完成！")
        else:
            self.set_error("indel_realn出错！", code = "32200906")


    def base_recalibrator(self, ref_fasta):
        cmd = "{}java -jar {}GenomeAnalysisTK.jar -T BaseRecalibrator -R {} -I {} -o {} -knownSites {}"\
            .format(self.java_path, self.gatk_path, ref_fasta, "indel_realn.bam", "bqsr.grp", self.known_vcf)
        self.logger.info("开始进行base_recalibrator")
        command = self.add_command("base_recalibrator", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("base_recalibrator完成！")
        else:
            self.set_error("base_recalibrator出错！", code = "32200907")


    def print_reads(self, ref_fasta):
        cmd = "{}java -jar {}GenomeAnalysisTK.jar -T PrintReads -R {} -I {} -o {}"\
            .format(self.java_path, self.gatk_path, ref_fasta, "indel_realn.bam", "realn.bam")
        self.logger.info("开始进行print_reads")
        command = self.add_command("print_reads", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("print_reads完成！")
        else:
            self.set_error("print_reads出错！", code = "32200908")


    def gatk_split(self, ref_fasta, bam):
        """
        step0: gatk split
        """

        cmd = "{}java -jar {}GenomeAnalysisTK.jar -T SplitNCigarReads " \
              "-R {} -I {} -o {} -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS"\
            .format(self.java_path, self.gatk_path, ref_fasta, bam, "split.bam")

        #cmd = "{}java -Dsnappy.disable=true -jar {}/gatk-package-4.beta.5-local.jar SplitNCigarReads -R {} -I {} -O {}"\
        #    .format(self.java_path, self.gatk_path, ref_fasta, bam, "split.bam")
        self.logger.info("开始进行split NCigar reads")
        command = self.add_command("splitncigarreads", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("split NCigar reads完成！")
        else:
            self.set_error("split NCigar reads出错！", code = "32200909")

        
    def gatk_vc(self, ref_fasta, split_bam):
        """
        step1：gatk variant calling
        """
        cmd = "{}java -jar {}GenomeAnalysisTK.jar -T HaplotypeCaller -R {} -I {} -o {}"\
            .format(self.java_path, self.gatk_path, ref_fasta, split_bam, "output.vcf")
        self.logger.info("开始进行variant calling")
        command = self.add_command("variant calling", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("variant calling完成！")
        else:
            self.set_error("variant calling出错！", code = "32200910")

            
    def gatk_vf(self, ref_fasta, variantfile):
        """
        step2：gatk variant filtering
        """
        cmd = self.config.SOFTWARE_DIR + "/"
        cmd += "{}java -jar {}GenomeAnalysisTK.jar -R {} -T VariantFiltration -V {} -window 35 -cluster 3 " \
               "-filterName FS -filter 'FS > 30.0' -filterName QD -filter 'QD < 2.0' -o {}"\
            .format(self.java_path, self.gatk_path, ref_fasta, variantfile, "filtered.vcf")
        self.logger.info("开始进行variant filtering")
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('variant filtering成功')
        except subprocess.CalledProcessError:
            self.logger.info('variant filtering失败')
            self.set_error('无法生成vcf文件', code = "32200911")

            
    def run(self):
        super(GatkTool, self).run()
        if self.option("ref_genome") == "customer_mode" and self.option("ref_fa").is_set:
            ref = self.option("ref_fa").prop["path"]
            ref_name = os.path.split(ref)[-1]
            if ref_name.endswith(".fasta") or ref_name.endswith(".fa"):
                shutil.copy(ref, self.work_dir)  # 将参考基因组移动到工作目录下
            else:
                ref_name += ".fa"
                shutil.copy(ref, self.work_dir + "/" + ref_name)
            self.logger.info("参考基因组移动至当前工作目录下！")
            
            ref_current = os.path.join(self.work_dir, ref_name)  # 当前工作目录下的参考基因组的路径
            self.logger.info(ref_current)  # 显示当前参考基因组路径
            dict_name = os.path.splitext(ref_name)[0] + ".dict"   # 所建字典的名字，是一个字符串
            self.logger.info("参考基因组为自定义模式的情况下建字典！")
            if not os.path.exists(dict_name):
                self.dict(ref_current, dict_name)  # 参考基因组
            self.logger.info("参考基因组字典建立完成！")
            ref_fai = ref_name + ".fai"
            if not os.path.exists(os.path.join(self.work_dir, ref_fai)):
                self.samtools_faidx(ref_current)
            self.gatk_split(ref_current, self.option("input_bam").prop["path"])
            split_bam = os.path.join(self.work_dir, "split.bam")
            self.gatk_vc(ref_current, split_bam)
            vcf_path = os.path.join(self.work_dir, "output.vcf")
            self.gatk_vf(ref_current, vcf_path)
        else:
            self.logger.info("在参考基因组从数据库中选择时，运行star")
            ref_genome_json = self.config.SOFTWARE_DIR + "/database/refGenome/scripts/ref_genome.json"
            with open(ref_genome_json, "r") as f:
                dict = json.loads(f.read())
                ref = dict[self.option("ref_genome")]["ref_dict"]  # 是ref_dict的路径，作为参数传给比对函数
                ref_current = ref
            self.gatk_split(ref, self.option("input_bam").prop["path"])
            if self.option("ref_genome") == "human":
                self.realignment(ref_current)
                self.indel_realn(ref_current)
                self.base_recalibrator(ref_current)
                self.print_reads(ref_current)
                split_bam = os.path.join(self.work_dir, "realn.bam")
            else:
                split_bam = os.path.join(self.work_dir, "split.bam")
            self.gatk_vc(ref, split_bam)
            vcf_path = os.path.join(self.work_dir, "output.vcf")
            self.gatk_vf(ref, vcf_path)
        outputs = os.listdir(os.getcwd())
        for i in outputs:
            if re.match(r"filtered.vcf*", i):
                shutil.copy(i, self.output_dir)
        self.end()

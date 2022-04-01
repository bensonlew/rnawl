#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'fwy'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import subprocess
import shutil
import json
import glob
import os
import pandas as pd
import collections
from collections import OrderedDict
import re
import datetime
import unittest
import random


class GatkToolsAgent(Agent):
    """
    gatk对vcf文件进行处理
    """

    def __init__(self, parent):
        super(GatkToolsAgent, self).__init__(parent)
        options = [
            {"name": "vcf_list", "type": "infile", "format":  "ref_rna_v2.common"},  # 输入文件
            {"name": "target", "type": "string","default": "merge_vcfs"},  # 处理方式
            {"name": "out_vcf", "type": "string"},  # 输出vcf文件

            {"name": "ref_genome", "type": "string"},  # 参考基因组类型
            {"name": "input_file", "type": "infile", "format": "gene_structure.vcf,gene_structure.vcf_dir"},  # 输入文件
            {"name": "ref_fasta", "type": "infile", "format": "ref_rna_v2.common"},  # 输入文件,参考基因组文件
            {"name": "combine_vcf", "type": "bool", "default": False},  # 输入文件

        ]
        self.add_option(options)
        self.step.add_steps('vcffiltergatk')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.vcffiltergatk.start()
        self.step.update()

    def step_end(self):
        self.step.vcffiltergatk.finish()
        self.step.update()

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option("input_file").is_set:
            raise OptionError("请输入VCF格式文件", code="33709802")
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 11
        self._memory = '100G'

    def end(self):
        super(GatkToolsAgent, self).end()


class GatkToolsTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(GatkToolsTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.gatk_path = self.config.SOFTWARE_DIR + "/bioinfo/ref_rna_v2/gatk-4.1.4.1/gatk-4.1.4.1/"
        # self.gatk_path = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/gatk_4.1.4test/software/gatk-4.1.4.1/"
        self.java_path = "/program/sun_jdk1.8.0/bin/"
        self.picard_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/"

    def dict(self):
        """
        使用picard对参考基因组构建字典
        """
        if 'path' in self.option("ref_fasta").prop.keys():
            ref_fasta = self.option("ref_fasta").prop["path"]
        else:
            ref_fasta = self.option("ref_fasta")
        dict_name = os.path.dirname(ref_fasta) + "/" + ".".join(os.path.basename(ref_fasta).split(".")[:-1]) + ".dict"
        cmd = "program/sun_jdk1.8.0/bin/java -jar {}picard.jar CreateSequenceDictionary R={} O={}"\
            .format(self.picard_path, ref_fasta, dict_name)
        if os.path.exists(dict_name):
            os.remove(dict_name)
        self.logger.info("开始用picard对参考基因组构建字典")
        self.logger.info("cmd")
        command = self.add_command("dict", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("参考基因组构建dict done!")
        else:
            self.set_error("构建dict过程error！", code = "33709817")

    def merge_vcfs(self):
        vcfs =[]
        with open(self.option("vcf_list").prop["path"]) as f:
            for line in f.readlines():
                vcf = line.strip()
                vcfs.append(vcf)
        cmd = "{}java -Xmx50G -jar {}gatk-package-4.1.4.1-local.jar MergeVcfs".format(self.java_path,self.gatk_path)
        for vcf in vcfs:
            cmd += " -I {}".format(vcf)
        cmd += ""

        self.logger.info(cmd)
        self.logger.info("开始进行vcf_merge")
        command1 = self.add_command("vcf_merge", cmd).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("vcf_merge完成！")
        else:
            self.set_error("vcf_merge出错！")
            self.set_error("vcf_merge出错！")

    def combine_vcf(self):
        samples_option = ''
        # 因为samtools产生的就是一个vcf,在module的参数,不能作为glob的输入
        if os.path.isdir(self.option("input_file").prop["path"]):
            vcf_files = glob.glob('{}/*.vcf'.format(self.option("input_file").prop["path"]))
            for vf in sorted(vcf_files):
                samples_option += ' --variant {}'.format(vf)
        else:
            vcf_files = self.option("input_file").prop["path"]
            samples_option += ' --variant {}'.format(vcf_files)

        cmd = '{}java -jar {}gatk-package-4.1.4.1-local.jar -R {} -T CombineVariants {} -o Combine_Variants.vcf -genotypeMergeOptions UNIQUIFY' \
            .format(self.java_path, self.gatk_path, self.option("ref_fasta").prop["path"], samples_option)
        command = self.add_command("combine_variants", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行combine_variants结束")
        else:
            self.set_error("运行combine_variants出错", code="33709818")

    def gatk_indel_run(self):
        """
        java -Xmx50G -jar gatk-package-4.1.4.1-local.jar -T SelectVariants  -R ref.fa  -V pop.variant.vcf -selectType INDEL
        -o pop.indel.vcf -log pop.selectINDEL.log --maxIndelSize 10 --selectTypeToExclude SYMBOLIC
        """
        cmd = "{}java -Xmx50G -jar {}gatk-package-4.1.4.1-local.jar SelectVariants -R {} -V {} --select-type-to-include INDEL --output {} " \
              "--select-type-to-exclude SYMBOLIC --max-indel-size 10" \
            .format(self.java_path, self.gatk_path, self.option("ref_fasta").prop['path'], self.vcf_path,
                    os.path.join(self.output_dir, "pop.indel.vcf"), os.path.join(self.work_dir, "pop.selectINDEL.log"))
        self.logger.info(cmd)
        self.logger.info("开始进行gatk_indel")
        command = self.add_command("gatk_indel", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gatk_indel完成！")
        else:
            self.set_error("gatk_indel出错！", code="33709819")
            self.set_error("gatk_indel出错！", code="33709820")

    def gatk_indel_filter_run(self):
        """
        java -Xmx50G -jar gatk-package-4.1.4.1-local.jar -T VariantFiltration -R ref.fa -V pop.indel.vcf  --filterExpression
        "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || InbreedingCoef < -0.5 ||SQR > 10.0" --filterName Failer
        -o pop.indel.filter.vcf -log pop.indel.filter.log
        :return:
        """
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/WGS/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "var_filter_{}.sh".format(now_time)
        cmd1 = "java -Xmx50G -jar {}gatk-package-4.1.4.1-local.jar VariantFiltration -R {} -V {} --filter-expression {} " \
               "--filter-name Failer --output {} " \
            .format(self.gatk_path, self.option("ref_fasta").prop['path'],
                    os.path.join(self.output_dir, "pop.indel.vcf"),
                    '"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || InbreedingCoef < -0.5 ||SQR > 10.0"',
                    os.path.join(self.output_dir, "pop.indel.filter.vcf"),
                    os.path.join(self.work_dir, "pop.indel.filter.log"))
        self.logger.info(cmd1)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd1)
        code = os.system('chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改%s为可执行文件失败！", variables=(file_path), code="33709821")
        shell = "/bioinfo/WGS/script_temp/{}".format(os.path.basename(file_path))
        self.logger.info("开始进行indel_filter")
        command1 = self.add_command("indel_filter", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("indel_filter完成！")
        else:
            self.set_error("indel_filter出错！", code="33709822")
            self.set_error("indel_filter出错！", code="33709823")
        os.system('rm {}'.format(file_path))

    def gatk_indel_filter_recode_run(self):
        """
        java -Xmx50G -jar gatk-package-4.1.4.1-local.jar -T SelectVariants  -R ref.fa  -V pop.indel.filter.vcf
        -o pop.indel.filter.recode.vcf --setFilteredGtToNocall --excludeFiltered --excludeNonVariants
        -log pop.selectINDEL.log
        :return:
        """
        cmd2 = "{}java -Xmx50G -jar {}gatk-package-4.1.4.1-local.jar SelectVariants -R {} -V {} --output {} " \
               "--set-filtered-gt-to-nocall --exclude-filtered --exclude-non-variants " \
            .format(self.java_path, self.gatk_path, self.option("ref_fasta").prop['path'],
                    os.path.join(self.output_dir, "pop.indel.filter.vcf"),
                    os.path.join(self.output_dir, "pop.indel.filter.recode.vcf"),
                    os.path.join(self.work_dir, "pop.selectINDEL.log"))
        self.logger.info(cmd2)
        self.logger.info("开始进行filter_recode")
        command2 = self.add_command("filter_recode_indel", cmd2).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("filter_recode完成！")
        else:
            self.set_error("filter_recode出错！", code="33709824")
            self.set_error("filter_recode出错！", code="33709825")

    def gatk_snp_run(self):
        """
        java -Xmx50G -jar gatk-package-4.1.4.1-local.jar -T SelectVariants
        -R ref.fa -V pop.variant.vcf  -selectType SNP  -o pop.snp.vcf -log pop.selectSNP.log --selectTypeToExclude
        SYMBOLIC
        """
        cmd = "{}java -Xmx50G -jar {}gatk-package-4.1.4.1-local.jar SelectVariants -R {} -V {} --select-type-to-include SNP --output {} " \
              "--select-type-to-exclude SYMBOLIC" \
            .format(self.java_path, self.gatk_path, self.option("ref_fasta").prop['path'], self.vcf_path,
                    os.path.join(self.output_dir, "pop.snp.vcf"), os.path.join(self.work_dir, "pop.selectSNP.log"))
        self.logger.info(cmd)
        self.logger.info("开始进行gatk_snp")
        command = self.add_command("gatk_snp", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gatk_snp完成！")
        else:
            self.set_error("gatk_snp出错！", code="33709826")
            self.set_error("gatk_snp出错！", code="33709827")

    def gatk_snp_filter_run(self):
        """
        java -Xmx50G -jar /mnt/ilustre/users/dna/.env/bin/gatk-package-4.1.4.1-local.jar -T VariantFiltration
        -R ref.fa -V pop.snp.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 ||
        ReadPosRankSum < -8.0 || SQR > 3.0" --filterName Failer -o pop.snp.filter.vcf -log pop.snp.filter.log
        :return:
        """
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/WGS/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "var_filter_{}.sh".format(now_time)
        cmd1 = "java -Xmx50G -jar {}gatk-package-4.1.4.1-local.jar VariantFiltration -R {} -V {} --filter-expression {} " \
               "--filter-name Failer --output {}" \
            .format(self.gatk_path, self.option("ref_fasta").prop['path'], os.path.join(self.output_dir, "pop.snp.vcf"),
                    '"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SQR > 3.0"',
                    os.path.join(self.output_dir, "pop.snp.filter.vcf"),
                    os.path.join(self.work_dir, "pop.snp.filter.log"))
        self.logger.info(cmd1)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd1)
        code = os.system('chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改%s为可执行文件失败！", variables=(file_path), code="33709828")
        shell = "/bioinfo/WGS/script_temp/{}".format(os.path.basename(file_path))
        self.logger.info("开始进行snp_filter")
        command1 = self.add_command("snp_filter", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("snp_filter完成！")
        else:
            self.set_error("snp_filter出错！", code="33709829")
            self.set_error("snp_filter出错！", code="33709830")
        os.system('rm {}'.format(file_path))

    def gatk_snp_filter_recode_run(self):
        """
        java -Xmx50G -jar /mnt/ilustre/users/dna/.env/bin/gatk-package-4.1.4.1-local.jar -T SelectVariants -R ref.fa
        -V pop.snp.filter.vcf -o pop.snp.filter.recode.vcf --setFilteredGtToNocall --excludeFiltered
        --excludeNonVariants -log pop.filterSNP.log
        :return:
        """
        cmd2 = "{}java -Xmx50G -jar {}gatk-package-4.1.4.1-local.jar SelectVariants -R {} -V {} --output {} " \
               "--set-filtered-gt-to-nocall --exclude-filtered --exclude-non-variants" \
            .format(self.java_path, self.gatk_path, self.option("ref_fasta").prop['path'],
                    os.path.join(self.output_dir, "pop.snp.filter.vcf"),
                    os.path.join(self.output_dir, "pop.snp.filter.recode.vcf"),
                    os.path.join(self.work_dir, "pop.filterSNP.log"))
        self.logger.info(cmd2)
        self.logger.info("开始进行filter_recode")
        command2 = self.add_command("filter_recode_snp", cmd2).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("filter_recode完成！")
        else:
            self.set_error("filter_recode出错！", code="33709831")
            self.set_error("filter_recode出错！", code="33709832")

    def combine_snp_indel(self):
        self.snp_vcf=os.path.join(self.output_dir,"pop.snp.filter.vcf")
        self.indel_vcf = os.path.join(self.output_dir, "pop.indel.filter.vcf")
        self.findel_vcf = os.path.join(self.output_dir, "final.pop.indel.filter.vcf")
        with open(self.indel_vcf,"r") as indel,open(self.findel_vcf,"w") as findel:
            headList = []
            for line in indel.readlines():
                if line.startswith("#"):
                    if not line.startswith("#CHROM"):  # 这个主要是一定要把#CHROM开始的这一行划归为后面的文件中
                        headList.append(line)
                        findel.write(line)
                    else:
                        findel.write("##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n")
                        findel.write(line)
                    # findel.write(line)
                else:
                    details=line.strip().split("\t")
                    INFO=details[7]
                    details[7]="INDEL;"+INFO
                    findel.write("\t".join(details)+"\n")
        cmd = "{}java -Xmx50G -jar {}gatk-package-4.1.4.1-local.jar MergeVcfs -I {} -I {} -O {} ".format(self.java_path, self.gatk_path,os.path.join(self.output_dir,"pop.snp.filter.vcf"),os.path.join(self.output_dir,"final.pop.indel.filter.vcf"),os.path.join(self.output_dir, "final.norecode.vcf"))
        self.logger.info(cmd)
        self.logger.info("开始进行vcf_merge")
        command1 = self.add_command("vcf_merge", cmd).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("vcf_merge完成！")
        else:
            self.set_error("vcf_merge出错！")
            self.set_error("vcf_merge出错！")
        # with open(os.path.join(self.output_dir,"final.vcf"),"w") as f:
        #     with open (self.snp_vcf,"r") as sn:
        #         with open (self.findel_vcf,"r") as ind:
        #             for i in sn.readlines():
        #                 f.write(i)
        #             for i in ind.readlines():
        #                 if not i.startswith("#"):
        #                     f.write(i)

    def gatk_filter_recode_run(self):
        cmd = "{}java -Xmx50G -jar {}gatk-package-4.1.4.1-local.jar SelectVariants -R {} -V {} --output {} " \
               "--set-filtered-gt-to-nocall --exclude-filtered --exclude-non-variants" \
            .format(self.java_path, self.gatk_path, self.option("ref_fasta").prop['path'],
                    os.path.join(self.output_dir, "final.norecode.vcf"),
                    os.path.join(self.output_dir, "final.vcf"),
                   )
        self.logger.info(cmd)
        self.logger.info("开始进行filter_recode")
        command2 = self.add_command("filter_recode_snp_indel", cmd).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("filter_recode完成！")
        else:
            self.set_error("filter_recode出错！")
            self.set_error("filter_recode出错！")

    def run(self):
        super(GatkToolsTool, self).run()
        if self.option("targer") == "merge_vcfs":
            self.merge_vcfs()

        self.dict()
        self.vcf_path = self.option("input_file").prop["path"]
        if self.option('combine_vcf'):
            self.combine_vcf()
            self.vcf_path = "Combine_Variants.vcf"
        self.gatk_indel_run()
        self.gatk_indel_filter_run()
        # self.gatk_indel_filter_recode_run()
        self.gatk_snp_run()
        self.gatk_snp_filter_run()
        # self.gatk_snp_filter_recode_run()
        self.combine_snp_indel()
        self.gatk_filter_recode_run()
        self.end()


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
            "id": "10w2100wfilter" + str(random.randint(1, 10000))+"-yyyyyy",
            "type": "tool",
            "name": "medical_transcriptome.somatic.gatk_tools",
            "instant": False,
            "options": dict(
              vcf_list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/vcf_list",
              out_vcf = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test1026/final.vcf"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180410

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import datetime
import random
import os


class GatkIndelAgent(Agent):
    """
    gatk snp 筛选出来
    """
    def __init__(self, parent):
        super(GatkIndelAgent, self).__init__(parent)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "pop_var_vcf", "type": "string"}   # pop.variant.vcf
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
        if not self.option("ref_fasta").is_set:
            raise OptionError("缺少ref_fasta参数", code="34502901")
        if not self.option("pop_var_vcf"):
            raise OptionError("缺少pop_var_vcf参数", code="34502902")
        ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        ref_file_name = os.path.basename(self.option("ref_fasta").prop['path']).split(".")[0]
        if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
                or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
            raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="34502903")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '55G'
        
    def end(self):
        super(GatkIndelAgent, self).end()


class GatkIndelTool(Tool):
    def __init__(self, config):
        super(GatkIndelTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.gatk_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/GenomeAnalysisTK-3.8-0-ge9d806836/"
        self.java_path = "/program/sun_jdk1.8.0/bin/"
        
    def gatk_indel_run(self):
        """
        java -Xmx50G -jar GenomeAnalysisTK.jar -T SelectVariants  -R ref.fa  -V pop.variant.vcf -selectType INDEL
        -o pop.indel.vcf -log pop.selectINDEL.log --maxIndelSize 10 --selectTypeToExclude SYMBOLIC
        """
        cmd = "{}java -Xmx50G -jar {}GenomeAnalysisTK.jar -T SelectVariants -R {} -V {} -selectType INDEL -o {} " \
              "-log {} --selectTypeToExclude SYMBOLIC --maxIndelSize 10"\
            .format(self.java_path, self.gatk_path, self.option("ref_fasta").prop['path'], self.option("pop_var_vcf"),
                    os.path.join(self.output_dir, "pop.indel.vcf"), os.path.join(self.work_dir, "pop.selectINDEL.log"))
        self.logger.info(cmd)
        self.logger.info("开始进行gatk_indel")
        command = self.add_command("gatk_indel", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gatk_indel完成！")
        else:
            self.set_error("gatk_indel出错！", code="34502901")
            self.set_error("gatk_indel出错！", code="34502911")
            
    def gatk_indel_filter_run(self):
        """
        java -Xmx50G -jar GenomeAnalysisTK.jar -T VariantFiltration -R ref.fa -V pop.indel.vcf  --filterExpression
        "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || InbreedingCoef < -0.5 ||SQR > 10.0" --filterName Failer
        -o pop.indel.filter.vcf -log pop.indel.filter.log
        :return:
        """
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/WGS/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "var_filter_{}.sh".format(now_time)
        cmd1 = "java -Xmx50G -jar {}GenomeAnalysisTK.jar -T VariantFiltration -R {} -V {} --filterExpression {} " \
               "--filterName Failer -o {} -log {}" \
            .format(self.gatk_path, self.option("ref_fasta").prop['path'],
                    os.path.join(self.output_dir, "pop.indel.vcf"),
                    '"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || InbreedingCoef < -0.5 ||SQR > 10.0"',
                    os.path.join(self.output_dir, "pop.indel.filter.vcf"),
                    os.path.join(self.work_dir, "pop.indel.filter.log"))
        self.logger.info(cmd1)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd1)
        code = os.system('/bin/chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改%s为可执行文件失败！",variables=(file_path), code="34502902")
        shell = "/bioinfo/WGS/script_temp/{}".format(os.path.basename(file_path))
        self.logger.info("开始进行indel_filter")
        command1 = self.add_command("indel_filter", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("indel_filter完成！")
        else:
            self.set_error("indel_filter出错！", code="34502903")
            self.set_error("indel_filter出错！", code="34502912")
        os.system('rm {}'.format(file_path))

    def gatk_indel_filter_recode_run(self):
        """
        java -Xmx50G -jar GenomeAnalysisTK.jar -T SelectVariants  -R ref.fa  -V pop.indel.filter.vcf
        -o pop.indel.filter.recode.vcf --setFilteredGtToNocall --excludeFiltered --excludeNonVariants
        -log pop.selectINDEL.log
        :return:
        """
        cmd2 = "{}java -Xmx50G -jar {}GenomeAnalysisTK.jar -T SelectVariants -R {} -V {} -o {} " \
               "--setFilteredGtToNocall --excludeFiltered --excludeNonVariants -log {}"\
            .format(self.java_path, self.gatk_path, self.option("ref_fasta").prop['path'],
                    os.path.join(self.output_dir, "pop.indel.filter.vcf"),
                    os.path.join(self.output_dir, "pop.indel.filter.recode.vcf"),
                    os.path.join(self.work_dir, "pop.selectINDEL.log"))
        self.logger.info(cmd2)
        self.logger.info("开始进行filter_recode")
        command2 = self.add_command("filter_recode", cmd2).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("filter_recode完成！")
        else:
            self.set_error("filter_recode出错！", code="34502904")
            self.set_error("filter_recode出错！", code="34502913")

    def run(self):
        super(GatkIndelTool, self).run()
        self.gatk_indel_run()
        self.gatk_indel_filter_run()
        self.gatk_indel_filter_recode_run()
        self.end()

# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180410

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import datetime
import random
import os


class GatkSnpAgent(Agent):
    """
    gatk snp 筛选出来
    """
    def __init__(self, parent):
        super(GatkSnpAgent, self).__init__(parent)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "pop_var_vcf", "type": "string"}
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
            raise OptionError("缺少ref_fasta参数", code="34503001")
        if not self.option("pop_var_vcf"):
            raise OptionError("缺少pop_var_vcf参数", code="34503002")
        ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        ref_file_name = os.path.basename(self.option("ref_fasta").prop['path']).split(".")[0]
        if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
                or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
            raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="34503003")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '55G'
        
    def end(self):
        super(GatkSnpAgent, self).end()


class GatkSnpTool(Tool):
    def __init__(self, config):
        super(GatkSnpTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin')
        self.gatk_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/GenomeAnalysisTK-3.8-0-ge9d806836/"
        self.java_path = "/program/sun_jdk1.8.0/bin/"
        
    def gatk_snp_run(self):
        """
        java -Xmx50G -jar GenomeAnalysisTK.jar -T SelectVariants
        -R ref.fa -V pop.variant.vcf  -selectType SNP  -o pop.snp.vcf -log pop.selectSNP.log --selectTypeToExclude
        SYMBOLIC
        """
        cmd = "{}java -Xmx50G -jar {}GenomeAnalysisTK.jar -T SelectVariants -R {} -V {} -selectType SNP -o {} " \
              "-log {} --selectTypeToExclude SYMBOLIC"\
            .format(self.java_path, self.gatk_path, self.option("ref_fasta").prop['path'], self.option("pop_var_vcf"),
                    os.path.join(self.output_dir, "pop.snp.vcf"), os.path.join(self.work_dir, "pop.selectSNP.log"))
        self.logger.info(cmd)
        self.logger.info("开始进行gatk_snp")
        command = self.add_command("gatk_snp", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gatk_snp完成！")
        else:
            self.set_error("gatk_snp出错！", code="34503001")
            self.set_error("gatk_snp出错！", code="34503011")
            
    def gatk_snp_filter_run(self):
        """
        java -Xmx50G -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T VariantFiltration
        -R ref.fa -V pop.snp.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 ||
        ReadPosRankSum < -8.0 || SQR > 3.0" --filterName Failer -o pop.snp.filter.vcf -log pop.snp.filter.log
        :return:
        """
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/WGS/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "var_filter_{}.sh".format(now_time)
        cmd1 = "java -Xmx50G -jar {}GenomeAnalysisTK.jar -T VariantFiltration -R {} -V {} --filterExpression {} " \
               "--filterName Failer -o {} -log {}" \
            .format(self.gatk_path, self.option("ref_fasta").prop['path'], os.path.join(self.output_dir, "pop.snp.vcf"),
                    '"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SQR > 3.0"',
                    os.path.join(self.output_dir, "pop.snp.filter.vcf"),
                    os.path.join(self.work_dir, "pop.snp.filter.log"))
        self.logger.info(cmd1)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd1)
        code = os.system('/bin/chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改%s为可执行文件失败！", variables=(file_path), code="34503002")
        shell = "/bioinfo/WGS/script_temp/{}".format(os.path.basename(file_path))
        self.logger.info("开始进行snp_filter")
        command1 = self.add_command("snp_filter", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("snp_filter完成！")
        else:
            self.set_error("snp_filter出错！", code="34503003")
            self.set_error("snp_filter出错！", code="34503012")
        os.system('rm {}'.format(file_path))

    def gatk_snp_filter_recode_run(self):
        """
        java -Xmx50G -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T SelectVariants -R ref.fa
        -V pop.snp.filter.vcf -o pop.snp.filter.recode.vcf --setFilteredGtToNocall --excludeFiltered
        --excludeNonVariants -log pop.filterSNP.log
        :return:
        """
        cmd2 = "{}java -Xmx50G -jar {}GenomeAnalysisTK.jar -T SelectVariants -R {} -V {} -o {} " \
               "--setFilteredGtToNocall --excludeFiltered --excludeNonVariants -log {}"\
            .format(self.java_path, self.gatk_path, self.option("ref_fasta").prop['path'],
                    os.path.join(self.output_dir, "pop.snp.filter.vcf"),
                    os.path.join(self.output_dir, "pop.snp.filter.recode.vcf"),
                    os.path.join(self.work_dir, "pop.filterSNP.log"))
        self.logger.info(cmd2)
        self.logger.info("开始进行filter_recode")
        command2 = self.add_command("filter_recode", cmd2).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("filter_recode完成！")
        else:
            self.set_error("filter_recode出错！", code="34503004")
            self.set_error("filter_recode出错！", code="34503013")

    def run(self):
        super(GatkSnpTool, self).run()
        self.gatk_snp_run()
        self.gatk_snp_filter_run()
        self.gatk_snp_filter_recode_run()
        self.end()

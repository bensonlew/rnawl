# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180403

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class HaplotypeAgent(Agent):
    """
    SNP工具haplotypeCaller功能
    """
    def __init__(self, parent):
        super(HaplotypeAgent, self).__init__(parent)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample_bam", "type": "infile", "format": "align.bwa.bam"},
            {"name": "sample_id", "type": "string"}  # 样本名编号, 样本名字中不能有.
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
            raise OptionError("缺少ref_fasta参数", code="34503601")
        if not self.option("sample_bam").is_set:
            raise OptionError("缺少sample_bam参数", code="34503602")
        if not self.option("sample_id"):
            raise OptionError("缺少sample_id参数", code="34503603")
        ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        ref_file_name = os.path.basename(self.option("ref_fasta").prop['path']).split(".")[0]
        if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
                or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
            raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="34503604")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = '100G'
        
    def end(self):
        super(HaplotypeAgent, self).end()


class HaplotypeTool(Tool):
    def __init__(self, config):
        super(HaplotypeTool, self).__init__(config)
        self.gatk_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/GenomeAnalysisTK-3.8-0-ge9d806836/"
        self.java_path = "/program/sun_jdk1.8.0/bin/"
        
    def gatk_vc(self):
        """
        step：gatk variant calling
        /mnt/ilustre/users/sanger-dev/app/program/sun_jdk1.8.0/bin/java -XX:+UseSerialGC -Xmx120G
        -Djava.io.tmpdir=/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/test/tmp/
        -jar /mnt/ilustre/users/sanger-dev/app/bioinfo/gene-structure/GenomeAnalysisTK.jar
        -T HaplotypeCaller -R /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/ref.fa
        -I /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/A8-10.mkdup.bam
        -o /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/test/A8-10.g.vcf
        -nct 8 --genotyping_mode DISCOVERY --emitRefConfidence GVCF
        -stand_call_conf 30 -variant_index_type LINEAR -variant_index_parameter 128000
        -filterNoBases -filterMBQ -filterRNC -dontUseSoftClippedBases
        -log /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/test/A8-10.gvcf.log
        """
        vcf_name = "{}.g.vcf".format(self.option("sample_id"))
        log_name = "{}.gvcf.log".format(self.option("sample_id"))
        cmd = "{}java -XX:+UseSerialGC -Xmx100G -Djava.io.tmpdir={}/tmp/ -jar {}GenomeAnalysisTK.jar -T HaplotypeCaller" \
              " -R {} -I {} -o {} -nct 8 --genotyping_mode DISCOVERY --emitRefConfidence GVCF -stand_call_conf 30 " \
              "-variant_index_type LINEAR -variant_index_parameter 128000 -filterNoBases -filterMBQ -filterRNC " \
              "-dontUseSoftClippedBases -log {}".format(self.java_path, self.work_dir, self.gatk_path,
                                                        self.option("ref_fasta").prop['path'],
                                                        self.option("sample_bam").prop['path'],
                                                        os.path.join(self.output_dir, vcf_name), log_name)
        self.logger.info(cmd)
        self.logger.info("开始进行HaplotypeCaller")
        command = self.add_command("haplotypecaller", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("HaplotypeCaller完成！")
        elif command.return_code == 1:
            if self.check_rerun():
                self.logger.info("返回码为1，重新运行一次")
                command.rerun()
                self.wait()
                if command.return_code == 0:
                    self.logger.info("重新运行一次成功！")
                else:
                    self.set_error("HaplotypeCaller出错！", code="34503601")
                    self.set_error("HaplotypeCaller出错！", code="34503607")
            else:
                self.set_error("HaplotypeCaller出错！", code="34503602")
                self.set_error("HaplotypeCaller出错！", code="34503608")
            
    def check_rerun(self):
        """
        "##### ERROR MESSAGE: Bad input: We encountered a non-standard non-IUPAC base in the provided reference: '10'"
        gatk 会随机报错，暂时找不到原因，所以就在这边对haplotypecaller.o结果进行检查如果有non-IUPAC base in the provided
        reference就重新运行一次
        :return:
        """
        with open(self.work_dir + "/haplotypecaller.o", 'r') as r:
            data = r.readlines()
            for line in data:
                if re.match(r".*non-IUPAC base in the provided reference.*", line):
                    self.logger.info("error：{}".format(line))
                    return True
        return False

    def run(self):
        super(HaplotypeTool, self).run()
        self.gatk_vc()
        self.end()

# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180403

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class SamtoolsCallAgent(Agent):
    """
    SNP工具使用Samtools进行call snp 首先通过samtools mpileup生成bcf文件，然后后面用bcftools call进行call snp与indel
    """
    def __init__(self, parent):
        super(SamtoolsCallAgent, self).__init__(parent)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "bed_file", "type": "string"},  # 切割后的bed文件
            {"name": "bam_list", "type": "string"}  # 所有样本的bam文件列表
        ]
        self.add_option(options)
        self.step.add_steps('call_snp')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.call_snp.start()
        self.step.update()

    def step_end(self):
        self.step.call_snp.finish()
        self.step.update()
        
    def check_options(self):
        if not self.option("ref_fasta").is_set:
            raise OptionError("缺少ref_fasta参数", code="34504801")
        if not self.option("bed_file"):
            raise OptionError("缺少bed_file参数", code="34504802")
        if not self.option("bam_list"):
            raise OptionError("缺少bam_list参数", code="34504803")
        ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        ref_file_name = os.path.basename(self.option("ref_fasta").prop['path']).split(".")[0]
        if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
                or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
            raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="34504804")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = '100G'
        
    def end(self):
        super(SamtoolsCallAgent, self).end()


class SamtoolsCallTool(Tool):
    def __init__(self, config):
        super(SamtoolsCallTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/align/samtools-1.4/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/bcftools-1.4/bin')
        self.samtools = 'miniconda2/bin/'
        self.bcftools = "bioinfo/seq/bcftools-1.7/"
        
    def samtools_mpileup(self):
        """
                samtools mpileup  -t DP,AD
                -uf /mnt/ilustre/users/qingmei.cui/newmdt/Project/BSA_Test_TAIR10/2018.3.2.var/02.ref-config/ref.fa
                -l /mnt/ilustre/users/qingmei.cui/newmdt/Project/BSA_Test_TAIR10/2018.3.2.var/step07.variant-
                call-samtools/1.bed -b /mnt/ilustre/users/qingmei.cui/newmdt/Project/BSA_Test_TAIR10/2018.3.2.v
                ar/step07.variant-call-samtools/bam.list|bcftools call -mv --format-fields GQ,GP
                --output-type z > /mnt/ilustre/users/qingmei.cui/newmdt/Project/BSA_Test_TAIR10/2018.3.2.
                var/step07.variant-call-samtools/1.vcf.gz
        """
        cmd = "{}samtools mpileup -t DP,AD -uf {} -l {} -b {} -o {}.bcf"\
            .format(self.samtools, self.option("ref_fasta").prop['path'], self.option("bed_file"),
                    self.option("bam_list"), os.path.join(self.work_dir, os.path.basename(self.option('bed_file'))))
        self.logger.info(cmd)
        self.logger.info("开始进行mpileup")
        command = self.add_command("cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("mpileup完成！")
        else:
            self.set_error("mpileup出错！", code="34504801")
            self.set_error("mpileup出错！", code="34504807")

    def bcftools_call(self):
        cmd1 = "{}bcftools call -mv --format-fields GQ,GP --output-type z {}.bcf -o {}.vcf.gz"\
            .format(self.bcftools, os.path.join(self.work_dir, os.path.basename(self.option('bed_file'))),
                    os.path.join(self.output_dir, os.path.basename(self.option("bed_file"))))
        self.logger.info(cmd1)
        self.logger.info("开始进行bcftools_call")
        command1 = self.add_command("cmd1", cmd1).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("bcftools_call完成！")
        else:
            self.set_error("bcftools_call出错！", code="34504802")
            self.set_error("bcftools_call出错！", code="34504808")

    def run(self):
        super(SamtoolsCallTool, self).run()
        self.samtools_mpileup()
        self.bcftools_call()
        self.end()

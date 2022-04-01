# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180404

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import time


class GvcfTypingAgent(Agent):
    """
    SNP工具GenotypeGVCFs功能, 将gvcf转换成vcf，这一步按照单个样品进行计算，先计算每个样本的vcf，然后后面在合并
    """
    def __init__(self, parent):
        super(GvcfTypingAgent, self).__init__(parent)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "gvcf_path", "type": "string"},  # gvcf文件路径，每个样本就有对应的gvcf文件
            {"name": "intervals_file", "type": "string"}  # 该文件是将基因组等分成25份
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
            raise OptionError("缺少ref_fasta参数", code="34503501")
        if not self.option("gvcf_path"):
            raise OptionError("缺少gvcf_path参数", code="34503502")
        if not self.option("intervals_file"):
            raise OptionError("缺少intervals_file参数", code="34503503")
        ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        ref_file_name = os.path.basename(self.option("ref_fasta").prop['path']).split(".")[0]
        # if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
        #         or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
        #     raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="34503504")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 16
        self._memory = '120G'
        
    def end(self):
        super(GvcfTypingAgent, self).end()


class GvcfTypingTool(Tool):
    def __init__(self, config):
        super(GvcfTypingTool, self).__init__(config)
        self.gatk_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/GenomeAnalysisTK-3.8-0-ge9d806836/"
        self.java_path = "/program/sun_jdk1.8.0/bin/"
        
    def gatk_vc(self):
        """
        step：gvcf转vcf
        /mnt/ilustre/users/sanger-dev/app/program/sun_jdk1.8.0/bin/java -Djava.io.tmpdir=$dOut/tmp/ -Xmx120G
        -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T GenotypeGVCFs -V $dOut/gvcf.list
         -o $dOut/pop.noid.vcf -R $ref -log $dOut/pop.variant.log -jdk_deflater -jdk_inflater( 后面两个jdk使用上会报错)
        """
        vcf_name = os.path.basename(self.option("intervals_file")).split(".")[0]
        cmd = "{}java -Djava.io.tmpdir={}/tmp/ -Xmx100G -jar {}GenomeAnalysisTK.jar -T GenotypeGVCFs -V {}/gvcf.list" \
              " -o {}/{}.noid.vcf -R {} -log {}/pop.variant.log -L {} --never_trim_vcf_format_field " \
              "--max_alternate_alleles 3 --disable_auto_index_creation_and_locking_when_reading_rods -jdk_inflater " \
              "-jdk_deflater -nt 15"\
            .format(self.java_path, self.output_dir, self.gatk_path, self.output_dir, self.output_dir, vcf_name,
                    self.option("ref_fasta").prop['path'], self.work_dir, self.option("intervals_file"))
        self.logger.info(cmd)
        self.logger.info("开始进行GenotypeGVCFs")
        # command = self.add_command("cmd", cmd).run()
        command = self.add_command("cmd", cmd, ignore_error=True).run()  # 上线后使用
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("GenotypeGVCFs完成！")
        elif command.return_code == 1:
            if self.check_rerun():
                self.logger.info("返回码为1，重新运行一次")
                command.rerun()
                self.wait()
                if command.return_code == 0:
                    self.logger.info("重新运行一次成功！")
                else:
                    self.set_error("GenotypeGVCFs出错！", code="34503501")
            else:
                self.set_error("GenotypeGVCFs出错！", code="34503502")

    def check_rerun(self):
        """
        "##### ERROR MESSAGE: Bad input: We encountered a non-standard non-IUPAC base in the provided reference: '10'"
        gatk 会随机报错，暂时找不到原因，所以就在这边对haplotypecaller.o结果进行检查如果有non-IUPAC base in the provided
        reference就重新运行一次
        :return:
        """
        with open(self.work_dir + "/cmd.o", 'r') as r:
            data = r.readlines()
            for line in data:
                if re.match(r".*non-IUPAC base in the provided reference.*", line):
                    self.logger.info("error：{}".format(line))
                    return True
        return False

    def set_gvcf_list(self):
        time.sleep(2)
        files = os.listdir(self.option("gvcf_path"))
        self.logger.info("开始设定gvcf文件列表！")
        if len(files) == 0:
            self.set_error("gvcf列表为空，不能进行后续分析！", code="34503503")
        else:
            gvcf = os.path.join(self.output_dir, "gvcf.list")
            with open(gvcf, 'w') as w:
                for m in files:
                    if re.match(r".*\.g\.vcf$", m):
                        w.write("{}\n".format(os.path.join(self.option("gvcf_path"), m)))
            self.logger.info("设定gvcf文件列表成功！")

    def run(self):
        super(GvcfTypingTool, self).run()
        self.set_gvcf_list()
        self.gatk_vc()
        self.end()

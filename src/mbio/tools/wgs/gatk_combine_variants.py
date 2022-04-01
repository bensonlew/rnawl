# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180411

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import datetime
import os
import re


class GatkCombineVariantsAgent(Agent):
    """
    gatk snp indel 注释文件合并
    """
    def __init__(self, parent):
        super(GatkCombineVariantsAgent, self).__init__(parent)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "snp_anno_primary_vcf", "type": "string"},  # snp.anno.primary.vcf
            {"name": "indel_anno_primary_vcf", "type": "string"}   # indel.anno.primary.vcf
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
            raise OptionError("缺少ref_fasta参数", code="34502801")
        if not self.option("snp_anno_primary_vcf"):
            raise OptionError("缺少snp_anno_primary_vcf参数", code="34502802")
        if not self.option("indel_anno_primary_vcf"):
            raise OptionError("缺少indel_anno_primary_vcf参数", code="34502803")
        ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        ref_file_name = os.path.basename(self.option("ref_fasta").prop['path']).split(".")[0]
        if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
                or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
            raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="34502804")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '55G'
        
    def end(self):
        super(GatkCombineVariantsAgent, self).end()


class GatkCombineVariantsTool(Tool):
    def __init__(self, config):
        super(GatkCombineVariantsTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin')
        self.gatk_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/GenomeAnalysisTK-3.8-0-ge9d806836/"
        self.java_path = "/program/sun_jdk1.8.0/bin/"
        
    def gatk_combine_run(self):
        """
        java -Djava.io.tmpdir=/tmp/ -Xmx50G -jar GenomeAnalysisTK.jar -T CombineVariants -R ref.fa  -V
        snp.anno.primary.vcf  -V indel.anno.primary.vcf  -o pop.final.vcf
        --genotypemergeoption UNSORTED -log pop.merge.log
        """
        cmd = "{}java -Djava.io.tmpdir={} -Xmx50G -jar {}GenomeAnalysisTK.jar -T CombineVariants -R {} -V {} " \
              "-V {} -o {} --genotypemergeoption UNSORTED -log {}"\
            .format(self.java_path, os.path.join(self.work_dir, "tmp"), self.gatk_path,
                    self.option("ref_fasta").prop['path'], self.option("snp_anno_primary_vcf"),
                    self.option("indel_anno_primary_vcf"), os.path.join(self.output_dir, "pop.final.vcf"),
                    os.path.join(self.work_dir, "pop.merge.log"))
        self.logger.info(cmd)
        self.logger.info("开始进行gatk_combine")
        command = self.add_command("gatk_combine", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gatk_combine完成！")
        else:
            self.set_error("gatk_combine出错！", code="34502801")
            self.set_error("gatk_combine出错！", code="34502805")

    def gzip_file(self):
        self.logger.info("开始压缩文件！")
        code = os.system("gzip {}".format(os.path.join(self.output_dir, "pop.final.vcf")))
        if code == 0:
            self.logger.info("压缩文件成功！")
        else:
            self.set_error("压缩文件失败！", code="34502802")

    def get_snp_indel_num(self):
        indel = self.get_file_lines(self.option("indel_anno_primary_vcf"))
        snp = self.get_file_lines(self.option("snp_anno_primary_vcf"))
        with open(self.output_dir + '/snp_indel_num.txt', 'w') as w:
            w.write('#snp_num\tindel_num\n')
            w.write('{}\t{}\n'.format(snp, indel))

    def get_file_lines(self, files):
        num = 0
        with open(files, 'r') as r:
            for line in r:
                if not re.match('#.*', line):
                    num += 1
        return num

    def run(self):
        super(GatkCombineVariantsTool, self).run()
        self.gatk_combine_run()
        self.get_snp_indel_num()
        # self.gzip_file()
        self.end()

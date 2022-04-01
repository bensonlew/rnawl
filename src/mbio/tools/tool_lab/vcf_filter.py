# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class VcfFilterAgent(Agent):
    """
    群体进化，vcftool过滤vcf文件
    """
    def __init__(self, parent):
        super(VcfFilterAgent, self).__init__(parent)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "max_missing", "type": "float"},   # --max-missing
            {"name": "maf", "type": "float"},       # --maf最小过滤值
            {"name":"chr","type":"string","default":"all"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("vcf_file").is_set:
            raise OptionError("请设置VCF文件")
        if not self.option("max_missing"):
            raise OptionError("请设置最大缺失率")
        if not self.option("maf"):
            raise OptionError("请设置maf最小过滤值")

    def set_resource(self):
        """
        运行所需资源
        vcftools一般不耗费内存资源，这里暂时不设置动态内存
        """
        self._cpu = 2
        self._memory = "3G"


    def end(self):
        super(VcfFilterAgent, self).end()


class VcfFilterTool(Tool):
    def __init__(self, config):
        super(VcfFilterTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.vcftools_path = "bioinfo/dna_evolution/vcftools"

    def VcfFilter(self):
        """
        vcftools --maf $maf --max-missing $mis --recode
         --out $output --vcf $input
        :return:
        """
        cmd = "{} --vcf {}".format(self.vcftools_path, self.option('vcf_file').prop['path'])
        if self.option('maf'):
            cmd += " --maf {}".format(float(self.option('maf')))
        if self.option('max_missing'):
            cmd += " --max-missing {}".format(1 - float(self.option('max_missing')))
            cmd += " %s" % '--recode'
            cmd += " --out {}".format(self.output_dir + "/pop")
        if self.option("chr") != "all":
            for i in self.option("chr").split(";"):
                cmd += " --chr {}".format(i)
        self.logger.info(cmd)
        self.logger.info("开始进行VcfFilter")
        command = self.add_command("vcftoolsfilter", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("VcfFilter完成！")
        else:
            self.set_error("VcfFilter出错！")
        # self.option("filter_vcf", self.output_dir + "/pop.recode.vcf")
        file_path = os.path.join(self.output_dir, "pop.recode.vcf")
        if os.path.exists(file_path):
            self.logger.info("filter成功")


    def run(self):
        super(VcfFilterTool, self).run()
        self.VcfFilter()
        self.end()

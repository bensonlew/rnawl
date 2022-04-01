# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180822

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class VcftoolsFilterAgent(Agent):
    """
    群体进化，vcftool过滤vcf文件
    """
    def __init__(self, parent):
        super(VcftoolsFilterAgent, self).__init__(parent)
        options = [
            {"name": "vcf_path", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "recode", "type": "bool", "default": False},
            {"name": "remove_indels", "type": "bool", "default": False},     # --remove-indels
            {"name": "remove_filtered_all", "type": "bool", "default": False},   # --remove-filtered-all
            {"name": "plink", "type": "bool", "default": False},  # --plink
            {"name": "minDP", "type": "int"},           # --minDP 2
            {"name": "maxDP", "type": "int"},           # --maxDP 6
            {"name": "max_missing", "type": "float"},   # --max-missing
            {"name": "min_maf", "type": "float"},       # --maf 0.05
            {"name": "max_maf", "type": "float"},       # --max-maf 1
            {"name": "filter_vcf", "type": "outfile", "format": "dna_evolution.vcf"},
            # {"name": "filter_vcf", "type": "string"},
            {"name": "keep", "type": "string"},  # 本参数后面添加分组类型 By Binbin Zhao
            {"name": "group_name", "type": "string", "default": "pop"}  # 本参数用于生成结果命名，连锁不平衡中使用。
        ]
        self.add_option(options)
        self.step.add_steps('VcftoolsFilter')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.VcftoolsFilter.start()
        self.step.update()

    def step_end(self):
        self.step.VcftoolsFilter.finish()
        self.step.update()

    def check_options(self):
        if not self.option("vcf_path").is_set:
            raise OptionError("请设置vcf_path")
        if type(self.option('recode')) is not bool:
            raise OptionError("recode参数类型错误")
        if type(self.option('remove_filtered_all')) is not bool:
            raise OptionError("remove_filtered_all参数类型错误")
        if not isinstance(self.option('remove_indels'), bool):
            raise OptionError("remove_indels参数类型错误")
        if not isinstance(self.option('plink'), bool):
            raise OptionError("plink参数类型错误")

    def set_resource(self):
        """
        运行所需资源
        vcftools一般不耗费内存资源，这里暂时不设置动态内存
        """
        self._cpu = 2
        self._memory = "3G"
        # size = os.path.getsize(self.option("vcf_path").prop['path'])
        # size = size / 1024 / 1024 / 1024
        # if size < 50:
        #     self._memory = "40G"
        # elif size < 80:
        #     self._memory = "70G"
        # elif size < 100:
        #     self._memory = "90G"
        # else:
        #     self._memory = "120G"

    def end(self):
        super(VcftoolsFilterAgent, self).end()


class VcftoolsFilterTool(Tool):
    def __init__(self, config):
        super(VcftoolsFilterTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.vcftools_path = "bioinfo/dna_evolution/vcftools"

    def VcftoolsFilter(self):
        """
        要重新写下！！！
        --remove-filtered-all --remove-indels --minDP 2
        --maxDP 6 --max-missing 0.7  --maf 0.05
        --max-maf 1 --vcf  --recode  --out ./pop
        :return:
        """
        cmd = "{} --vcf {}".format(self.vcftools_path, self.option('vcf_path').prop['path'])
        if self.option('remove_filtered_all') is True:
            cmd += " %s" % '--remove-filtered-all'
        if self.option('remove_indels') is True:
            cmd += " %s" % '--remove-indels'
        if self.option('minDP'):
            cmd += " --minDP {}".format(int(self.option('minDP')))
        if self.option('maxDP'):
            cmd += " --maxDP {}".format(int(self.option('maxDP')))
        if self.option('min_maf'):
            cmd += " --maf {}".format(float(self.option('min_maf')))
        if self.option('max_maf'):
            cmd += " --max-maf {}".format(float(self.option('max_maf')))
        if self.option('max_missing'):
            cmd += " --max-missing {}".format(1 - float(self.option('max_missing')))
        if self.option('recode') is True:
            cmd += " %s" % '--recode'
        if self.option('plink') is True:
            cmd += " %s" % '--plink'
        if self.option('keep'):
            cmd += " --keep {}".format(self.option('keep'))
        if self.option('group_name'):
            cmd += " --out {}".format(os.path.join(self.output_dir, self.option("group_name")))
        else:
            cmd += " --out {}".format(self.output_dir + "/pop")
        self.logger.info(cmd)
        self.logger.info("开始进行VcftoolsFilter")
        command = self.add_command("vcftoolsfilter", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("VcftoolsFilter完成！")
        else:
            self.set_error("VcftoolsFilter出错！")
        # self.option("filter_vcf", self.output_dir + "/pop.recode.vcf")
        file_path = os.path.join(self.output_dir, self.option("group_name") + ".recode.vcf")
        if os.path.exists(file_path):
            self.logger.info("设置filter_vcf成功")
            self.option("filter_vcf", file_path)

    def run(self):
        super(VcftoolsFilterTool, self).run()
        self.VcftoolsFilter()
        self.end()

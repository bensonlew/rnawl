# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class VcftoolsFilterAgent(Agent):
    """
    群体进化，vcftool过滤vcf文件
    """
    def __init__(self, parent):
        super(VcftoolsFilterAgent, self).__init__(parent)
        options = [
            {"name": "vcf_path", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "max_missing", "type": "float", 'default': 0.3},   # --max-missing
            {"name": "maf", "type": "float", 'default': 0.05},       # --maf 0.05
            {"name": "filtered_vcf", "type": "outfile", "format": "dna_evolution.vcf"},
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
        return True

    def set_resource(self):
        """
        运行所需资源
        vcftools一般不耗费内存资源
        vcftools --plink 较耗资源
        """
        self._cpu = 8
        size = os.path.getsize(self.option("vcf_path").prop['path'])
        size = size / 1024 / 1024 / 1024
        if size < 50:
            self._memory = "40G"
        elif size < 80:
            self._memory = "70G"
        elif size < 100:
            self._memory = "90G"
        else:
            self._memory = "120G"

    def end(self):
        super(VcftoolsFilterAgent, self).end()


class VcftoolsFilterTool(Tool):
    def __init__(self, config):
        super(VcftoolsFilterTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.program = {
            'vcftools': 'bioinfo/dna_evolution/vcftools',
        }
        self.file = {
            'pop_out': os.path.join(self.output_dir, 'pop'),
            'pop_recode': os.path.join(self.output_dir, 'pop.recode.vcf')
        }

    def vcftools_filter(self):
        """
        """
        cmd = "{} --remove-filtered-all ".format(self.program['vcftools'])
        cmd += "--vcf {} ".format(self.option('vcf_path').prop['path'])
        cmd += "--out {} ".format(self.file['pop_out'])
        cmd += "--max-missing {} --maf {} ".format(self.option('max_missing'), self.option('maf'))
        cmd += "--minDP 2 --recode --recode-INFO-all"
        self.logger.info(cmd)
        self.logger.info("开始进行VcftoolsFilter")
        command = self.add_command("vcftoolsfilter", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("VcftoolsFilter完成！")
        else:
            self.set_error("VcftoolsFilter出错！")

    def vcftools_plink(self):
        """
        """
        cmd = "{} --vcf {} ".format(self.program['vcftools'], self.file['pop_recode'])
        cmd += "--plink --out {} ".format(self.file['pop_out'])
        self.logger.info(cmd)
        self.logger.info("开始进行Vcftools plink")
        command = self.add_command("vcftoolsplink", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Vcftools plink完成！")
        else:
            self.set_error("Vcftools plink出错！")
        if os.path.exists(self.file['pop_recode']):
            self.logger.info("设置filtered_vcf成功")
            self.option("filtered_vcf", self.file['pop_recode'])

    def run(self):
        super(VcftoolsFilterTool, self).run()
        self.vcftools_filter()
        self.vcftools_plink()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "selective_sweep_filter_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.selective_sweep.vcftools_filter",
            "instant": False,
            "options": dict(
                vcf_path='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/selective_sweep/final.vcf',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
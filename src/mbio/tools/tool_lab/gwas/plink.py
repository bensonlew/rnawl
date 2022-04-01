# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class PlinkAgent(Agent):
    """
    群体进化，vcftool过滤vcf文件
    """
    def __init__(self, parent):
        super(PlinkAgent, self).__init__(parent)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_evolution.vcf"},
            {'name': 'chrom_map', "type": 'infile', 'format': 'ref_rna_v2.common'},
        ]
        self.add_option(options)
        self.step.add_steps('plink')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.plink.start()
        self.step.update()

    def step_end(self):
        self.step.plink.finish()
        self.step.update()

    def check_options(self):
        if not self.option("vcf_file").is_set:
            raise OptionError("请设置vcf_file")
        return True

    def set_resource(self):
        """
        运行所需资源
        vcftools一般不耗费内存资源
        vcftools --plink 较耗资源
        """
        self._cpu = 8
        size = os.path.getsize(self.option("vcf_file").prop['path'])
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
        super(PlinkAgent, self).end()


class PlinkTool(Tool):
    def __init__(self, config):
        super(PlinkTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.program = {
            'vcftools': 'bioinfo/dna_evolution/vcftools',
            'plink': 'bioinfo/dna_evolution/plink/plink',
        }
        self.file = {
            'pop_out': os.path.join(self.output_dir, 'pop'),
        }

    def vcftools_plink(self):
        """
        """
        cmd = "{} --vcf {} ".format(self.program['vcftools'], self.option('vcf_file').prop['path'])
        cmd += "--plink --chrom-map {} ".format(self.option('chrom_map').prop['path'])
        cmd += "--out {}".format(self.file['pop_out'])
        self.logger.info(cmd)
        self.logger.info("开始进行Vcftools plink")
        command = self.add_command("vcftoolsplink", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Vcftools plink完成！")
        else:
            self.set_error("Vcftools plink出错！")
        self.run_plink()

    def run_plink(self):
        cmd = '{} --file {} '.format(self.program['plink'], self.file['pop_out'])
        cmd += '--out {} '.format(self.file['pop_out'])
        cmd += '--make-bed --allow-extra-chr'
        self.logger.info(cmd)
        self.logger.info("开始进行plink")
        command = self.add_command("plink", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("plink运行完成！")
        else:
            self.set_error("plink运行出错！")

    def run(self):
        super(PlinkTool, self).run()
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
            "id": "gwas_plink_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.gwas.plink",
            "instant": False,
            "options": dict(
                vcf_file='/mnt/ilustre/users/sanger-dev/workspace/20210519/Single_gwas_filter_3931_6429/VcftoolsFilter/output/pop.recode.vcf',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import pandas as pd


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
            {"name": "pop_recode", "type": "outfile", "format": "dna_evolution.vcf"},
            {"name": "chrom_map", "type": "outfile", "format": "ref_rna_v2.common"},
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
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(VcftoolsFilterAgent, self).end()


class VcftoolsFilterTool(Tool):
    def __init__(self, config):
        super(VcftoolsFilterTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.program = {
            'vcftools': 'bioinfo/dna_evolution/vcftools',
            'bcftools': 'bioinfo/seq/bcftools-1.7/bcftools',
        }
        self.file = {
            'pop_out': os.path.join(self.output_dir, 'pop'),
            'pop_recode': os.path.join(self.output_dir, 'pop.recode.vcf'),
            'chrom_map': os.path.join(self.output_dir, 'chrom_map.txt'),
            'chrom': os.path.join(self.work_dir, 'chrom.txt'),
        }

    def vcftools_filter(self):
        """
        """
        cmd = "{} --recode ".format(self.program['vcftools'])
        cmd += "--vcf {} ".format(self.option('vcf_path').prop['path'])
        cmd += "--out {} ".format(self.file['pop_out'])
        cmd += "--max-missing {} --maf {} ".format(self.option('max_missing'), self.option('maf'))
        self.logger.info(cmd)
        self.logger.info("开始进行VcftoolsFilter")
        command = self.add_command("vcftoolsfilter", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("VcftoolsFilter完成！")
        else:
            self.set_error("VcftoolsFilter出错！")
        if os.path.exists(self.file['pop_recode']):
            self.logger.info("设置pop_recode成功")
            self.option("pop_recode", self.file['pop_recode'])
            self.get_chrom_map()

    def get_chrom_map(self):
        """
        """
        cmd = "{} ".format(os.path.join(self.config.SOFTWARE_DIR, self.program['bcftools']))
        cmd += "view -H {} ".format(self.file['pop_recode'])
        cmd += "| cut -f 1 | sort | uniq | awk \'{print $0}\' "
        cmd += "> {}".format(self.file['chrom'])
        self.logger.info(cmd)
        self.logger.info("开始获取chromosome信息")
        command = self.add_command("get_chrom_info", cmd, shell=True)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("获取chromosome信息完成！")
        else:
            self.set_error("获取chromosome信息出错！")
        if os.path.exists(self.file['chrom']):
            with open(self.file['chrom'], 'r') as c, open(self.file['chrom_map'], 'w') as m:
                count = 1
                for line in c:
                    m.write(line.strip() + '\t' + str(count) + '\n')
                    count += 1
        if os.path.exists(self.file['chrom_map']):
            self.logger.info("设置chrom_map成功")
            self.option("chrom_map", self.file['chrom_map'])

    def run(self):
        super(VcftoolsFilterTool, self).run()
        self.vcftools_filter()
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
            "id": "gwas_filter_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.gwas.vcftools_filter",
            "instant": False,
            "options": dict(
                vcf_path='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/gwas/chr.vcf',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
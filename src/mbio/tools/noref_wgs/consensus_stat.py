# -*- coding: utf-8 -*-
# __author__ = 'wentianliu'
# last modify 20190116

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import shutil


class ConsensusStatAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(ConsensusStatAgent, self).__init__(parent)
        options = [
            {"name": "catalog_path", "type": "string"},  # 04.cstacks  catalog.tags.tsv.gz
            {"name": "sample_num", "type": "int"},  # 样本数量
            {"name": "populations_snps_vcf", "type": "string"},  # 06.genotype  populations.snps.vcf
            {"name": "ustacks_path", "type": "string"},  # 04.ustacks output
        ]
        self.add_option(options)
        self.step.add_steps('ConsensusStat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.ConsensusStat.start()
        self.step.update()

    def step_end(self):
        self.step.ConsensusStat.finish()
        self.step.update()

    def check_options(self):
        if not self.option('catalog_path'):
            raise OptionError('必须输入:catalog_path', code="35500613")
        if not os.path.exists(self.option('catalog_path')):
            raise OptionError('catalog.tags.tsv.gz文件不存在', code="35500614")
        if not self.option('sample_num'):
            raise OptionError('必须输入:sample_num', code="35500615")
        if not self.option('populations_snps_vcf'):
            raise OptionError('必须输入:populations_snps_vcf', code="35500616")
        if not os.path.exists(self.option('populations_snps_vcf')):
            raise OptionError('populations.snps.vcf文件不存在', code="35500617")
        if not self.option('ustacks_path'):
            raise OptionError('必须输入:ustacks_path', code="35500618")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 12
        self._memory = "5G"

    def end(self):
        super(ConsensusStatAgent, self).end()


class ConsensusStatTool(Tool):
    def __init__(self, config):
        super(ConsensusStatTool, self).__init__(config)
        self.Python_path = '/program/Python/bin/python'
        self.Python_all_path = self.config.SOFTWARE_DIR + '/program/Python/bin/python'
        self.consensus_stat_path = self.config.PACKAGE_DIR + "/noref_wgs/consensus_stat.py"
        self.snp_depth_density_path = self.config.PACKAGE_DIR + "/noref_wgs/snp_depth_density.py"
        self.tag_dep_path = self.config.PACKAGE_DIR + "/noref_wgs/tag_dep.py"
        self.para_fly = 'program/parafly-r2013-01-21/bin/bin/ParaFly'

    def run_consensus_stat(self):
        """
        """
        cmd = "{} {} -i {} -o {} -n {}".format(self.Python_path, self.consensus_stat_path, self.option('catalog_path'),
                                               self.output_dir, self.option('sample_num'))
        self.logger.info("开始进行consensus_stat")
        command = self.add_command("consensus_stat", cmd).run()  # 必须小写，
        self.wait()
        if command.return_code == 0:
            self.logger.info("consensus_stat完成！")
        else:
            self.set_error("consensus_stat出错！", code="35500607")

    def run_snp_depth_density(self):
        """
        """
        cmd = "{} {} -i {} -o {}".format(self.Python_path, self.snp_depth_density_path,
                                         self.option('populations_snps_vcf'), self.output_dir,)
        self.logger.info("开始进行snp_depth_density")
        command = self.add_command("snp_depth_density", cmd).run()  # 必须小写，
        self.wait()
        if command.return_code == 0:
            self.logger.info("snp_depth_density完成！")
        else:
            self.set_error("snp_depth_density出错！", code="35500608")

    def run_tag_dep(self):
        """
        并行投递
        :return:
        """
        tag_dep = os.path.join(self.output_dir, "tag_dep")
        if os.path.exists(tag_dep):
            shutil.rmtree(tag_dep)
        os.mkdir(tag_dep)
        cmd_list = []
        allfiles = os.listdir(self.option('ustacks_path'))
        for i in allfiles:
            if i.endswith(".tags.tsv.gz"):
                input = os.path.join(self.option('ustacks_path'), i)
                sample = i.strip().split(".tags.tsv.gz")
                output = os.path.join(tag_dep, sample[0])
                cmd = "{} {} -i {} -o {}".format(self.Python_all_path, self.tag_dep_path, input, output)
                cmd_list.append(cmd)
        self.run_cmd_more(cmd_list, "tag_dep", 5)

    def run_cmd_more(self, cmd_list, cmd_name, cpu):
        """
        将多个cmd命令并行执行
        """
        cmd_file = os.path.join(self.work_dir, "cmd_list_{}.txt".format(cmd_name))
        wrong_cmd = os.path.join(self.work_dir, "failed_cmd_{}.txt".format(cmd_name))
        with open(cmd_file, "w") as f:
            for cmd in cmd_list:
                f.write(cmd + "\n")
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.para_fly, cmd_file, cpu, wrong_cmd)
        self.run_cmd(cmd_more, "more_" + cmd_name)

    def run_cmd(self, cmd, cmd_name):
        """
        执行cmd
        """
        self.logger.info(cmd)
        command = self.add_command(cmd_name, cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd_name))
        else:
            self.set_error("%s运行失败", code="35500609")

    def run(self):
        super(ConsensusStatTool, self).run()
        self.run_consensus_stat()
        self.run_snp_depth_density()
        self.run_tag_dep()
        self.end()

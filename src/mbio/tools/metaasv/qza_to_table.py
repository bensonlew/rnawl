# -*- coding: utf-8 -*-
#__author__: qingchen.zhang

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir, link_file


class QzaToTableAgent(Agent):
    """
    qiime2 格式转换将qza格式转换为table
    """
    def __init__(self, parent):
        super(QzaToTableAgent, self).__init__(parent)
        options = [
            {"name": "input_file", "type": "infile", "format": "metaasv.qza"}, ##输入的qza文件
            {"name": "type", "type": "string"}, ## 输入的文件类型
            {"name": "prefix", "type": "string"},  ## 文件的前缀
        ]
        self.add_option(options)
        self.step.add_steps('convert')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.convert.start()
        self.step.update()

    def step_end(self):
        self.step.convert.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('input_file').is_set:
            raise OptionError('必须提供输入的fasta序列')

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        file_size = os.path.getsize(self.option("input_file").prop['path']) / (1024 * 1024)
        if file_size < 10:
            memory = 10
        else:
            memory = int(file_size)
        self._memory = '{}G'.format(str(memory))

    def end(self):
        super(QzaToTableAgent, self).end()


class QzaToTableTool(Tool):
    """
    Tool 运行
    """
    def __init__(self, config):
        super(QzaToTableTool, self).__init__(config)
        self.qiime_path = "program/Python/bin/python"
        self.shell = "program/sh"
        self.shell_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/qza2table.sh")
        self.miniconda3 = os.path.join(self.config.SOFTWARE_DIR, "program/miniconda3/bin")
        self.set_environ(PATH=self.miniconda3)
        self.biom_path = "program/Python/bin/"

    def run_qiime2(self):
        """
        输入qiime数据
        :return:
        """
        self.logger.info("开始运行qiime2软件进行降噪！")
        qiime2_env = "qiime2-2020.2" ## 以便以后进行更换版本
        self.infile = self.option("input_file").prop['path']
        cmd = '{} {} {} {} {}'.format(self.shell, self.shell_path, qiime2_env, self.infile, self.output_dir)
        self.logger.info(cmd)
        command = self.add_command('qza2table', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("qiime2转格式成功！")
        else:
            self.set_error("qiime2转格式失败失败！")

    def biom_to_table(self):
        """
        将biom文件转为txt
        :return:
        """
        self.logger.info("开始根据biom文件格式转为TXT文件")
        self.asv_table = os.path.join(self.output_dir, "ASV_table.xls")
        if os.path.exists(self.asv_table):
            os.remove(self.asv_table)
        input_file = os.path.join(self.output_dir, "feature-table.biom")
        cmd = '{}biom convert -i {} -o {} --table-type "OTU table" --to-tsv'.format(self.biom_path, input_file, self.asv_table)
        self.logger.info(cmd)
        command = self.add_command('biom2file', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("biom2file转格式成功！")
        else:
            self.set_error("biom2file转格式失败！")

    def run(self):
        """
        运行
        :return:
        """
        super(QzaToTableTool, self).run()
        if len(os.listdir(self.output_dir)) != 0:
            self.end()
        else:
            if self.option("type") in ['qza']:
                self.run_qiime2()
            elif self.option("type") in ['qza_biom']:
                self.run_qiime2()
                self.biom_to_table()
            self.set_output()
            self.end()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        allfiles = os.listdir(self.output_dir)
        if len(allfiles) != 0:
            self.logger.info("设置结果文件目录成功！")
        else:
            self.set_error("未能正确的生成结果文件！")




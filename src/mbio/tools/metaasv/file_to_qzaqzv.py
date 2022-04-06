# -*- coding: utf-8 -*-
#__author__: qingchen.zhang

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
from mbio.packages.metaasv.common_function import link_dir, link_file


class FileToQzaqzvAgent(Agent):
    """
    qiime2 格式转换
    将fastq、fasta、taxonomy转为qza
    """
    def __init__(self, parent):
        super(FileToQzaqzvAgent, self).__init__(parent)
        options = [
            {"name": "input_file", "type": "string"}, ##输入的fastq文件或者二维表格式的文件
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
        if not self.option('input_file'):
            raise OptionError('必须提供输入的fasta序列')

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        input_size = float(os.path.getsize(self.option("input_file"))) /(1024*1024*100)
        if input_size >= 1:
            memory = 10 + int(10 * input_size)
        else:
            memory = 10
        self._memory = '{}G'.format(str(memory))

    def end(self):
        super(FileToQzaqzvAgent, self).end()


class FileToQzaqzvTool(Tool):
    """
    Tool 运行
    """
    def __init__(self, config):
        super(FileToQzaqzvTool, self).__init__(config)
        self.qiime_path = "miniconda2/bin/python"
        self.shell = "program/sh"
        self.shell_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/file2qza.sh")
        self.fasta_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/fasta2qzaqzv.sh")
        self.taxonomy_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/biom2qzaqzv.sh")
        self.miniconda3 = os.path.join(self.config.SOFTWARE_DIR, "program/miniconda3/bin")
        self.set_environ(PATH=self.miniconda3)
        self.biom_path = "program/Python/bin/"

    def file_to_biom(self):
        """
        将file文件转为biom格式的文件
        :return:
        """
        self.logger.info("开始生成biom文件格式")
        self.asv_table = os.path.join(self.work_dir, self.option("prefix")+".biom")
        cmd = '{}biom convert -i {} -o {} --table-type "OTU table" --to-hdf5'.format(self.biom_path, self.option("input_file"), self.asv_table)
        self.logger.info(cmd)
        command = self.add_command('file2biom', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("file2biom转格式成功！")
        else:
            self.set_error("file2biom转格式失败！")

    def run_fasta(self):
        """
        输入fasta数据
        :return:
        """
        self.logger.info("开始运行qiime2软件进行转格式！")
        qiime2_env = "qiime2-2020.2" ## 以便以后进行更换版本
        input_file = self.option("input_file")
        self.out_file = os.path.join(self.work_dir, self.option("prefix")+".qza")
        self.out_qzv = os.path.join(self.work_dir, self.option("prefix"))
        cmd = '{} {} {} {} {} {}'.format(self.shell, self.fasta_path, qiime2_env, input_file, self.out_file, self.out_qzv) #1
        self.logger.info(cmd)
        command = self.add_command('file2qza', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("qiime2转格式成功！")
        else:
            self.set_error("qiime2转格式失败失败！")

    def run_taxonmy(self):
        """
        输入taxonomy
        :return:
        """
        self.logger.info("开始运行qiime2软件进行转格式！")
        qiime2_env = "qiime2-2020.2" ## 以便以后进行更换版本
        self.out_file = os.path.join(self.work_dir, self.option("prefix")+".qza")
        self.out_qzv = os.path.join(self.work_dir, self.option("prefix"))
        cmd = '{} {} {} {} {} {}'.format(self.shell, self.taxonomy_path, qiime2_env, self.asv_table, self.out_file, self.out_qzv) #1
        self.logger.info(cmd)
        command = self.add_command('biom2qza', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("qiime2转格式成功！")
        else:
            self.set_error("qiime2转格式失败失败！")


    def run(self):
        """
        运行
        :return:
        """
        super(FileToQzaqzvTool, self).run()
        if len(os.listdir(self.output_dir)) != 0:
            self.end()
        else:
            if self.option("type") in ['fasta']:
                self.run_fasta()
            elif self.option("type") in ['table']:
                self.file_to_biom()
                self.run_taxonmy()
            self.set_output()
            self.end()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        if os.path.isdir(self.out_qzv):
            for file in os.listdir(self.out_qzv):
                file_path = os.path.join(self.out_qzv, file)
                link_file(file_path, os.path.join(self.output_dir, self.option("prefix") + ".qzv"))
        else:
            file_path = self.out_qzv + ".qzv"
            if os.path.isfile(file_path):
                link_file(file_path, os.path.join(self.output_dir, self.option("prefix") + ".qzv"))
        prefix_file = os.path.basename(self.out_file)
        out_file = os.path.join(self.output_dir, prefix_file)
        if os.path.exists(out_file):
            os.remove(out_file)
        os.link(self.out_file, out_file)

        self.logger.info("设置结果文件目录成功！")
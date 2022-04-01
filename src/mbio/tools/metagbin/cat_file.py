# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
from biocluster.core.exceptions import OptionError
import subprocess


class CatFileAgent(Agent):
    """
    catreads:通过cat命令将多个序列文件合并为一个
    version 1.0
    author: hgaohao
    last_modify: 2019.01.16
    """

    def __init__(self, parent):
        super(CatFileAgent, self).__init__(parent)
        options = [
            {"name": "fa_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 输入合并序列fa路径
        ]
        self.add_option(options)
        self.step.add_steps('cat_reads')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cat_reads.start()
        self.step.update()

    def step_end(self):
        self.step.cat_reads.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("fa_dir").is_set:
            raise OptionError("请传入fa_dir序列路径")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(CatFileAgent, self).end()

class CatFileTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(CatFileTool, self).__init__(config)
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'

    def cat_seq(self):
        """
        合并序列
        :return:
        """
        file_list=os.listdir(self.option('fa_dir').prop['path'])
        cmd = self.sh_path + 'cat_seq.sh'
        for file in file_list:
            cmd += ' ' + self.option('fa_dir').prop['path'] + '/' + file
        cmd += ' ' + self.work_dir + '/all.scaf.fa'
        self.logger.info('运行cat_seq，将sequence进行合并')
        command = self.add_command("cat_seq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cat_seq合并完成")
        else:
            self.set_error("cat_seq合并失败！")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        if os.path.exists(self.output_dir + '/all.scaf.fa'):
            os.remove(self.output_dir + '/all.scaf.fa')
        os.link(self.work_dir + '/all.scaf.fa',self.output_dir + '/all.scaf.fa')

    def run(self):
        super(CatFileTool, self).run()
        self.cat_seq()
        self.set_output()
        self.end()
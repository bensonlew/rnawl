# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' @20200528
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import shutil
from biocluster.core.exceptions import OptionError
import subprocess


class MergeFileAgent(Agent):
    """
    功能：1. 对reads进行改名操作
    2. 对改名后的reads进行合并
    catreads:通过cat命令将多个序列文件合并为一个
    version 1.0
    """

    def __init__(self, parent):
        super(MergeFileAgent, self).__init__(parent)
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
        super(MergeFileAgent, self).end()

class MergeFileTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(MergeFileTool, self).__init__(config)
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'

    def change_reads_name(self):
        """
        对混拼结果进行reads的重命名
        :return:
        """
        reads_number = 1
        file_dir = self.option("fa_dir").prop["path"]
        self.outfile = os.path.join(self.work_dir, "change_mix")
        if os.path.exists(self.outfile):
            shutil.rmtree(self.outfile)
        os.mkdir(self.outfile)
        for file in os.listdir(file_dir):
            file_path = os.path.join(file_dir, file)
            new_file_path = os.path.join(self.outfile, file)
            with open(file_path, "r") as f, open(new_file_path, "w") as w:
                for line in f:
                    if line[0] == ">":
                        name = "Mix_" + str(reads_number)
                        w.write(">{}\n".format(name))
                        reads_number += 1
                    else:
                        w.write(line)


    def cat_seq(self):
        """
        合并序列
        :return:
        """
        file_list=os.listdir(self.outfile)
        if os.path.exists(os.path.join(self.work_dir, 'all.scaf.fa')):
            os.remove(os.path.join(self.work_dir, 'all.scaf.fa'))
        cmd = self.sh_path + 'cat_seq.sh'
        for file in file_list:
            cmd += ' ' + self.outfile + '/' + file
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
        super(MergeFileTool, self).run()
        self.change_reads_name()
        self.cat_seq()
        self.set_output()
        self.end()
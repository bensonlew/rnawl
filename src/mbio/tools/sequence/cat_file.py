# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import subprocess


class CatFileAgent(Agent):
    """
    将同时样本的多分数据经过cat合并成一个文件，输出文件用样本名称命名
    """

    def __init__(self, parent):
        super(CatFileAgent, self).__init__(parent)
        options = [
            {"name": "input_dir", "type": "infile", "format": "sequence.fastq_dir2"},  # 输入合并序列fq路径
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("input_dir").is_set:
            raise OptionError("请传入fq序列路径")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = int(self.option("input_dir").sups) + 1
        self._memory = str(self._cpu * 5) + 'G'


class CatFileTool(Tool):
    def __init__(self, config):
        super(CatFileTool, self).__init__(config)
        self.list_txt = []

    def link(self, src, dest):
        """
        非不测序列直接链接改名
        """
        dest = os.path.join(self.output_dir, dest)
        if os.path.exists(dest):
            os.remove(dest)
        os.link(src, dest)

    def cat(self, files, dest):
        """
        合并补测的序列
        """
        dest = os.path.join(self.output_dir, dest)
        cmd = '/bin/cat {} > {}'.format(' '.join(files), dest)
        commands = self.add_command('cat_' + str(len(self.cat_cmds)), cmd, shell=True)
        return commands.run()

    def run(self):
        super(CatFileTool, self).run()
        new_list = open(os.path.join(self.output_dir, "list.txt"), 'w')
        self.cat_cmds = []
        for sample, infos in self.option("input_dir").sample_info.items():
            print("{}\t{}".format(sample, infos))
            for direction, files in infos.items():
                if direction == 'l':
                    di = '1'
                elif direction == 'r':
                    di = '2'
                else:
                    di = direction
                dest = '{}.{}.fq.gz'.format(sample, di)
                if len(files) > 1:
                    self.cat_cmds.append(self.cat(files, dest))
                else:
                    self.link(files[0], dest)
                self.list_txt.append([dest, sample, direction])
                new_list.write("{}\t{}\t{}\n".format(dest, sample, direction))
        if self.cat_cmds:
            self.wait()
        self.end()

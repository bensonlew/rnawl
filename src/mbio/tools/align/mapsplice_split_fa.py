# -*- coding: utf-8 -*-

# __author__ = 'linfang.jin'
# time: 2017/1/25 10:54

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError


class MapspliceSplitFaAgent(Agent):
    '''
    mapsplice 第1步：调用fsata.split_single_seq方法劈开ref.fa
    author: linfang.jin
    last_modify: 2017/1/25 10:54
    '''

    def __init__(self, parent):
        super(MapspliceSplitFaAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "infile", 'format': 'sequence.fasta'},
        ]
        self.add_option(options)
        self.step.add_steps('mapsplice_split')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.mapsplice_split.start()
        self.step.update()

    def step_end(self):
        self.step.mapsplice_split.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        :return:
        """
        if not self.option('ref_fa'):
            raise OptionError("必须设置输入的fasta文件")

        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '100G'

    def end(self):
        """
        agent结束后一些文件德操作

        :return:"""
        super(MapspliceSplitFaAgent, self).end()


class MapspliceSplitFaTool(Tool):
    def __init__(self, config):
        super(MapspliceSplitFaTool, self).__init__(config)
        self._version = "v2.1.8"
        self.cmd_path = self.config.SOFTWARE_DIR+"/bioinfo/rna/MapSplice-v2.1.8/"
        self.Python_path = 'program/Python/bin/python '

    def run_mapsplice_split(self):
        """
        运行mapsplice_split
        :return:
        """

        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)

        cmd = "{} {}split_fa.py  -i {}  -o {}".format(self.Python_path, self.cmd_path, self.option('ref_fa').path,
                                                   self.output_dir)
        self.logger.info('运行mapsplce_index')
        command = self.add_command("mapsplce_index_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("mapsplice_split运行完成")
        else:
            self.set_error("mapsplice_split运行出错!")

    def run(self):
        """
        运行
         :return:
        """
        super(MapspliceSplitFaTool, self).run()
        self.run_mapsplice_split()
        self.end()

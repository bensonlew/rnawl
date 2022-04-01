# -*- coding: utf-8 -*-

# __author__ = 'linfang.jin'
# time: 2017/1/25 9:33

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
import re


class MapspliceIndexAgent(Agent):
    '''
    mapsplice 第一步：调用bowtie做index
    author: linfang.jin
    last_modify: 2017/1/25 9:33
    '''

    def __init__(self, parent):
        super(MapspliceIndexAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "infile", 'format': 'sequence.fasta'},
            # {"name": "index_dir", "type": "outfile", "format": 'ref_rna.gene_structure.index_dir',
            #  'default': self.output_dir},
            {"name": "index_prefix", "type": "string", "default": None}
        ]
        self.add_option(options)
        self.step.add_steps('mapsplice_index')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.mapsplice_index.start()
        self.step.update()

    def step_end(self):
        self.step.mapsplice_index.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        :return:
        """
        if not self.option('ref_fa'):
            raise OptionError("必须设置输入的fasta文件")
        # if not self.option('index_dir'):
        #     raise OptionError("必须设置索引文件所在文件夹")
        if not self.option('index_prefix'):
            raise OptionError("必须设置索引文件名称前缀")

        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '100G'

    def end(self):
        """
        agent结束后一些文件德操作

        :return:"""
        super(MapspliceIndexAgent, self).end()


class MapspliceIndexTool(Tool):
    def __init__(self, config):
        super(MapspliceIndexTool, self).__init__(config)
        self._version = "v2.1.8"
        self.cmd_path = "bioinfo/rna/MapSplice-v2.1.8/bin/"

    def run_mapsplce_index(self):
        """
        运行rmats
        :return:
        """
        # if not os.path.exists(self.option('index_dir').path):
        #     os.mkdir(self.option('index_dir').path)
        # self.logger.info(self.option('index_dir').path)
        # self.logger.info(self.option('index_prefix'))
        # index_str = os.path.join(self.option('index_dir').path, self.option('index_prefix'))
        # self.logger.info(index_str)
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        cmd = "{}bowtie-build  {}  {}/{}".format(self.cmd_path, self.option('ref_fa').path, self.output_dir,
                                                 self.option('index_prefix'))
        self.logger.info('运行mapsplce_index')
        command = self.add_command("mapsplce_index_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("mapsplce_index运行完成")
        else:
            self.set_error("mapsplce_index运行出错!")

    def run(self):
        """
        运行
         :return:
        """
        super(MapspliceIndexTool, self).run()
        self.run_mapsplce_index()
        self.end()

# -*- coding: utf-8 -*-


import os, shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import unittest

class Table2biomAgent(Agent):
    def __init__(self, parent):
        super(Table2biomAgent, self).__init__(parent)
        options = [
            {"name": "table", "type": "infile", "format": "tool_lab.simple"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        return True

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(Table2biomAgent, self).end()


class Table2biomTool(Tool):
    def __init__(self, config):
        super(Table2biomTool, self).__init__(config)
        self.src_path =  "/program/Python/bin/biom"


    def run(self):
        """
        运行
        :return:
        """
        super(Table2biomTool, self).run()
        self.logger.info("开始运行命令！")
        self.start_fun()

        self.set_output()

    def start_fun(self):

        infile = self.option("table").path
        outfile = self.output_dir +'/out_json.biom'
        cmd = self.src_path + ' convert  -i {} -o {} --table-type="OTU table" --to-json'.format( infile, outfile)

        self.logger.info(cmd)
        command = self.add_command('convert', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("convert biom succeed")
        else:
            self.set_error("convert biom  failed")
            raise Exception("convert biom  failed")

        outfile = self.output_dir +'/out_hdf5.biom'

        with open(infile) as fr:
            heads = fr.readline().strip().split('\t')
            if 'taxonomy' in heads:
                cmd = self.src_path + ' convert  -i {} -o {} --table-type="OTU table" --to-hdf5 --process-obs-metadata taxonomy'.format( infile, outfile)
            else:
                cmd = self.src_path + ' convert  -i {} -o {} --table-type="OTU table" --to-hdf5'.format( infile, outfile)

        self.logger.info(cmd)
        command2 = self.add_command('convert2', cmd, ignore_error=True).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("convert biom hdf5 succeed")
        else:
            self.set_error("convert biom hdf5 failed")


    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """

        self.end()



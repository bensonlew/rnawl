# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180417

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SelectMatrixAgent(Agent):
    """
    利用select-marker.pl对matrix矩阵进行筛选。
    """
    def __init__(self, parent):
        super(SelectMatrixAgent, self).__init__(parent)
        options = [
            {"name": "depth", "type": "string"},  # "infile","format": ""
            {"name": "mark", "type": "string"},
            {"name": "sample", "type": "string"},       # 样品,分割
            {"name": "xy", "type": "string"},       # 坐标轴，分割
        ]
        self.add_option(options)
        self.step.add_steps('SelectMatrix')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.SelectMatrix.start()
        self.step.update()

    def step_end(self):
        self.step.SelectMatrix.finish()
        self.step.update()

    def check_options(self):
        if not self.option("depth"):
            raise OptionError("请设置depth参数", code="34801801")
        if not self.option("mark"):
            raise OptionError("请设置mark参数", code="34801802")
        if not self.option("xy"):
            raise OptionError("请设置横纵坐标轴位置", code="34801803")
        if not self.option("sample"):
            raise OptionError("请设置sample位置", code="34801804")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '30G'

    def end(self):
        super(SelectMatrixAgent, self).end()
########

class SelectMatrixTool(Tool):
    def __init__(self, config):
        super(SelectMatrixTool, self).__init__(config)
        self.select_path = self.config.PACKAGE_DIR + "/dna_gmap/select-marker.pl"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '

    def Selectmatrix(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{}{} -depth {} -mark {} -y {} -x {} -out {}".format(self.perl_path,
            self.select_path, self.option("depth"), self.option("mark"),
            self.option("xy"), self.option("sample"), self.output_dir + "/select.matrix.xls")		
        self.logger.info(cmd)
        self.logger.info("开始进行SelectMatrix")
        command = self.add_command("selectmatrix", cmd).run()  # SelectMatrix必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("SelectMatrix完成！")
        else:
            self.set_error("SelectMatrix出错！", code="34801801")
            self.set_error("SelectMatrix出错！", code="34801804")

    def run(self):
        super(SelectMatrixTool, self).run()
        self.Selectmatrix()
        self.end()

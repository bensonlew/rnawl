# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180910

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import datetime
import random


class CverrorAgent(Agent):
    """
    群体结构子模块，structure中计算cverror的结果
    """
    def __init__(self, parent):
        super(CverrorAgent, self).__init__(parent)
        options = [
            {"name": "structure_list", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps('script')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.script.start()
        self.step.update()

    def step_end(self):
        self.step.script.finish()
        self.step.update()

    def check_options(self):
        if not self.option("structure_list"):
            raise OptionError("缺少structure_list参数", code="11111")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 3
        self._memory = '10G'

    def end(self):
        super(CverrorAgent, self).end()


class CverrorTool(Tool):
    def __init__(self, config):
        super(CverrorTool, self).__init__(config)
        self.perl_path = 'miniconda2/bin/perl'
        self.script_path = self.config.PACKAGE_DIR + "/dna_evolution/CVerror.pl"

    def cverror(self):
        """
        perl CVerror.pl -i structure.list -o ./
        out_file: best.2.xls cv.error group.list
        :return:
        """
        cmd1 = "{} {} -i {} -o {}" \
            .format(self.perl_path, self.script_path, self.option("structure_list"), self.output_dir)
        self.logger.info(cmd1)
        self.logger.info("开始进行CVerror")
        command = self.add_command("cverror", cmd1).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("CVerror完成！")
        else:
            self.set_error("CVerror出错！")

    def rename_group_name(self):
        """
        这里生成的分组方案名字是数字，但是这样会有问题，所以我们这里前面加上字符串Q
        :return:
        """
        with open(self.output_dir + "/group.list", 'r') as r, open(self.output_dir + "/group.txt", 'w') as w:
            data = r.readlines()
            for line in data:
                temp = line.strip().split('\t')
                w.write("{}\tQ{}\n".format(temp[0], temp[1]))
        os.remove(self.output_dir + "/group.list")
        self.logger.info("删除group.list成功！")
        os.rename(self.output_dir + "/group.txt", self.output_dir + "/group.list")
        self.logger.info("删除group.txt为group.list成功！")

    def run(self):
        super(CverrorTool, self).run()
        self.cverror()
        self.rename_group_name()
        self.end()

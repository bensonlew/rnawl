# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 2018080907

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import datetime
import random
import os
import re


class AdmiStructureAgent(Agent):
    """
    群体结构子模块，计算Structure  --该步一个k值 计算一个tool，看标记个数，几万标记不到半小时，标记越少越快，申请的线程8个
    """
    def __init__(self, parent):
        super(AdmiStructureAgent, self).__init__(parent)
        options = [
            {"name": "pop_bed", "type": "infile", "format": "dna_evolution.bed", "required": True},
            {"name": "k_value", "type": "int", "default": 2},
            {"name": "sample_list", "type": "string"}
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
        if type(self.option("k_value")) != int:
            raise OptionError("k_value一定要是整型数据", code="123345")
        else:
            if self.option("k_value") > 100 or self.option("k_value") < 1:
                raise OptionError("k值要大于%s且小于%s", variables=(1, 100), code='12345')

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 9
        self._memory = '15G'

    def end(self):
        super(AdmiStructureAgent, self).end()


class AdmiStructureTool(Tool):
    def __init__(self, config):
        super(AdmiStructureTool, self).__init__(config)
        self.script_path = self.config.SOFTWARE_DIR + "/bioinfo/dna_evolution/admixture"

    def admistructure_run(self):
        """
        admixture pop.bed 2 --cv -j8 > pop.2.log &&
        paste sample.list pop.2.Q > pop.2.xls
        :return:
        """
        cmd = "{} {} {} --cv -j8 > {}".format(self.script_path, self.option("pop_bed").prop['path'],
                                              self.option("k_value"),
                                              self.output_dir + "/pop.{}.log".format(self.option("k_value")))
        self.logger.info(cmd)
        self.sbatch_task(cmd, "admixture")

    def paste(self):
        cmd = "paste {} {} > {}".format(self.option("sample_list"),
                                        self.work_dir + "/pop.{}.Q".format(self.option("k_value")),
                                        self.output_dir + "/pop.{}.xls".format(self.option("k_value")))
        self.logger.info(cmd)
        code = os.system(cmd)
        if code == 0:
            self.logger.info("执行paste命令成功！")
        else:
            self.set_error("执行paste命令失败！")

    def sbatch_task(self, cmd, cmd_name="script"):
        """
        该模块是为了解决add_command中有> |等定向符号导致的报错，所以这里就写一个外包程序
        :params: cmd 正常的tool中组建的运行shell命令，这里可以包含有定向符，具体执行命令一定要写成绝对路径，或者加载对应的环境变量
        :params: cmd_name 同一个tool中 如果要执行多个cmd要保证cmd_name不唯一
        :return:
        """
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/dna_evolution/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "script_{}.sh".format(now_time)
        self.logger.info("执行脚本路径:{}".format(file_path))
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd)
        code = os.system('/bin/chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改{}为可执行文件失败！".format(file_path))
        shell = "/bioinfo/dna_evolution/script_temp/{}".format(os.path.basename(file_path))
        self.logger.info("开始执行:{}".format(shell))
        command1 = self.add_command(cmd_name, shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("脚本{}执行完成！".format(shell))
        else:
            self.set_error("脚本{}执行出错！".format(shell))
        os.system('rm {}'.format(file_path))

    def set_output(self):
        """
        将结果文件link到output文件夹下面 pop.3.P pop.3.Q -- 暂时不用
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        file1 = "/pop.{}.Q".format(self.option("k_value"))
        file2 = "/pop.{}.P".format(self.option("k_value"))
        if not os.path.exists(self.output_dir + file1):
            os.link(self.work_dir + file1, self.output_dir + file1)
        if not os.path.exists(self.output_dir + file2):
            os.link(self.work_dir + file2, self.output_dir + file2)
        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(AdmiStructureTool, self).run()
        self.admistructure_run()
        self.paste()
        self.end()

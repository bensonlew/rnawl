# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180425

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import datetime
import random
import os


class SoapDenovoAgent(Agent):
    """
    用于数据组装与转基因中的拼接功能--两者共用
    """
    def __init__(self, parent):
        super(SoapDenovoAgent, self).__init__(parent)
        options = [
            {"name": "config_file", "type": "string"}  # A8_10.config
        ]
        self.add_option(options)
        self.step.add_steps('snpeff')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.snpeff.start()
        self.step.update()

    def step_end(self):
        self.step.snpeff.finish()
        self.step.update()

    def check_options(self):
        if not self.option("config_file"):
            raise OptionError("缺少config_file参数", code="34506001")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '230G'
        
    def end(self):
        super(SoapDenovoAgent, self).end()


class SoapDenovoTool(Tool):
    def __init__(self, config):
        super(SoapDenovoTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/SOAPdenovo2-src-r240')
        # self.script_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/SOAPdenovo2-src-r240/SOAPdenovo-127mer'

    def snp_eff(self):
        """
        SOAPdenovo-127mer all -s 02.config/A8-10.config -K 63 -R -o A8-10.denovo -p 20 1> A8-10.ass.log 2> A8-10.ass.err
        :return:
        """
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/WGS/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        sample = os.path.basename(self.option("config_file")).split('.')[0]
        file_path = script_path + "script_{}.sh".format(now_time)
        cmd1 = "SOAPdenovo-127mer all -s {} -K 63 -R -o {} -p 20 1> {} 2> {}"\
            .format(self.option("config_file"), os.path.join(self.output_dir, "{}.denovo".format(sample)),
                    os.path.join(self.work_dir, "{}.ass.log".format(sample)),
                    os.path.join(self.work_dir, "{}.ass.err".format(sample)))
        self.logger.info(cmd1)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd1)
        code = os.system('/bin/chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改%s为可执行文件失败！",variables=(file_path), code="34506001")
        shell = "/bioinfo/WGS/script_temp/{}".format(os.path.basename(file_path))
        self.logger.info("开始进行scrpit")
        command1 = self.add_command("scrpit", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("scrpit运行完成！")
        else:
            self.set_error("scrpit运行出错！", code="34506002")
            self.set_error("scrpit运行出错！", code="34506005")
        os.system('rm {}'.format(file_path))

    def run(self):
        super(SoapDenovoTool, self).run()
        self.snp_eff()
        self.end()

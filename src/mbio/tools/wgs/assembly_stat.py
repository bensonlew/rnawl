# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180427

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import datetime
import random


class AssemblyStatAgent(Agent):
    """
    数据组装模块中对拼接结果进行统计
    """
    def __init__(self, parent):
        super(AssemblyStatAgent, self).__init__(parent)
        options = [
            {"name": "denovo_scafseq", "type": "string"},
            {"name": "denovo_scafstatistics", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps('stat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.stat.start()
        self.step.update()

    def step_end(self):
        self.step.stat.finish()
        self.step.update()
        
    def check_options(self):
        if not self.option("denovo_scafseq"):
            raise OptionError("缺少denovo_scafseq参数", code="34500301")
        if not self.option("denovo_scafstatistics"):
            raise OptionError("缺少denovo_scafstatistics参数", code="34500302")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '10G'
        
    def end(self):
        super(AssemblyStatAgent, self).end()


class AssemblyStatTool(Tool):
    def __init__(self, config):
        super(AssemblyStatTool, self).__init__(config)
        self.scrpit_path = self.config.PACKAGE_DIR + "/wgs/03.denovo.stat.pl"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '
        
    def cnv_diff(self):
        """
        perl 03.denovo.stat.pl -seq Lands.denovo.scafSeq -stat Lands.denovo.scafStatistics -o .
        :return:
        """
        cmd = "{}{} -seq {} -stat {} -o {}".format(self.perl_path, self.scrpit_path, self.option("denovo_scafseq"),
                                                   self.option("denovo_scafstatistics"), self.work_dir)
        self.logger.info(cmd)
        self.logger.info("开始进行scrpit")
        command = self.add_command("scrpit", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("scrpit完成！")
        else:
            self.set_error("scrpit出错！", code="34500301")

    def set_output(self):
        """
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info('开始设置文件夹路径')
        sample_id = os.path.basename(self.option("denovo_scafseq")).split(".")[0]
        time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3] + "_" + str(random.randint(1, 10000))
        results = os.listdir(self.work_dir)
        for f in results:
            if f in ["denovo.stat.xls"]:
                os.link(os.path.join(self.work_dir, f),
                        self.output_dir + "/{}.{}.denovo.stat.xls".format(sample_id, time))
        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(AssemblyStatTool, self).run()
        self.cnv_diff()
        self.set_output()
        self.end()

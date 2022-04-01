# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# last_modify: 2018.08.29

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import time
from bson.objectid import ObjectId


class DrawHeatmapAgent(Agent):
    """
    用R脚本画热图，输入区域的vcf文件（chr3_1_200017536.LD.vcf）
    """
    def __init__(self, parent):
        super(DrawHeatmapAgent, self).__init__(parent)
        options = [
            {"name": "vcf_path", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "area", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps("region_hapmap")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.region_hapmap.start()
        self.step.update()

    def stepfinish(self):
        self.step.region_hapmap.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('vcf_path'):
            raise OptionError('必须输入:vcf_path')
        if not self.option('area'):
            raise OptionError('必须输入:area')

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "50G"

    def end(self):
        super(DrawHeatmapAgent, self).end()


class DrawHeatmapTool(Tool):
    """
    """
    def __init__(self, config):
        super(DrawHeatmapTool, self).__init__(config)
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl'
        self.R_path = "program/R-3.3.3/bin/Rscript"
        self.LDheatmap_vcf2matrix_path = self.config.PACKAGE_DIR + "/dna_evolution/LDheatmap.vcf2matrix.pl"
        self.LDheatmap_path = self.config.PACKAGE_DIR + "/dna_evolution/LDheatmap.R"
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')

    def run_get_data(self):
        """
        获得生成单倍体图的数据
        """
        cmd = "{} {} -i {} -o {}".format(self.perl_path, self.LDheatmap_vcf2matrix_path,
                                         self.option("vcf_path").prop['path'],
                                         os.path.join(self.output_dir, (self.option("area")+".data")))
        self.logger.info('开始运行get_data')
        self.run_cmd(cmd, "get_data")

    def run_heatmap(self):
        """
        生成热图。
        """
        cmd = "{} {} --input {} --output {}".format(self.R_path,
                                                    self.LDheatmap_path,
                                                    os.path.join(self.output_dir, (self.option("area")+".data")),
                                                    os.path.join(self.output_dir, self.option("area")))
        self.logger.info('开始运行heatmap')
        self.run_cmd(cmd, "heatmap")

    def run_cmd(self, cmd, cmd_name):
        """
        执行cmd
        """
        command = self.add_command(cmd_name, cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd_name))
        else:
            self.set_error("{}运行失败".format(cmd_name))
            raise Exception("{}运行失败".format(cmd_name))

    def run(self):
        super(DrawHeatmapTool, self).run()
        self.run_get_data()
        self.run_heatmap()
        self.end()

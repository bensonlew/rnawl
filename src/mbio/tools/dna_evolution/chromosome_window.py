# -*- coding: utf-8 -*-
# __author__ = 'binbin zhao'
# modified 2018.09.06

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import datetime
starttime = datetime.datetime.now()


class ChromosomeWindowAgent(Agent):
    """
    变异位点染色体分布图接口Tool
    """

    def __init__(self, parent):
        super(ChromosomeWindowAgent, self).__init__(parent)
        options = [
            {"name": "pop_table", "type": "infile", "format": "dna_evolution.pop_table"},  # 变异位点比较分析产生的结果文件之一。
            {"name": "step_num", "type": "int", "default": 10},
            {"name": "project_sn", "type": "string"},
            {"name": "task_id", "type": "string"},
            # 用于存放展现的表格是SNP，indel还是ALL。默认为“ALL”
            {"name": "variant_type", "type": "string", "default": "all"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("pop_table"):
            raise OptionError("请pop.table路径不正确", code="34902501")
        if not self.option("step_num") in [5, 10, 50, 100, 200, 500]:
            raise OptionError("step_num 只能是5K，10K，100K，200K，500K", code="34902502")
        if not self.option("variant_type"):
            raise OptionError("请输入variant_type", code="34902503")
        else:
            if not self.option("variant_type") in ["snp", "indel", "all"]:
                raise OptionError("variant_type必须是SNP, INDEL 和 ALL 其中的一种", code="34902504")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(ChromosomeWindowAgent, self).end()


class ChromosomeWindowTool(Tool):
    def __init__(self, config):
        super(ChromosomeWindowTool, self).__init__(config)
        self.python_path = 'miniconda2/bin/python'
        self.chromosome_window_path = self.config.PACKAGE_DIR + "/dna_evolution/step_calc.py"

    def run_chromosome_windows(self):
        """
        运行chromosome_windows

        """
        cmd = "{} {} -p {} -s {} -o {} -v {}".format(
            self.python_path,
            self.chromosome_window_path,
            self.option("pop_table").prop["path"],
            1000*self.option("step_num"),
            self.output_dir,
            self.option("variant_type"))
        command = self.add_command("chromosome_windows", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("chromosome_windows运行完成")
        else:
            self.set_error("chromosome_windows运行失败", code="34902501")

    def set_db(self):
        """
        保存结果到mongo
        """
        if self.option("main_id"):
            api_chromosome = self.api.api("dna_evolution.chromosome_window")
            api_chromosome.add_sg_distribution(self.output_dir, self.option("main_id"), self.option("variant_type"))

    def run(self):
        super(ChromosomeWindowTool, self).run()
        self.run_chromosome_windows()
        self.set_db()
        self.end()

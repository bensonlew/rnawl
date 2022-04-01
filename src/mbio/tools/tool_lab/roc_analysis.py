# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'binbin.zhao'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class RocAnalysisAgent(Agent):
    """
    ROC分析
    """
    def __init__(self, parent):
        super(RocAnalysisAgent, self).__init__(parent)
        options = [
            {"name": "roc_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "main_id", "type": "string"},
            {"name": "smooth", "type": "bool", "default": False},
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("roc_table"):
            raise OptionError("必须输入snp_table文件")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        """
        计算结束
        """
        super(RocAnalysisAgent, self).end()


class RocAnalysisTool(Tool):
    """
    用于相关性分析
    """
    def __init__(self, config):
        super(RocAnalysisTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.R_path = "program/R-3.3.1/bin/Rscript"
        self.script_path = self.config.PACKAGE_DIR + "/tool_lab/roc_calculation.R"

    def calculate_roc(self):
        """
        开始准备输出数据
        """
        self.logger.info(self.script_path)
        self.logger.info("?????????????????????")
        cmd = self.R_path + ' {} -i {} -o {}'.format(
        self.script_path, self.option('roc_table').prop["path"],
        self.output_dir)
        self.logger.info('运行roc程序')
        command = self.add_command('roc', cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("roc分析成功")
        else:
            self.set_error("roc分析失败")
            raise Exception("roc分析失败")

    def run(self):
        """
        运行
        """
        super(RocAnalysisTool, self).run()
        self.calculate_roc()# 对输入数据进行处理，得到new_corr_table.csv
        self.linkfile()
        self.set_db()
        self.end()

    def linkfile(self):
        files = os.listdir(self.output_dir)
        for file in files:
            basename = os.path.basename(file)
            in_link = self.output_dir + '/' + basename
            out_link = self.output_dir + '/' + basename + '.txt'
            if os.path.exists(out_link):
                os.remove(out_link)
            os.link(in_link, out_link)
            if os.path.exists(out_link):
                os.remove(in_link)

    def set_db(self):
        self.logger.info("tool中导表结束")
        api_roc = self.api.api('tool_lab.roc_analysis')
        api_roc.add_sg_roc(self.option('main_id'), self.output_dir, self.option('smooth'))
        self.logger.info("tool中导表结束")
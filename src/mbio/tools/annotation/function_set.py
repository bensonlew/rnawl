# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class FunctionSetAgent(Agent):
    """
    生成功能集的详情表
    """

    def __init__(self, parent):
        super(FunctionSetAgent, self).__init__(parent)
        options = [
            {"name": "database", "type": "string"},  # database
            {"name": "level", "type": "string"},  # 创建功能集的功能水平
            {"name": "member", "type": "string"},# ,号分割
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("database"):
            raise OptionError("database option not set ", code="31204201")
        if not self.option("level"):
            raise OptionError("level option not set ", code="31204202")
        if not self.option("member"):
            raise OptionError("member option not set ", code="31204203")

    def set_resource(self):
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(FunctionSetAgent, self).end()

class FunctionSetTool(Tool):
    def __init__(self, config):
        super(FunctionSetTool, self).__init__(config)
        #self.python_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/python"
        self.python_path = "miniconda2/bin/python"
        self.python_script = self.config.PACKAGE_DIR + '/annotation/function_set_produce.py'

    def run(self):
        """
        运行
        :return:
        """
        super(FunctionSetTool, self).run()
        self.run_funset()
        #self.set_output()
        self.end()

    def run_funset(self):
        #self.out = self.work_dir + '/function_detail.xls'
        member = '"'+self.option('member')+'"'
        cmd = '{} {} -d {} -l {} -m {}  -t {}'.format(self.python_path, self.python_script, self.option("database"),self.option("level"),
                                         member,  self.option('main_id'))
        if self.config.DBVersion:
            cmd += ' -v {}'.format(self.config.DBVersion)
        self.logger.info(cmd)
        command = self.add_command("runfunc", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('runfunc SUCCESS')
        else:
            self.set_error('runfunc Failed', code="31204201")

    def set_output(self):
        pass
        # if os.path.exists(self.output_dir + '/function_detail.xls'):
        #     os.remove(self.output_dir + '/function_detail.xls')
        # os.link(self.out, self.output_dir + '/function_detail.xls')

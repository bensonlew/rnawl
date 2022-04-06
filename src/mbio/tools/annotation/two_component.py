# -*- coding: utf-8 -*-
#20190328
#zouguanqing

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os




class TwoComponentAgent(Agent):
    """
    双组分调控分析
    """

    def __init__(self, parent):
        super(TwoComponentAgent, self).__init__(parent)
        options = [
            {"name": "pfam_anno", "type": "infile", "format": "sequence.profile_table"},   # 基因具体的注释信息
            {"name": "comp_anno", "type": "outfile", "format": "sequence.profile_table"},
            {'name': 'sample_name', "type": "string",'default':'out'}   # 样本名

        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("pfam_anno").is_set:
            raise OptionError("必须设置输入文件")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(TwoComponentAgent, self).end()


class TwoComponentTool(Tool):
    def __init__(self, config):
        super(TwoComponentTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = "/miniconda2/bin/python"
        self.python_script = self.config.PACKAGE_DIR + '/annotation/senser_regulator.py'
        self.database = self.config.SOFTWARE_DIR + '/database/bacgenome/Pfam_two_component/senser_regulator.database'

    def run(self):
        """
        运行
        :return:
        """
        super(TwoComponentTool, self).run()
        self.run_comp_anno()
        self.set_output()
        self.end()

    def run_comp_anno(self):
        table = self.option('pfam_anno').path
        cmd = '{} {} {} {}'.format(self.python_path, self.python_script, self.database, table)
        self.logger.info(cmd)
        command = self.add_command("two_component", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("two_component succeed")
        else:
            self.set_error("two_component failed")


    def set_output(self):
        self.logger.info("set_output")
        xls = self.output_dir + "/{}.senser_regulator.xls".format(self.option('sample_name'))
        if os.path.exists(xls):
            os.remove(xls)
        os.link(self.work_dir + "/senser_regulator.xls",xls)

        stat = self.output_dir + "/{}.senser_regulator.stat".format(self.option('sample_name'))
        if os.path.exists(stat):
            os.remove(stat)
        os.link(self.work_dir + "/senser_regulator.stat",stat)

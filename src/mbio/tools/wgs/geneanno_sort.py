# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modify 20180713

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class GeneannoSortAgent(Agent):
    """
    anno.summary的排序
    """
    def __init__(self, parent):
        super(GeneannoSortAgent, self).__init__(parent)
        options = [
            {"name": "anno_summary", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps('GeneannoSort')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.GeneannoSort.start()
        self.step.update()

    def step_end(self):
        self.step.GeneannoSort.finish()
        self.step.update()

    def check_options(self):
        if not self.option("anno_summary"):
            raise OptionError("请设置anno_summary")    # 必须有

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(GeneannoSortAgent, self).end()


class GeneannoSortTool(Tool):
    def __init__(self, config):
        super(GeneannoSortTool, self).__init__(config)
        self.sort_path = self.config.PACKAGE_DIR + "/wgs/sort_anno.summary.pl"
        self.perl_path = 'miniconda2/bin/perl'

    def geneannosort(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{} {} -i {} -o {}"\
            .format(self.perl_path, self.sort_path, self.option("anno_summary"),
                    self.output_dir + "/anno.summary")
        self.logger.info(cmd)
        self.logger.info("开始进行GeneannoSort")
        command = self.add_command("geneannosort", cmd).run()  # nranno必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("GeneannoSort完成！")
        else:
            self.set_error("GeneannoSort出错！")

    def run(self):
        super(GeneannoSortTool, self).run()
        self.geneannosort()
        self.end()

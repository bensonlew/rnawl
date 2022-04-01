# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os


class FakefakeAgent(Agent):
    """
    基因集IPATH 分析
    last_modify: 2018.3.13
    """
    def __init__(self, parent):
        super(FakefakeAgent, self).__init__(parent)
        options = [
            {"name": "fake", "type": "string", "default": 'feng'},  
        ]
        self.add_option(options)
        self.step.add_steps("fake")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.fake.start()
        self.step.update()

    def stepfinish(self):
        self.step.fake.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '1G'

    def end(self):
        super(FakefakeAgent, self).end()


class FakefakeTool(Tool):
    def __init__(self, config):
        super(FakefakeTool, self).__init__(config)

    def run(self):
        """
        运行
        :return:
        """
        super(FakefakeTool, self).run()
        # self.get_kegg_pics()
        os.system('sleep 9')
        self.logger.info("fake运行完毕")
        self.set_output()
        self.end()

    def set_output(self):
        # all_files = ['gene_fake_input.xls', 'Metabolic_pathways.svg', 'Regulatory_pathways.svg', 'Biosynthesis_of_secondary_metabolities.svg']
        pass


# -*- coding: utf-8 -*-


import os
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.tool_lab.anno_keggp import AnnoKeggp



class AnnoKeggpAgent(Agent):
    """
    代谢组通路注释
    """

    def __init__(self, parent):
        super(AnnoKeggpAgent, self).__init__(parent)
        options = [
            {"name": "metab_name", "type": "string", "default" : ""},
            {"name": "detail_out", "type": "outfile", "format": "sequence.profile_table"},  # 注释层级结果
            {"name": "stat_out", "type": "outfile", "format": "sequence.profile_table"},  # 注释统计结果
        ]
        self.add_option(options)
        self.step.add_steps("keggp")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.keggp.start()
        self.step.update()

    def stepfinish(self):
        self.step.keggp.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('metab_name'):
            raise OptionError('必须输入代谢物详情')
        return True

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):

        super(AnnoKeggpAgent, self).end()


class AnnoKeggpTool(Tool):
    def __init__(self, config):
        super(AnnoKeggpTool, self).__init__(config)

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoKeggpTool, self).run()
        obj = AnnoKeggp()
        #self.logger.info(self.option("metab_name").prop["path"])
        try:
            self.detial_out = self.output_dir+'/detail.xls'
            self.stat_out = self.output_dir+'/stat.xls'
            obj.run(self.option('metab_name'), self.detial_out, self.stat_out)
        except Exception as e:
            self.logger.error(e)
            self.set_error(e)
        self.set_output()
        self.end()



    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.option('detail_out').set_path(self.detial_out)
        self.option('stat_out').set_path(self.stat_out)
        self.logger.info("设置anno_keggc分析结果目录成功")

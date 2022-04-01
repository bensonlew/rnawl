# -*- coding: utf-8 -*-
# __author__ = "shenghe"
# last_modify:20160809

import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.denovo_rna.express.go_gene_regulate import GO_level_2_regulate


class GoRegulateAgent(Agent):
    """
    version v1.0
    author: hesheng
    last_modify: 2016.08.10
    """
    def __init__(self, parent):
        super(GoRegulateAgent, self).__init__(parent)
        options = [
            {"name": "diff_stat", "type": "infile", "format": "rna.diff_stat_table"},
            {"name": "go_level_2", "type": "infile", "format": "annotation.go.level2"}
        ]
        self.add_option(options)
        self.step.add_steps("goregulate")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.goregulate.start()
        self.step.update()

    def stepfinish(self):
        self.step.goregulate.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        for item in self._options.values():
            if not item.is_set:
                raise OptionError('缺少输入文件')


    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = ''

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["GO_regulate.xls", "xls", "基因上下调在GO2level层级分布情况表"],
        ])
        super(GoRegulateAgent, self).end()


class GoRegulateTool(Tool):
    """
    """
    def __init__(self, config):
        super(GoRegulateTool, self).__init__(config)

    def run_go_regulate(self):
        try:
            self.logger.info(self.output_dir + '/GO_regulate.xls')
            GO_level_2_regulate(self.option('diff_stat').path, self.option('go_level_2').path, self.output_dir + '/GO_regulate.xls')
            self.end()
        except Exception:
            self.set_error('计算错误')

    def run(self):
        super(GoRegulateTool, self).run()
        self.run_go_regulate()

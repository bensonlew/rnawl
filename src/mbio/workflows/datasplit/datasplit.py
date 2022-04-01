# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""数据拆分"""

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class DatasplitWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DatasplitWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'board_path', 'type': "string"}  # 下机数据路径
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.bcl2fastq = self.add_tool("datasplit.test_workflow")

    def check_options(self):
        """
        检查参数设置
        """
        return True

    def run_bcl2fastq(self):
        self.logger.info("开始运行bcl2fastq")
        self.bcl2fastq.set_options({
            "sample_info": self.option('board_path')
        })
        self.bcl2fastq.run()


    def run(self):
        self.run_bcl2fastq()
        self.end()
        super(DatasplitWorkflow, self).run()

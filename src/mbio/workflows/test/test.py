# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

"""测试childend的问题"""

from biocluster.workflow import Workflow


class TestWorkflow(Workflow):
    """
    报告中调用组间差异性分析检验时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TestWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "test", "type": "string", "default": 'qiuqiu'}

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tools = []

    def run_multiple(self):
        for i in range(16):
            self.test = self.add_tool("test.test_qiu")
            self.test.set_options({'test': self.option('test')})
            self.test.run()
            self.tools.append(self.test)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(TestWorkflow, self).end()

    def set_out(self):
        print "tools all is end, start workflow end"
        self.end()

    def run(self):
        self.run_multiple()
        self.on_rely(self.tools, self.set_out)
        super(TestWorkflow, self).run()

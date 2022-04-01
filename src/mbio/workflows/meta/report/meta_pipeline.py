# -*- coding: utf-8 -*-
# __author__ = 'hongdong.xuan'
from biocluster.workflow import Workflow


class MetaPipelineWorkflow(Workflow):
    """
    用于meta一键化交互的计算
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetaPipelineWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "data", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "pipe_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.pipe_submit_all = self.add_tool('meta.pipe_submit')

    def run_all(self):
        """
        运行子分析接口投递的tool
        :return:
        """
        options = {
            'data': self.option("data"),
            'pipe_id': self.option("pipe_id"),
            'task_id': self.get_task_id(),
            'batch_task_id': self._sheet.id
        }
        self.pipe_submit_all.set_options(options)
        self.pipe_submit_all.on('end', self.end)
        self.pipe_submit_all.run()

    def get_task_id(self):
        split_id = self._sheet.id.split('_')
        split_id.pop()
        split_id.pop()
        task_id = '_'.join(split_id)
        return task_id

    def run(self):
        self.run_all()
        super(MetaPipelineWorkflow, self).run()

    def end(self):
        super(MetaPipelineWorkflow, self).end()

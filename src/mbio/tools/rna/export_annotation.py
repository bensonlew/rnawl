# -*- coding: utf-8 -*-
# __author__ = 'shijin'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os

class ExportAnnotationAgent(Agent):
    """
    将数据库中存在的注释文件导入到注释结果目录文件夹中
    version v1.0
    author: shijin
    last_modify: 2017.7.14
    """
    def __init__(self, parent):
        super(ExportAnnotationAgent, self).__init__(parent)
        options = [
            {"name": "anno_path", "type": "string", "default": ""},  # 参考基因组数据库
            {"name": "target_path", "type": "string", "default": ""}
        ]
        self.add_option(options)
        self.step.add_steps("export_annotation")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.export_annotation.start()
        self.step.update()

    def stepfinish(self):
        self.step.export_annotation.finish()
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
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(ExportAnnotationAgent, self).end()


class ExportAnnotationTool(Tool):
    """
    将数据库中存在的注释文件导入到注释结果目录文件夹中
    """
    def __init__(self, config):
        super(ExportAnnotationTool, self).__init__(config)
        self._version = '1.0.1'

    def start_move(self):
        os.system("cp -r {} {}".format(self.option("anno_path"), self.option("target_path")))

    def run(self):
        super(ExportAnnotationTool, self).run()
        self.start_move()
        self.end()

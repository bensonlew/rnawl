# -*- coding: utf-8 -*-
# __author__ = 'shijin'

from __future__ import division
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from bson import ObjectId
import shutil
import os

class GetSampleAgent(Agent):
    """
    通过样本集id获取到样本文件夹
    """
    def __init__(self, parent):
        super(GetSampleAgent, self).__init__(parent)
        options = [
            {"name": "samplebase_id", "type": "string", "default": ""},
        ]
        self.add_option(options)
        self.step.add_steps("get_sample")
        self.on('start', self.start_get_sample)
        self.on("end", self.end_get_sample)

    def start_get_sample(self):
        self.step.get_sample.start()
        self.step.update()

    def end_get_sample(self):
        self.step.get_sample.finish()
        self.step.update()

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 4
        self._memory = "4G"

class GetSampleTool(Tool):
    def __init__(self, config):
        super(GetSampleTool, self).__init__(config)
        self.uri = config.uri
        self.database = self.uri['samplebase']
        self.sample_id_list = list()
        self.sample_name_dict = dict()
        self.sample_path_dict = dict()

    def run(self):
        super(GetSampleTool, self).run()
        self.get_sample()
        self.get_new_name()
        self.get_sample_path()
        self.move_sample_to_workdir()
        self.end()

    def get_sample(self):  # 获取到sample的id号
        col = self.database["sg_test_batch_specimen"]
        results = col.find({"test_batch_id": ObjectId(self.option("samplebase_id"))})
        for result in results:
            self.sample_id_list.append(result["test_specimen_id"])

    def get_new_name(self):  # 获取到sample id和改名后名字的对应关系
        col = self.database["sg_test_batch_task_specimen"]
        results = col.find({"test_batch_id": ObjectId(self.option("samplebase_id"))})
        for result in results:
            self.sample_name_dict[result["test_specimen_id"]] = result["alias_name"]

    def get_sample_path(self):  # 获取到sample在本地存放的路径
        col = self.database["sg_test_specimen"]
        for sample_id in self.sample_id_list:
            result = col.find_one({"_id": sample_id})  # _id为sample所带的固有id
            file_path = result["file_path"]
            self.sample_path_dict[sample_id] = file_path

    def move_sample_to_workdir(self):  # 将sample移动到工作目录下
        for sample_id in self.sample_id_list:
            target_path = os.path.join(self.output_dir, self.sample_name_dict[sample_id])
            if not os.path.exists(target_path):
                shutil.copy(self.sample_path_dict[sample_id], target_path)
            else:
                source = self.sample_path_dict[sample_id]
                target = target_path
                self.cat_together(source, target)

    @staticmethod
    def cat_together(source, target):
        with open(target, "a") as a:
            file = open(source, "r")
            for line in file:
                a.write(line)
            file.close()

# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
import os, re, subprocess, shutil, time
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class ModelAgent(Agent):
    """
    先将输入的fq进行整合，然后用bowtie2进行比对，最后计算丰度
    author: guhaidong
    last_modified: 20180521
    """

    def __init__(self, parent):
        super(ModelAgent, self).__init__(parent)
        options = [
            {"name": "train_file", "type": "infile", "format": "sequence.profile_table"},
            {"name": "model", "type": "string"},
            {"name": "model_type", "type": "string"},
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)
        self.step.add_steps("model")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.model.start()
        self.step.update()

    def stepfinish(self):
        self.step.model.finish()
        self.step.update()

    def check_option(self):
        if not self.option('train_file').is_set:
            raise OptionError('必须输入预测文件')
        if not self.option('model').is_set:
            raise OptionError('必须输入疾病类型')
        if not self.option('model_type').is_set:
            raise OptionError('必须输入模型方法')
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"


class ModelTool(Tool):
    def __init__(self, config):
        super(ModelTool, self).__init__(config)
        if self.option('model_type') == 'randomforest':
            tmp_model_type = 'rf'
        else:
            tmp_model_type = self.option('model_type')
        self.file_name = self.option('model') + '_' + tmp_model_type + '_model'
        self.model_path = self.config.SOFTWARE_DIR + '/bioinfo/model/' + self.file_name
        self.script = self.config.PACKAGE_DIR + '/hmdb/scripts/train.py'
        self.python_path = '/miniconda2/bin/python'

    def run_model(self):
        self.logger.info("运行模型预测")
        cmd = "%s %s -m %s -i %s -t txt -model %s -mode class -o %s -disease %s" % (
        self.python_path, self.script, self.option('model_type'), self.option('train_file').path, self.model_path,
        self.output_dir + '/model_result.xls', self.option('model'))
        self.logger.info("model command: %s" % cmd)
        command = self.add_command("model", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("模型预测完成")
        else:
            self.logger.info("模型预测出错")

    def set_output(self):
        self.logger.info("设置结果目录")
        self.option("out").set_path(self.output_dir + '/model_result.xls')
        self.logger.info("设置结果目录成功")

    def run(self):
        super(ModelTool, self).run()
        self.run_model()
        self.set_output()
        self.end()

# -*- coding: utf-8 -*-
# __author__ :zhouxuan

import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import re

class MergeFastqAgent(Agent):
    """
    医学流程数据拆分部分fastq.gz合并
    author：xuan.zhou
    last_modify by hongdong
    last_modify: 2017.1020
    """

    def __init__(self, parent):
        super(MergeFastqAgent, self).__init__(parent)
        options = [
            {"name": "sample_dir_name", "type": "string"},
            {"name": "data_dir", "type": "infile", "format": "paternity_test.data_dir"},
            {"name": "result_dir", "type": "string"},
            {"name": "ws_single", "type": "string", "default": "false"},
        ]
        self.add_option(options)
        self.step.add_steps("merge_fastq")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self.r1_path = ''
        self.r2_path = ''
        self.new_name = ''

    def stepstart(self):
        self.step.merge_fastq.start()
        self.step.update()

    def stepfinish(self):
        self.step.merge_fastq.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('sample_dir_name'):
            raise OptionError("必须输入指定样本文件夹名称")
        if not self.option('data_dir'):
            raise OptionError("必须提供样本序列文件夹")
        return True

    def set_resource(self):  # 后续需要测试确认
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 5
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            ["", "", ""],
        ])
        super(MergeFastqAgent, self).end()

class MergeFastqTool(Tool):
    def __init__(self, config):
        super(MergeFastqTool, self).__init__(config)
        self._version = "v1.0"

    def run(self):
        """
        运行
        :return:
        """
        super(MergeFastqTool, self).run()
        self.run_mf()
        self.set_output()
        self.end()

    def run_mf(self):
        file_path = os.path.join(self.option("data_dir").prop['path'], self.option('sample_dir_name'))
        file_name = os.listdir(file_path)
        r1_list = []
        r2_list = []
        for p in file_name:
            m = re.match('(.*)_R1_([0-9].*).fastq.gz', p)
            if m:
                r1_list.append(p)
            else:
                r2_list.append(p)
        r1_list.sort()
        r2_list.sort()
        self.logger.info('r1_list:{}'.format(r1_list))
        self.logger.info('r2_list:{}'.format(r2_list))
        sample_name_ = self.option('sample_dir_name').split("_")
        self.new_name = ("_").join(sample_name_[1:])
        self.logger.info(self.new_name)
        self.r1_path = os.path.join(file_path, self.new_name + "_R1.fastq")
        self.r2_path = os.path.join(file_path, self.new_name + "_R2.fastq")
        r1_path = []
        for q in r1_list:
            new_file_path = os.path.join(file_path, q)
            r1_path.append(new_file_path)
        result1_1 = os.system('cat {} {} {} {} > {}'.format(r1_path[0], r1_path[1], r1_path[2], r1_path[3], self.r1_path + ".gz"))
        self.logger.info('result1_1:{}'.format(result1_1))
        # result1_2 = os.system('gzip {}'.format(self.r1_path))
        # self.logger.info('result1_2:{}'.format(result1_2))
        if self.option("ws_single") == 'false':
            r2_path = []
            for l in r2_list:
                new_file_path_2 = os.path.join(file_path, l)
                r2_path.append(new_file_path_2)
            result2_1 = os.system('cat {} {} {} {} > {}'.format(r2_path[0], r2_path[1], r2_path[2], r2_path[3], self.r2_path + ".gz"))
            self.logger.info('result2_1:{}'.format(result2_1))
            # result2_2 = os.system('gzip {}'.format(self.r2_path))
            # self.logger.info('result2_2:{}'.format(result2_2))

    def set_output(self):
        """
        把合并成功的所有双端fastq序列放在结果文件夹中
        :return:
        """
        gz_r1_file_path = self.r1_path + ".gz"
        self.logger.info(gz_r1_file_path)
        gz_r2_file_path = self.r2_path + ".gz"
        self.logger.info(gz_r2_file_path)
        if os.path.exists(gz_r1_file_path):
            # os.link(gz_r1_file_path, self.output_dir + '/' + self.new_name + "_R1.fastq.gz")
            os.link(gz_r1_file_path, self.option("result_dir") + '/' + self.new_name + "_R1.fastq.gz")
        else:
            self.set_error("no {}_R1_fastq.gz file".format(self.new_name))
            raise Exception("no {}_R1_fastq.gz file".format(self.new_name))
        if self.option("ws_single") == 'false':
            if os.path.exists(gz_r2_file_path):
                # os.link(gz_r2_file_path, self.output_dir + '/' + self.new_name + "_R2.fastq.gz")
                os.link(gz_r2_file_path, self.option("result_dir") + '/' + self.new_name + "_R2.fastq.gz")
            else:
                self.set_error("no {}_R2_fastq.gz file".format(self.new_name))
                raise Exception("no {}_R2_fastq.gz file".format(self.new_name))
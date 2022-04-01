#-*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
import json
from biocluster.core.exceptions import OptionError

class FastpStatAgent(Agent):
    """
    对fastp软件计算结果进行统计
    """
    def __init__(self, parent):
        super(FastpStatAgent, self).__init__(parent)
        options = [
            {"name": "json_dir", "type": "infile", "format": "meta.qc.json_dir"},  # fastp软件计算json结果文件夹
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("json_dir").is_set:
            raise OptionError("必须提供fastp软件计算json结果文件夹")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "2G"  # 内存10G增加到15G by GHD @20180427

    def end(self):
        super(FastpStatAgent, self).end()


class FastpStatTool(Tool):
    def __init__(self, config):
        super(FastpStatTool, self).__init__(config)

    def stat_info(self):
        """
        统计结果
        :return:
        """
        json_dir = self.option('json_dir').prop['path']
        raw_stat = self.work_dir + "/reads.rawData.stat"
        clean_stat = self.work_dir + "/reads.cleanData.stat"
        with open(raw_stat, 'w') as w, open(clean_stat, 'w') as m:
            w.write('#Sample\tReadsNum\tBasesNum\tAverageLength\n')
            m.write('#Sample\tReadsNum\tBasesNum\tAverageLength\n')
            for sample_file in os.listdir(json_dir):
                json_path = os.path.join(json_dir, sample_file)
                json_dict = json.loads(open(json_path, "r").read())
                summary = json_dict["summary"]
                specimen_name = sample_file.split('.')[0]
                raw_average = (int(summary["before_filtering"]["read1_mean_length"]) + int(summary["before_filtering"]["read2_mean_length"])) / 2
                raw_read_num = int(summary["before_filtering"]["total_reads"])
                raw_base = int(summary["before_filtering"]["total_bases"])
                clean_average = (int(summary["after_filtering"]["read1_mean_length"]) + int(summary["after_filtering"]["read2_mean_length"])) / 2
                clean_read_num = int(summary["after_filtering"]["total_reads"])
                clean_base = int(summary["after_filtering"]["total_bases"])
                w.write('{}\t{}\t{}\t{}\n'.format(specimen_name, raw_read_num, raw_base, raw_average))
                m.write('{}\t{}\t{}\t{}\n'.format(specimen_name, clean_read_num, clean_base, clean_average))

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        self.logger.info('正在生成文件结果目录')
        if os.path.exists(self.output_dir + "/reads.rawData.stat"):
            os.remove(self.output_dir + "/reads.rawData.stat")
        os.link(self.work_dir + "/reads.rawData.stat", self.output_dir + "/reads.rawData.stat")
        if os.path.exists(self.output_dir + "/reads.cleanData.stat"):
            os.remove(self.output_dir + "/reads.cleanData.stat")
        os.link(self.work_dir + "/reads.cleanData.stat", self.output_dir + "/reads.cleanData.stat")
        self.logger.info('生成文件结果目录成功！')

    def run(self):
        super(FastpStatTool, self).run()
        self.stat_info()
        self.set_output()
        self.end()


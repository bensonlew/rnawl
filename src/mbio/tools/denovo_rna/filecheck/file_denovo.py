# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.sequence.fastq import FastqFile
from mbio.files.sequence.file_sample import FileSampleFile
import os
import re


class FileDenovoAgent(Agent):
    """
    version 1.0
    author: qindanhua
    last_modify: 2016.06.22
    用于在denovoRNA的workflow开始之前对所输入的文件进行详细的内容检测
    """

    def __init__(self, parent):
        super(FileDenovoAgent, self).__init__(parent)
        options = [
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # fastq文件夹
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 有生物学重复的时候的分组文件
            {"name": "control_file", "type": "infile", "format": "denovo_rna.express.control_table"}  # 对照组文件，格式同分组文件
        ]
        self.add_option(options)
        self.step.add_steps("file_check")
        self.on('start', self.start_file_check)
        self.on('end', self.end_file_check)

    def start_file_check(self):
        self.step.file_check.start()
        self.step.update()

    def end_file_check(self):
        self.step.file_check.finish()
        self.step.update()

    def check_option(self):
        if not self.option('fastq_dir'):
            raise OptionError("必须输入fastq文件参数")
        if not self.option('fastq_dir').prop['has_list_file']:
            raise OptionError('fastq文件夹中必须含有一个名为list.txt的文件名--样本名的对应文件')
        if not self.option('control_file').is_set:
            raise OptionError("必须输入对照组文件，用于查找上下调基因")
        if not self.option('fq_type').is_set:
            raise OptionError("必须设置测序类型：PE or SE")
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError("测试类型只能是PE或者SE")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 10
        self._memory = ''


class FileDenovoTool(Tool):
    """
    检查denovo rna输入文件的格式是否符合要求
    """
    def __init__(self, config):
        super(FileDenovoTool, self).__init__(config)
        self.samples = list()

    def check_fastq(self):
        self.logger.info("正在检测fastq_dir文件")
        if not self.option("fastq_dir").prop["has_list_file"]:
            raise OptionError('fastq文件夹中必须含有一个名为list.txt的文件名--样本名的对应文件')
        self.samples = self.option("fastq_dir").prop["samples"]
        col_num = self.get_list_info()
        if self.option("fq_type") in ["PE"] and col_num != 3:
            raise OptionError("PE文件夹的list应该包含三行信息，文件名-样本名-左端OR右端")
        file_list = FileSampleFile()
        list_txt = os.path.join(self.option('fastq_dir').prop['path'], "list.txt")
        file_list.set_path(list_txt)
        file_sample = file_list.get_list()
        # self.logger.info('%s' % file_sample)
        if self.option('fq_type') == 'PE':
            for i in file_sample.keys():
                if len(file_sample[i]) != 2:
                    raise OptionError("PE测序时，每个样本至少有一个左端fq和右端fq文件")
        files = self.option('fastq_dir').prop['fastq_basename']
        # self.logger.info('%s' % files)
        for f in files:
            fq_path = os.path.join(self.option('fastq_dir').prop['path'], f)
            my_fastq = FastqFile()
            my_fastq.set_path(fq_path)
            if re.search('\.gz$', f) or re.search('\.gzip$', f):
                my_fastq.is_gz = True
            my_fastq.check_content()
        self.logger.info("fastq文件检测完毕")

    def get_list_info(self):
        list_path = self.option("fastq_dir").prop["path"] + "/list.txt"
        with open(list_path, "r") as l:
            col_num = len(l.readline().strip().split())
        return col_num

    def check_group(self):
        if self.option('group_table').is_set:
            self.logger.info("正在检测group文件")
            self.option("group_table").get_info()
            gp_sample = self.option("group_table").prop["sample"]
            for gp in gp_sample:
                if gp not in self.samples:
                    raise Exception("group表出错, 样本{}在fastq文件中未出现".format(gp))
        else:
            self.logger.info("未检测到group文件， 跳过...")
        self.logger.info("group文件检测完毕")

    def check_control(self):
        self.logger.info("正在检测control文件")
        vs_list = self.option("control_file").prop["vs_list"]
        con_samples = []
        if self.option('group_table').is_set:
            group_scheme = self.option('group_table').prop['group_scheme'][0]
            group_name = self.option('group_table').get_group_name(group_scheme)
        for vs in vs_list:
            for i in vs:
                if i not in con_samples:
                    con_samples.append(i)
            for cp in con_samples:
                if self.option('group_table').is_set:
                    if cp not in group_name:
                        raise Exception("control表出错，分组{}在fastq文件中未出现".format(cp))
                else:
                    if cp not in self.samples:
                        raise Exception("control表出错，样本{}在fastq文件中未出现".format(cp))
        self.logger.info("control文件检测完毕")

    def run(self):
        super(FileDenovoTool, self).run()
        self.check_fastq()
        self.check_group()
        self.check_control()
        self.end()

# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""解压bcl2fastq的结果"""
import os
import re
import errno
import gzip
import time
import multiprocessing
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class UnzipAgent(Agent):
    def __init__(self, parent=None):
        super(UnzipAgent, self).__init__(parent)
        self._run_mode = "ssh1"
        options = [
            {'name': 'sample_info', 'type': "infile", 'format': 'datasplit.miseq_split'},  # 样本拆分信息表
            {'name': 'data_path', 'type': "string"}  # bcl2fastq的下机输出目录
        ]
        self.add_option(options)

    def check_option(self):
        if not self.option('sample_info').is_set:
            raise OptionError("参数sample_info不能为空")
        if not self.option('data_path'):
            raise OptionError("参数data_path不能为空")
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = ''


class UnzipTool(Tool):
    def __init__(self, config):
        super(UnzipTool, self).__init__(config)
        self._version = 1.0
        self.fastqs = list()

    def find_parent_sample_id(self, filename):
        """
        用一个文件名查找一个样本的id
        :param filename: 文件名
        """
        for p in self.option('sample_info').prop["parent_sample"]:
            sample_id = re.sub(r'-', r'_', p["sample_id"])
            if re.search(sample_id, filename):
                return p["sample_id"]
            if re.search(p["sample_id"], filename):
                return p["sample_id"]
        raise Exception("没有找到对应的样本: {}".format(filename))

    def make_ess_dir(self):
        unzip_dir = os.path.join(self.work_dir, "unzip")
        dir_list = [unzip_dir]
        for name in dir_list:
            try:
                os.makedirs(name)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(name):
                    pass
                else:
                    raise OSError("创建目录失败")

    def unzip(self):
        """
        解压和重命名bcl2fastq的结果，供后续的fastxtoolkit使用
        """
        p = list()
        for library_type in self.option('sample_info').prop["library_type"]:
            fq_dir = os.path.join(self.option('data_path'), library_type)
            fq_list = os.listdir(fq_dir)
            for fq in fq_list:
                sample_id = self.find_parent_sample_id(fq)
                if re.search(r'.+_R1_001.+', fq):
                    unzip_name = sample_id + "_r1.fastq"
                elif re.search(r'.+_R2_001.+', fq):
                    unzip_name = sample_id + "_r2.fastq"
                else:
                    raise Exception("错误的文件名")
                fq = os.path.join(fq_dir, fq)
                self.logger.debug("开始解压文件" + fq)
                unzip_name = os.path.join(self.work_dir, "unzip", unzip_name)
                self.fastqs.append(unzip_name)
                t = multiprocessing.Process(target=self.unzip_file, args=(fq, unzip_name))
                p.append(t)
        for my_p in p:
            my_p.daemon = True
            my_p.start()
            time.sleep(2)
        for my_p in p:
            my_p.join()

    def unzip_file(self, infile, outfile):
        """
        输入一个压缩文件名和一个输出文件名，进行文件的解压
        """
        with gzip.open(infile, 'rb') as r, open(outfile, 'wb') as w:
            for line in r:
                w.write(line)

    def run(self):
        super(UnzipTool, self).run()
        self.make_ess_dir()
        self.unzip()
        self.end()

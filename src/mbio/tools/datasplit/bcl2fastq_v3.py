# -*_ coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171218
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class Bcl2fastqV3Agent(Agent):
    """
    bcl2fastq2-v2.17.1.14
    文库拆分
    """
    def __init__(self, parent=None):
        super(Bcl2fastqV3Agent, self).__init__(parent)
        options = [
            {"name": "data_path", "type": "string"},  # 下机数据路径
            {"name": "sample_sheet", "type": "infile", "format": "datasplit.sample_sheet"},  # csv文件
            {"name": "bases_mask", "type": "string"},  # 测序模式，hiseq4000 为y151,i6nn,y151; miseq为y301,i6,y301
            {"name": "barcode_mismatch", "type": "int", "default": 0},  # --barcode-mismatches，barcode错配数
        ]
        self.add_option(options)

    def check_options(self):
        """参数检查"""
        if not self.option("data_path"):
            raise OptionError("缺少下机数据路径")
        # base_path = self.option("data_path") + "/Data/Intensities/BaseCalls/"
        # if not os.path.exists(base_path):
        #     raise OptionError("下机数据里缺少/Data/Intensities/BaseCalls/文件夹")
        if not self.option("sample_sheet").is_set:
            raise OptionError("缺少csv表")
        if not self.option("bases_mask"):
            raise OptionError("缺少测序模式")
        if self.option("barcode_mismatch") < 0:
            raise OptionError("参数barcode_mismatch:{}不能小于0，请检查".format(self.option("barcode_mismatch")))
        if re.match(r"\d+", self.option("bases_mask")):
            raise OptionError("测序模式有问题，请检查")

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 5
        self._memory = "30G"


class Bcl2fastqV3Tool(Tool):
    def __init__(self, config):
        super(Bcl2fastqV3Tool, self).__init__(config)
        self._version = 2.0
        # self.bcl2fastq_path = "bioinfo/seq/bcl2fastq2-v2.17.1.14/bin/bcl2fastq"
        self.bcl2fastq_path = "bioinfo/seq/bcl2fastq2-v2.20/bin/bcl2fastq"
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/bin')

    def bcl2fastq(self):
        """
        运行bcl2fastq
        """
        self.bcl_path = self.option("data_path") + "/Data/Intensities/BaseCalls/"
        # cmd = self.bcl2fastq_path + ' -i ' + self.bcl_path + ' -o ' + self.output_dir + ' --sample-sheet '\
        #       + self.option("sample_sheet").prop["path"] + ' --use-bases-mask ' + self.option("bases_mask")\
        #       + ' --ignore-missing-bcl ' + '-R ' + self.option("data_path") + ' -r 4 -w 4 -d 2 -p 10 --barcode-mismatches 0'
        cmd = self.bcl2fastq_path + ' -i ' + self.bcl_path + ' -o ' + self.output_dir + ' --sample-sheet '\
              + self.option("sample_sheet").prop["path"] + ' --use-bases-mask ' + self.option("bases_mask")\
              + ' --ignore-missing-bcl ' + '-R ' + self.option("data_path") + ' -r 4 -w 4 -p 10 '\
              + '--barcode-mismatches ' + str(self.option("barcode_mismatch"))
        # command = self.add_command('bcl2fastq', cmd).run()
        command = self.add_command('bcl2fastq', cmd, script_dir=False, default_return_code=0, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行bcl2fastq完成")
        else:
            self.set_error("运行bcl2fastq出错,可能是测序模式有问题，请检查")
            raise Exception("运行bcl2fastq出错,可能是测序模式有问题，请检查")

    def run(self):
        super(Bcl2fastqV3Tool, self).run()
        self.bcl2fastq()
        self.end()

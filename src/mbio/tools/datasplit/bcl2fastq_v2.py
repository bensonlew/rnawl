# -*_ coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171116
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class Bcl2fastqV2Agent(Agent):
    """
    bcl2fastq2-v2.17.1.14
    文库拆分
    """
    def __init__(self, parent=None):
        super(Bcl2fastqV2Agent, self).__init__(parent)
        options = [
            {"name": "sample_info", "type": "infile", "format": "datasplit.sample_info"},  # 下机数据板信息
        ]
        self.add_option(options)

    def check_options(self):
        """参数检查"""
        if not self.option("sample_info").is_set:
            raise OptionError("缺少下机数据的信息")

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 5
        self._memory = "30G"


class Bcl2fastqV2Tool(Tool):
    """
    """
    def __init__(self, config):
        super(Bcl2fastqV2Tool, self).__init__(config)
        self._version = 2.0
        self.bcl2fastq_path = "bioinfo/seq/bcl2fastq2-v2.17.1.14/bin/bcl2fastq"
        self.option('sample_info').get_info()
        self.data_path = self.option('sample_info').prop["data_path"]
        self.run_info = os.path.join(self.data_path, "RunInfo.xml")
        self.bcl_path = os.path.join(self.data_path, "Data/Intensities/BaseCalls/")
        self.seq_mode = self.option("sample_info").prop["seq_mode"]
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/bin')

    def get_sample_sheet_csv(self):
        """
        生成csv表，供bcl2fastq使用
        """
        self.sample_sheet_csv = os.path.join(self.work_dir, "sample_sheet.csv")
        seq_type = self.option("sample_info").prop["seq_type"]
        self.option("sample_info").create_sample_sheet(seq_type, self.sample_sheet_csv)

    def bcl2fastq(self):
        """
        运行bcl2fastq
        """
        cmd = self.bcl2fastq_path + ' -i ' + self.bcl_path + ' -o ' \
              + self.output_dir + ' --sample-sheet ' + self.sample_sheet_csv + ' --use-bases-mask ' + self.seq_mode\
              + ' --ignore-missing-bcl ' + '-R ' + self.data_path + ' -r 4 -w 4 -d 2 -p 10 --barcode-mismatches 0'
        self.logger.info(cmd)
        command = self.add_command('bcl2fastq', cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行bcl2fastq完成")
        else:
            self.set_error("运行bcl2fastq出错")

    def run(self):
        super(Bcl2fastqV2Tool, self).run()
        self.get_sample_sheet_csv()
        self.bcl2fastq()
        self.end()

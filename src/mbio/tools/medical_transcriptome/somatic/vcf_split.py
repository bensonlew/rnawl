# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# modified 2020.10.26

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import os
import subprocess
from biocluster.config import Config
import datetime

class VcfSplitAgent(Agent):
    """
    基因比对接口
    """

    def __init__(self, parent):
        super(VcfSplitAgent, self).__init__(parent)
        options = [
            {"name": "combine_vcf", "type": "infile", "format": "ref_rna_v2.common"},  # 输入的bam
            {"name": "line_num", "type": "int", "default": 10000},  # 输入格式  bam/cram 20191231
        ]
        self.add_option(options)
        self._memory_increase_step = 10

    def check_options(self):
        # if not self.option("bam_file"):
        #     raise OptionError("请设置bam路径")
        # if not self.option("fa_file"):
        #     raise OptionError("请设置ref.fa路径")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(VcfSplitAgent, self).end()


class VcfSplitTool(Tool):
    def __init__(self, config):
        super(VcfSplitTool, self).__init__(config)

    def split_vcf(self):
        self.logger.info('开始分割文件')
        self.logger.info('每个子文件{}行'.format(self.option("line_num")))
        flag = 0
        name = 1
        # 存放数据
        dataList = []
        header = []
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        with open(self.option("combine_vcf").prop["path"], 'r') as f_source:
            for line in f_source:
                if line.startswith("#"):
                    header.append(line)
                else:
                    flag += 1
                    dataList.append(line)
                    if flag == 10000:
                        with open(os.path.join(self.output_dir, "split_{}".format(str(name)) + ".vcf"), 'w+') as f_target:
                            for head in header:
                                f_target.write(head)
                            for data in dataList:
                                f_target.write(data)
                        name += 1
                        flag = 0
                        dataList = []

        # 处理最后一批行数少于10000行的
        with open(os.path.join(self.output_dir, "split_{}".format(str(name)) + ".vcf"), 'w+') as f_target:
            for head in header:
                f_target.write(head)
            for data in dataList:
                f_target.write(data)

        print("完成。。。。。")
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    def run(self):
        super(VcfSplitTool, self).run()
        self.split_vcf()
        self.end()



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "split_vcf" + str(random.randint(1, 10000))+"yyyy",
            "type": "tool",
            "name": "medical_transcriptome.somatic.vcf_split",
            "instant": False,
            "options": dict(
                combine_vcf="/mnt/ilustre/users/sanger-dev/workspace/20201105/MedicalTranscriptome_v6dpivfr84k4ooq2lmvh47dpfs/CallSnpIndelSentieon/VcfFilterGatk/output/final.vcf",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
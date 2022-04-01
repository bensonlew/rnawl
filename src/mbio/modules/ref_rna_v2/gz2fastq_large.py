# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
import shutil
import unittest

import pandas as pd
from Bio import SeqIO

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile
from multiprocessing import Pool
import  glob
from collections import OrderedDict
import re


FAKE_FASTQ = """@ST-E00575:252:HJ7HCCCXY:7:1101:8674:1713 1:N:0:CGTACG
GNCCCAACTTGCCATCAAGGATATCTATCTCGGCAACCGCTTCGTTAAATGTCTCTTCGTGGTCAGCTTCAATAGCCAATTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTGAAA
+
A#AF-7FFJFJJAJJJJJJJJJJJJJJ-AJFJJJJJJJJJJJJJ<JJJJAJJ7FJJFFFJFFFAFJJJJJAJJJJFFJJJJJJFJFFJJFJJJJFJJJ<FJJJJFJFJJJJJFJJJJJJJJJJJJJJ7FJJJJJJJJ7AJJJJJJJFJJJJ
@ST-E00575:252:HJ7HCCCXY:7:1101:9810:1713 1:N:0:CGTACG
CNGTAATCGTTTGTGGCGTTAGAAATAAAGCCTCAGCCGCCCCGACGACAGAGCCTTCCTTACAAACTTGCCAAAAATAATAAAGATGATTGAAATTGATGTGCGACATTCGCATGTTGTTATCCCCAGATCGGAAGAGCCACACGTCTGA
+
A#AFFJJJJJJJJJJJJJJJJJFJJJJ<JJJJJJJJJJJJAJJJJJJJJJJJJJJJ<-JJJJJJJJJJJJAAJJJJJJJJJJJJJJJJFFJJJJJJJJFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJA
@ST-E00575:252:HJ7HCCCXY:7:1101:12408:1713 1:N:0:CGTACG
CNCACCAATCATCCTGGACTGGCTCTCAATCTCCATCCTGGAGGTGTCCTTTGTTTCTTCCTGAAACATCCCTTCACTCATCCTAAGCAGTCCCTGAGTCCTTCATCCTGAAGTGGCACCATCCTGATACCGTCCTTTAGATCGGAAGAGC
+
A#AFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJFFJFJJJJFJJJJJFFJJJJJ<AJJJJJJJJJJJJJJJJJJJFJFJJJ<AFFJFJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJF
"""


class Gz2fastqLargeModule(Module):
    """
    针对1000个样本以上的大样本项目进行解压缩操作,通过并行解决时间过长的问题
    """

    def __init__(self, work_id):
        super(Gz2fastqLargeModule, self).__init__(work_id)
        options = [
            {'name': 'fastq_path', 'type': 'string', 'default': ''}
        ]
        self.add_option(options)
        self.samples = {}
        self.tools = []
        self.list_info_new = ""

    def check_options(self):
        """
        检查参数
        """
        # if self.option('fq_type') not in ['PE', 'SE']:
        #     raise OptionError("请说明序列类型，PE or SE?", code = "33705201")
        return True

    def gz2fastq_run(self):
        file_list = self.get_list_file()
        for n,file in enumerate(file_list):
            fastq_stat = self.add_tool("ref_rna_v3.large.pigzfastq2fastq")
            options = {
                "fastq_path" : self.option("fastq_path"),
                'fastq_info': file
            }
            fastq_stat.set_options(options)
            self.tools.append(fastq_stat)
        self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def get_list_file(self):
        sample_file =OrderedDict()
        fastq_file_list = []
        fq_dir = self.option("fastq_path")
        list_info = os.path.join(fq_dir, "list.txt")
        with open(list_info, 'r') as list_r:
            list_info = list_r.read()
            list_info_new = ''
            for line in list_info.split('\n'):
                line1 = line.strip().split('\t')
                if line1[0].split(".")[-1] in ["gz"]:
                    gzed = ".".join(line1[0].split(".")[:-2]) + ".fastq"
                    list_info_new += gzed + '\t' + '\t'.join(line1[1:]) + '\n'
                    fastq_file_list.append((line))
                else:
                    if line:
                        list_info_new += '\t'.join(line) + '\n'
        self.list_info_new = list_info_new
        self.logger.info(fastq_file_list)

        def yield_fastq_list(fastq_file_list,chunksize =20):
            file_info = ""
            n = 0
            for sample in fastq_file_list:
                if  n == chunksize:
                    yield file_info
                    file_info = ""
                    n=0
                file_info += "{}\n".format(sample)
                n+=1
            if file_info:
                yield file_info
        if os.path.exists(os.path.join(self.work_dir,"fq_list")):
            shutil.rmtree(os.path.join(self.work_dir,"fq_list"))
        os.makedirs(os.path.join(self.work_dir,"fq_list"))
        for n,file_info in enumerate(yield_fastq_list(fastq_file_list)):
            with open(os.path.join(self.work_dir,"fq_list","fq_list"+"_{}".format(str(n))), "wb") as w:
                    w.write(file_info)
        file_list = glob.glob(self.work_dir +"/fq_list/fq_list_*")
        return sorted(file_list)


    def run(self):
        super(Gz2fastqLargeModule, self).run()
        self.gz2fastq_run()


    def set_output(self):
        self.logger.info("set output")
        os.system('rm -rf ' + self.output_dir)
        os.system('mkdir ' + self.output_dir)
        os.remove(self.option("fastq_path") + '/list.txt')
        with open(self.option("fastq_path") + '/list.txt', 'w') as list_w:
            list_w.write(self.list_info_new.strip('\n') + '\n')

        self.logger.info("done")
        self.end()

    def end(self):
        super(Gz2fastqLargeModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        data = {
            "id": "large_fastq_stat" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "module",
            "name": "ref_rna_v2.gz2fastq_large",
            "instant": False,
            "options": dict(
                # fastq_path='/mnt/ilustre/users/isanger/workspace/20210315/Refrna_n34u_49qekvdhgok1ed33c1m5hf/remote_input/fastq_dir/rawdata',
                fastq_path = "/mnt/ilustre/users/isanger/workspace/20210301/Denovorna_majorbio_324951/remote_input/fastq_dir/Rawdata_FX2021022200180"
                # fq_type = "PE"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

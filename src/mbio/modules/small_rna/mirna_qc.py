# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20171213

import os
import re
import time
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import unittest
import gevent

class MirnaQcModule(Module):
    """
    miRNA质控的module,用于miRNA质控，多个fastq文件进行质控
    是否去除前三个碱基，去接头，去低值，去未知碱基序列、去过长过短序列
    """
    def __init__(self, work_id):
        super(MirnaQcModule, self).__init__(work_id)
        options = [
            {"name": "list_file", "type": "infile", "format": "datasplit.list_file"},  # fastq文件及对应的样本
            {"name": "length", "type": "int", "default": 18},  # 最小序列长度，丢弃比此值更短的序列
            {"name": "skip_qc", "type": "string", "default": "no"}, # other 跳过质控

            {"name": "config_file", "type": "outfile", "format": "small_rna.ini"},  # 质控序列配置文件
            {"name": "adapter", "type": "string", "default": 'AGATCGGAAGAGCACACGTC'},  # 接头序列，如果微量建库，则改为AAAAAA
            {"name": "phred_score", "type": "string", "default": "20"},  # Phred得分（在0到40之间）
            {"name": "minlen", "type": "string", "default": "18"},  # 最短序列
            {"name": "fastq_format", "type": "string", "default": "sanger"}, # 质量体系"sanger", "solexa", "illumina", "Q33", "Q64"
            {"name": "max_length", "type": "string", "default": "32"},  # 最长序列
            {"name": "extract_length", "type": "int", "default": 75},  # 保留原始序列的长度
            {"name": "cut_left", "type": "int", "default": 0}, #miRNA是否要切除前3bp, 或4bp, 保留前51bp
            {"name": "cut_tail", "type": "int", "default": 0}, #miRNA是否要切除去接头尾部的4bp
            {"name": "rawdata", "type": "outfile", "format": "sequence.fastq_dir"}, # 输出原始数据文件路径
            {"name": "cleandata", "type": "outfile", "format": "sequence.fastq_dir"} # 输出质控数据文件路径
        ]
        self.add_option(options)
        self.step.add_steps("single_mirna_qc")
        self.fastqs = {}

    def set_step(self, event):
        if "start" in event["data"].keys():
            event["data"]["start"].start()
        if "end" in event["data"].keys():
            event["data"]["end"].finish()
        self.step.update()

    def check_options(self):
        if not self.option("list_file").is_set:
            raise OptionError("必须设置list文件")
        '''
        with open(self.option("list_file").prop['path'], 'r') as f, open(self.option("list_file").prop['path'] + '.clean', 'w') as list_clean:
            for line in f.readlines():
                line = line.strip()
                if len(line.split("\t")) >= 3:
                    if line.split("\t")[2] == 'r':
                        pass
                    else:
                        list_clean.write("{}\t{}\n".format(line.split("\t")[0], line.split("\t")[1]))
        self.option("list_file", self.option("list_file").prop['path'] + '.clean')
        '''

    def get_sample_list(self):
        '''
        获取样本列表
        '''
        sample_list = []
        with open(self.option("list_file").prop['path'], 'r') as list_f:
            for line in list_f.readlines():
                item = line.strip().split()
                if item[1] in sample_list:
                    pass
                else:
                    sample_list.append(item[1])
        return sample_list


    def run_mirna_qc(self):
        options = {
            "length": self.option("length"),
            "skip_qc": self.option("skip_qc"),
            "adapter": self.option("adapter"),
            "phred_score": self.option("phred_score"),
            "min_length": self.option("minlen"),
            "max_length": self.option("max_length"),
            "fastq_format": self.option("fastq_format"),
            "length_contain": self.option("extract_length"),
            "cut_left": self.option("cut_left"),
            "cut_tail": self.option("cut_tail")
        }

        modules = []
        for s in self.fastqs.keys():
            options["fastq"] = self.fastqs[s][0]
            options["sample_name"] = s
            self.single_mirna_qc  = self.add_module("small_rna.single_mirna_qc")
            self.single_mirna_qc.set_options(options)
            self.single_mirna_qc.on('start', self.set_step, {'start': self.step.single_mirna_qc})
            self.single_mirna_qc.on('end', self.set_step, {'end': self.step.single_mirna_qc})
            self.single_mirna_qc.on('end', self.set_output, 'mirna_qc_{}'.format(s))
            modules.append(self.single_mirna_qc)
        if len(modules) == 1:
            modules[0].on("end", self.end)
        else:
            self.on_rely(modules, self.end)
        for m in modules:
            m.run()

    def set_output(self, event):
        obj = event["bind_object"]
        sample = event["data"].split("mirna_qc_")[1]
        olddir = obj.output_dir
        if not os.path.isdir(olddir):
            raise Exception("需要移动到output目录的文件夹不存在")
        if os.path.exists(self.output_dir + "/raw_data"):
            pass
        else:
            os.mkdir(self.output_dir + "/raw_data")

        if os.path.exists(self.output_dir + "/clean_data"):
            pass
        else:
            os.mkdir(self.output_dir + "/clean_data")

        if os.path.exists(self.output_dir + "/raw_data"):
            pass
        else:
            os.mkdir(self.output_dir + "/raw_data")

        if os.path.exists(os.path.join(self.output_dir, 'raw_data', sample + ".fq")):
            os.remove(os.path.join(self.output_dir, 'raw_data', sample + ".fq"))
        os.link(obj.option("raw_fastq").prop['path'] , os.path.join(self.output_dir, 'raw_data', sample + ".fq"))
        for f in os.listdir(olddir):
            if re.search(r".+.fasta", f):
                f1 = sample + "_clean.fasta"
            elif re.search(r".+.xls", f):
                f1 = sample + "_qc_stat.xls"
            elif re.search(r".+.trimmed", f):
                f1 = sample + "_clip_s.fastq.trimmed"
            elif re.search(r".+.length.txt", f):
                f1 = sample + "_clean.length.txt"
            elif re.search(r".+.fq", f):
                continue
            else:
                continue
            self.logger.info("link {} to {}".format(f, f1))
            if os.path.exists(os.path.join(self.output_dir, 'clean_data', f1)):
                os.remove(os.path.join(self.output_dir, 'clean_data', f1))
            os.link(os.path.join(olddir, f), os.path.join(self.output_dir, 'clean_data', f1))


    def set_qc_file(self):
        '''
        生成config文件
        '''
        sample_list = self.get_sample_list()
        sample_num = {k:v for v,k in enumerate(sample_list)}
        with open(os.path.join(self.work_dir, 'qc_file.config'), 'w') as f:
            f.write("[FASTA]\n")
            i = 0
            id2name = dict()
            ids = list()
            for s in sorted(self.fastqs.keys(), key=lambda x:sample_num[x]):
                i = i + 1
                sample_id = 'S' + '{:x}'.format(i).zfill(2)
                id2name[sample_id] = s
                ids.append(sample_id)
                f.write("{}={}\n".format(sample_id, os.path.join(self.output_dir, 'clean_data', s + "_clean.fasta")))
            f.write("\n[NAME]\n")
            for i in ids:
                f.write("{}={}\n".format(i, id2name[i]))
        if os.path.exists(os.path.join(self.output_dir, 'clean_data', 'qc_file.config')):
            os.remove(os.path.join(self.output_dir, 'clean_data', 'qc_file.config'))
        os.link(os.path.join(self.work_dir, 'qc_file.config'), os.path.join(self.output_dir, 'clean_data', 'qc_file.config'))
        self.option('config_file', os.path.join(self.output_dir,'clean_data', 'qc_file.config'))
        # 生成trim后的fq列表
        with open(os.path.join(self.work_dir, 'list.txt'), 'w') as f:
            for s in sorted(self.fastqs.keys(), key=lambda x:sample_num[x]):
                f.write("{}\t{}\n".format(os.path.join(self.output_dir, 'clean_data', s + "_clip_s.fastq.trimmed"), s))
        if os.path.exists(os.path.join(self.output_dir, 'clean_data', 'list.txt')):
            os.remove(os.path.join(self.output_dir, 'clean_data', 'list.txt'))
        os.link(os.path.join(self.work_dir, 'list.txt'), os.path.join(self.output_dir, 'clean_data', 'list.txt'))
        with open(os.path.join(self.work_dir, 'list.raw.txt'), 'w') as f:
            for s in sorted(self.fastqs.keys(), key=lambda x:sample_num[x]):
                f.write("{}\t{}\n".format(os.path.join(self.output_dir, 'raw_data', s + ".fq"), s))
        if os.path.exists(os.path.join(self.output_dir, 'raw_data', 'list.txt')):
            os.remove(os.path.join(self.output_dir, 'raw_data', 'list.txt'))
        os.link(os.path.join(self.work_dir, 'list.raw.txt'), os.path.join(self.output_dir, 'raw_data', 'list.txt'))

        self.option("rawdata", os.path.join(self.output_dir, 'raw_data'))
        self.option("cleandata", os.path.join(self.output_dir, 'clean_data'))

    def run(self):
        super(MirnaQcModule, self).run()
        self.fastqs = self.option("list_file").prop["samples"]
        self.run_mirna_qc()

    def end(self):
        self.set_qc_file()
        super(MirnaQcModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime

        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA'
        data = {
            "id": "mirna_qc" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "module",
            "name": "small_rna.mirna_qc",
            "instant": False,
            "options": dict(
                list_file = test_dir + "/" + "data0/list.txt",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()

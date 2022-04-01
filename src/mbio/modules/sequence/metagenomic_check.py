#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == shaohua.yuan

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import unittest
# from biocluster.wpm.client import *
import pandas as pd


class MetagenomicCheckModule(Module):
    """
    宏基因组工作流2非冗余基因集输入时
    若文件为压缩文件，解压缩，若为非压缩文件，直接链接
    将上传list文件改名为list.txt并放入解压后的fq_dir文件夹
    严格检查list.txt文件和上传的的非冗余基因集文件
    """
    def __init__(self, work_id):
        super(MetagenomicCheckModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},
            {'name': 'clean_list', 'type': 'infile', 'format': 'sequence.profile_table'}, # 工作流2使用的clean_data的list.txt
            {'name': 'insertsize', 'type': 'infile', 'format': 'sample.insertsize_table'},  # 插入片段长度表
            {"name": "out_fastq_dir", "type": "outfile", "format": "sequence.fastq_dir"},
        ]
        self.add_option(options)
        self.base_info = self.add_tool("meta.qc.base_info")
        self.qc_stat = self.add_tool("sequence.raw_qc_stat")
        self.tools = []
        self.link_fq = {}

    def check_options(self):
        """
        检查参数
        """
        if not self.option("fastq_dir").is_set:
            raise OptionError("must input fastq_dir!", code="24000701")
        if not self.option("clean_list").is_set:
            raise OptionError("must input clean_list!", code="24000702")
        return True

    def run(self):
        super(MetagenomicCheckModule, self).run()
        self.run_ungiz()

    def check_list(self):
        list_path = self.option("clean_list").prop["path"]
        sample = {}
        with open(list_path, "rb") as l:
            for line in l:
                line = line.strip().split()
                if len(line) == 3:
                    sample_path = os.path.join(self.option("fastq_dir").prop['path'], line[0])
                    if line[1] not in sample:
                        sample[line[1]] = {line[2]: sample_path}
                    elif line[2] in sample[line[1]].keys():
                        sample[line[1]][line[2]] = sample[line[1]][line[2]] + " " + sample_path
                    else:
                        sample[line[1]][line[2]] = sample_path
                else:
                     raise OptionError('list.txt file format is error', code="24000703")
        return sample

    def run_ungiz(self):
        self.samples = self.check_list()
        samples = self.samples
        if not os.path.exists(self.output_dir + "/clean_dir"):
            os.mkdir(self.output_dir + "/clean_dir")
        new_list = os.path.join(self.output_dir + "/clean_dir","list.txt")
        outf = open(new_list,"w")
        self.result_path = os.path.join(self.work_dir, "ungiz_dir")
        if not os.path.exists(self.result_path):
            os.mkdir(self.result_path)
        for sample in samples:
            for d in samples[sample]:
                direct = ""
                if d == "r":
                    direct = "2"
                elif d == "l":
                    direct = "1"
                elif d == "s":
                    direct = "s"
                else:
                    raise OptionError("序列的方向不对，必须为：l/r/s", code="24000704")
                fq_final_path_name = sample + "." + direct + ".fq"
                #fq_final_path = os.path.join(self.output_dir, fq_final_path_name)
                outf.write(fq_final_path_name + "\t" + sample + "\t" + d + "\n")
                last = samples[sample][d].split(".")
                if last[-1] == "gz":
                    gunzip_fastq = self.add_tool('sequence.fastq_ungz')
                    gunzip_fastq.set_options({
                        "fastq": samples[sample][d],
                        "sample_name": sample,
                        "direction": direct,
                        "result_path": self.result_path
                    })
                    self.tools.append(gunzip_fastq)  # 修复缩进 guhaidong 20171204
                elif last[-1] == "fq":
                    self.link_fq[sample][direct] = samples[sample][d]
        if len(self.tools) == 0:
            self.set_output()
        elif len(self.tools) > 1:
            self.on_rely(self.tools, self.run_base_info)
        else:
            self.tools[0].on('end', self.run_base_info)
        for tool in self.tools:
            tool.run()

    def run_base_info(self):
        self.set_clean()
        self.base_info.set_options({"fastq_path":  os.path.join(self.work_dir, "ungiz_dir")})
        self.base_info.on('end', self.run_qc_stat)
        self.base_info.run()

    def run_qc_stat(self):
        options = {
            "base_info_dir": os.path.join(self.base_info.output_dir, "base_info"),
            "sickle_dir": os.path.join(self.output_dir, "clean_dir")
        }
        self.qc_stat.set_options(options)
        self.qc_stat.on('end', self.set_output)
        self.qc_stat.run()

    def set_output(self):
        base_info_dir = os.path.join(self.base_info.output_dir,"base_info")
        self.linkdir(base_info_dir,"base_info")
        qc_stat = os.path.join(self.qc_stat.output_dir,"reads.cleanData.stat")
        new_qc = os.path.join(self.output_dir, "optimize_reads.stat.xls")
        if os.path.exists(new_qc):
            os.remove(new_qc)
        if self.option("insertsize").is_set:
            self.inset = self.option("insertsize").prop["path"]
            inset_table = pd.read_table(self.inset, sep="\t",header=None)
            inset_table.columns = ["Sample","Insertsize"]
            qc_stat_table = pd.read_table(qc_stat, sep="\t",header=0)
            new_qc_table = pd.merge(qc_stat_table,inset_table,left_on = "#Sample",right_on="Sample",how="left")
            new_qc_table = new_qc_table.drop(columns=["Sample"])
            new_qc_table.to_csv(new_qc,sep="\t",index=False)
        else:
            os.link(qc_stat, new_qc)
        self.end()

    def set_clean(self):
        self.logger.info("start set clean dir")
        if len(self.link_fq) >0:
            clean_dir = os.path.join(self.output_dir, "clean_dir")
            if not os.path.exists(clean_dir):
                os.mkdir(clean_dir)
            for eachsample in self.link_fq:
                for direct in self.link_fq[eachsample]:
                    old_file = self.link_fq[eachsample][direct]
                    new_name = eachsample + "." + direct + ".fq"
                    new_path = os.path.join(clean_dir, new_name)
                if os.path.exists(new_path):
                    os.remove(new_path)
                os.link(name, new_path)
        old_dir = self.result_path
        self.linkdir(old_dir, "clean_dir")

    def end(self):
        super(MetagenomicCheckModule, self).end()

    def linkdir(self, dirpath, new_dir=None):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param new_dir: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        if not new_dir:
            newdir = self.output_dir
        else:
            newdir = os.path.join(self.output_dir, new_dir)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                file_name = os.listdir(oldfiles[i])
                os.mkdir(newfiles[i])
                for file_name_ in file_name:
                    os.link(os.path.join(oldfiles[i], file_name_), os.path.join(newfiles[i], file_name_))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            #"id": "metag_checkfile" + str(random.randint(1, 10000)),
            "id": "metag_checkfile",
            "type": "module",
            "name": "sequence.metagenomic_check",
            "instant": True,
            "options": dict(
                fastq_dir="/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/project/metagenomic_v2/workflow2_data/onesample/Clean_gzip",
                clean_list="/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/project/metagenomic_v2/workflow2_data/onesample/Clean_gzip/list.txt",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()
        '''
        worker = worker_client()
        result = worker.add_task(data)
        print result
        '''


if __name__ == '__main__':
    unittest.main()
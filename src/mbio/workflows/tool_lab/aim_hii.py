# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.tool_lab.common_function import fasta_dir_sort
import os
import re, shutil
import subprocess
import gevent

class AimHiiWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AimHiiWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"}, # 参考基因组
            {"name": "insert_fa", "type": "infile", "format": "sequence.fasta"}, # 插入序列
            {'name': 'in_fastq', 'type': 'infile', 'format': 'tool_lab.raw_dir'},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.tools = []

    def check_options(self):
        if not self.option("ref_fasta").is_set:
            raise OptionError("请传入参考基因组！")
        if not self.option("insert_fa").is_set:
            raise OptionError("请传入插入序列！")
        if not self.option("in_fastq").is_set:
            raise OptionError("请传入原始测序序列文件夹！")

    def run(self):

        self.run_ungiz()
        self.run_aimhii()
        super(AimHiiWorkflow, self).run()

    def run_aimhii(self):
        """
        目的片段查找
        :return:
        """
        reslut_path = os.path.join(self.work_dir, "ungiz_dir")
        for sample in self.sample_dict:
            self.aimhii = self.add_tool("tool_lab.aimhii")
            if len(self.sample_dict[sample]) == 2:
                opts = {
                    "ref_fa": self.option("ref_fasta"),
                    "insert_fa": self.option("insert_fa"),
                    "clean_fq_1": reslut_path + "/" + self.sample_dict[sample][0].strip(".gz"),
                    "clean_fq_2": reslut_path + "/" + self.sample_dict[sample][1].strip(".gz"),
                    "sample": sample,
                }
            else:
                opts = {
                    "ref_fa": self.option("ref_fasta"),
                    "insert_fa": self.option("insert_fa"),
                    "clean_fq_1": reslut_path + "/" + self.sample_dict[sample][0].strip(".gz"),
                    "sample": sample,
                }
            self.aimhii.set_options(opts)
            self.tools.append(self.aimhii)
            self.logger.info(self.tools)
            if self.tools:
                if len(self.tools) > 1:
                    self.on_rely(self.tools, self.set_output)
                elif len(self.tools) == 1:
                    self.tools[0].on('end', self.set_db)
                for tool in self.tools:
                    self.logger.info(tool)
                    gevent.sleep(1)
                    tool.run()

    def run_ungiz(self):
        self.sample_dict = {}
        with open(self.option("in_fastq").prop["path"] + "/list.txt", "r") as f:
            data = f.readlines()
            for i in data[1:]:
                if i:
                    if len(i.split("\t")) < 2:
                        self.set_error("list文件格式有误！")
                    if len(i.strip().split("\t")[1].strip().split(";")) >= 2:
                        self.sample_dict[i.split("\t")[0]] = [i.strip().split("\t")[1].split(";")[0],
                                                              i.strip().split("\t")[1].split(";")[1]]
                    else:
                        self.sample_dict[i.split("\t")[0]] = [i.strip().split("\t")[1].split(";")[0]]
        reslut_path = os.path.join(self.work_dir, "ungiz_dir")
        if not os.path.exists(reslut_path):
            os.mkdir(reslut_path)
        for sample in self.sample_dict:
            shutil.copy(self.option("in_fastq").prop["path"] + "/" + self.sample_dict[sample][0],reslut_path+"/"+self.sample_dict[sample][0])
            os.system("gzip -d {}".format(reslut_path+"/"+self.sample_dict[sample][0]))
            if len(self.sample_dict[sample]) >= 2:
                shutil.copy(self.option("in_fastq").prop["path"] + "/" + self.sample_dict[sample][1],
                        reslut_path + "/" + self.sample_dict[sample][1])
                os.system("gzip -d {}".format(reslut_path + "/" + self.sample_dict[sample][1]))

    def set_db(self):
        with open(self.output_dir + "/InsertionResult.csv","w") as t:
            t.write("sample,AIM Num,type,ref_chrom,left_start_b1,left_junction_b1,left_numreads"
                    ",gap_length,insert_length,insert_chrom,insert_start,insert_end,insert_strand"
                    ",right_junction_b1,right_end_b1,right_numreads\n")
            for tool in self.tools:
                num = 0
                for file in os.listdir(tool.output_dir):
                    if file.endswith("result.csv"):
                        a = open(tool.output_dir + "/" + file)
                        data = a.readlines()
                        for i in data[1:]:
                            num += 1
                            tmp = i.rstrip().split(",")
                            t.write(file.rstrip(".result.csv") + "," + "AIM" + str(num).zfill(3) + "," + ",".join(tmp) + "\n")
                    if file.endswith(".pdf"):
                        new_name = file.split("readplot_")[0] + "AIM" + str(int(file.split("readplot_")[1].split(".")[0]) + 1).zfill(3) + ".pdf"
                        if os.path.exists(self.output_dir + "/" + new_name):
                            os.remove(self.output_dir + "/" + new_name)
                        os.link(tool.output_dir + "/" + file,self.output_dir + "/" + new_name)
        remote_dir = self._sheet.output
        #remote_dir = self.output_dir
        api_main = self.api.api("tool_lab.aimhii")
        if self.option('main_id'):
            main_id =  self.option('main_id')
        else:
            main_id = api_main.add_aimhii_main()
        api_main.add_aimhii_detail(main_id=main_id, file_path=self.output_dir + "/InsertionResult.csv",picture_path=remote_dir)
        api_main.add_aimhii_main(main_id)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "插入位点结果目录", 0],
            ["./*.pdf", "pdf", "插入位点示意图", 0],
            ["./*.csv", "csv", "插入位点鉴定结果表", 0]
        ])
        super(AimHiiWorkflow, self).end()
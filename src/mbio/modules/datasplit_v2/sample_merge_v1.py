# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190125

import os
import re
import json
import subprocess
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class SampleMergeModule(Module):
    """
    合并fastq文件
    """
    def __init__(self, work_id):
        super(SampleMergeModule, self).__init__(work_id)
        options = [
            {"name": "sample_list", "type": "infile", "format": "datasplit.sample_merge_list"},
            {"name": "merge_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.raw_r1, self.raw_r2 = [], []
        self.clean_r1, self.clean_r2 = [], []
        self.new_paths, self.download_files = {}, {}
        self.sample_name = ""

    def check_options(self):
        if not self.option("sample_list").is_set:
            raise OptionError("缺少sample_list,请设置")

    def get_sample_path(self):
        """
        获取raw和clean对应的r1、r2,若是s3开头换成要下载的对应文件
        """
        with open(self.option("sample_list").prop["path"], "rb") as f:
            lines = f.readlines()
            for i in range(1, len(lines)):
                item = lines[i].strip().split("\t")
                self.sample_name = item[2]
                raw_path = self.get_new_files(item[3], i, "raw")
                # clean_path = self.get_new_files(item[4], i, "clean")  # 20201221 clean数据暂时不处理
                self.raw_r1.append(raw_path[0])
                if len(raw_path) == 2:
                    self.raw_r2.append(raw_path[1])
                # self.clean_r1.append(clean_path[0])
                # if len(clean_path) == 2:
                #     self.clean_r2.append(clean_path[1])

    def get_new_files(self, paths, i, type):
        new_paths_key = []
        old_paths = paths.split(";")
        for j in range(len(old_paths)):
            new_paths_key.append(type+"_"+str(i)+"_"+str(j))
            self.new_paths[type+"_"+str(i)+"_"+str(j)] = old_paths[j]
            if old_paths[j].startswith("s3"):
                self.download_files[type+"_"+str(i)+"_"+str(j)] = old_paths[j]
        return new_paths_key

    def run_fastq_download(self):
        self.tools = []
        for key in self.download_files.keys():
            options = ({
                "fastq": self.download_files[key]
            })
            self.fastq_download = self.add_tool("datasplit_v2.fastq_download")
            self.fastq_download.set_options(options)
            self.fastq_download.on("end", self.get_fastq_path, key)
            self.tools.append(self.fastq_download)
        if len(self.tools) == 0:
            self.run_merge_fastq()
        elif len(self.tools) == 1:
            self.tools[0].on("end", self.run_merge_fastq)
            self.tools[0].run()
        else:
            self.on_rely(self.tools, self.run_merge_fastq)
            for tool in self.tools:
                tool.run()

    def get_fastq_path(self, event):
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            self.new_paths[event["data"]] = os.path.join(obj.output_dir, f)

    def run_merge_fastq(self):
        """
        分别合并raw和clean对应的r1、r2
        """
        self.merge_fastq(self.raw_r1, self.sample_name + ".R1.raw.fastq.gz")
        if self.raw_r2:
            self.merge_fastq(self.raw_r2, self.sample_name + ".R2.raw.fastq.gz")
        # self.merge_fastq(self.clean_r1, self.sample_name + ".clean.1.fastq.gz")
        # if self.clean_r2:
        #     self.merge_fastq(self.clean_r2, self.sample_name + ".clean.2.fastq.gz")
        self.set_db()

    def merge_fastq(self, key_list, file_name):
        """
        合并fastq文件
        """
        path_list = []
        for key in key_list:
            path_list.append(self.new_paths[key])
        out_file = os.path.join(self.output_dir, file_name)
        cmd = "cat " + " ".join(path_list) + " > " + out_file
        self.logger.info("cmd:%s", cmd)
        command = subprocess.Popen(cmd, shell=True)
        command.communicate()
        if command.returncode == 0:
            self.logger.info("%s cat完成" % " ".join(path_list))
        else:
            self.set_error("%s cat失败" % " ".join(path_list))

    def set_db(self):
        """
        结果导表
        """
        if self.option("merge_id"):
            self.logger.info("开始进行结果导表")
            datasplit_api = self.api.api("datasplit.datasplit_new")
            data_json = os.path.dirname(self.work_dir) + "/data.json"
            s3_output_dir = json.loads(open(data_json).read())["output"]
            datasplit_api.sg_split_specimen_merge(self.option("merge_id"), self.output_dir, s3_output_dir)
        self.end()

    def run(self):
        super(SampleMergeModule, self).run()
        self.get_sample_path()
        self.run_fastq_download()

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SampleMergeModule, self).end()

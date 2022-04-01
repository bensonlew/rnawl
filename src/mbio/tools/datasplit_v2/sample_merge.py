# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modified: 20201105

import os
import json
import subprocess
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from s3transfer.s3 import S3TransferManager


class SampleMergeAgent(Agent):
    """
    合并fastq文件
    """
    def __init__(self, parent=None):
        super(SampleMergeAgent, self).__init__(parent)
        options = [
            {"name": "sample_list", "type": "infile", "format": "datasplit.sample_merge_list"},
            {"name": "merge_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sample_list").is_set:
            raise OptionError("缺少sample_list,请设置")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SampleMergeAgent, self).end()


class SampleMergeTool(Tool):
    def __init__(self, config):
        super(SampleMergeTool, self).__init__(config)
        self.raw_r1, self.raw_r2 = [], []
        self.clean_r1, self.clean_r2 = [], []
        self.download_files = {}
        self.sample_name = ""

    def get_new_files(self, paths, id):
        new_paths = []
        old_paths = paths.split(";")
        for i in range(len(old_paths)):
            f = old_paths[i]
            if f.startswith("s3:"):
                path = os.path.join(self.work_dir, os.path.basename(f))
                if path in new_paths:
                    path = os.path.join(self.work_dir, id[i]+"--"+os.path.basename(f))
                self.download_files[f] = path
                new_paths.append(path)
            else:
                new_paths.append(f)
        return new_paths

    def get_sample_path(self):
        """
        获取raw和clean对应的r1、r2,若是s3开头换成要下载的对应文件
        """
        with open(self.option("sample_list").prop["path"], "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                self.sample_name = item[2]
                raw_path = self.get_new_files(item[3], item[1])
                clean_path = self.get_new_files(item[4], item[1])
                self.raw_r1.append(raw_path[0])
                if len(raw_path) == 2:
                    self.raw_r2.append(raw_path[1])
                self.clean_r1.append(clean_path[0])
                if len(clean_path) == 2:
                    self.clean_r2.append(clean_path[1])

    def run_download_files(self):
        """
        下载对象存储文件
        """
        if self.download_files:
            bcl_type = "tsg"
            if self.download_files.keys()[0].startswith("s3://datasplit"):
                bcl_type = "sanger"
            transfer = S3TransferManager(base_path=None, use_db=False, overwrite=True, log=True, etag=0, bcl_type=bcl_type)
            for s3_path in self.download_files:
                transfer.add(s3_path, self.download_files[s3_path])
            transfer.wait()

    def merge_fastq(self, path_list, file_name):
        """
        合并fastq文件
        """
        out_file = os.path.join(self.output_dir, file_name)
        cmd = "cat " + " ".join(path_list) + " > " + out_file
        self.logger.info("cmd:%s", cmd)
        command = subprocess.Popen(cmd, shell=True)
        command.communicate()
        if command.returncode == 0:
            self.logger.info("%s cat完成" % " ".join(path_list))
        else:
            self.set_error("%s cat失败" % " ".join(path_list))

    def run_merge_fastq(self):
        """
        分别合并raw和clean对应的r1、r2
        """
        self.merge_fastq(self.raw_r1, self.sample_name + ".R1.raw.fastq.gz")
        if self.raw_r2:
            self.merge_fastq(self.raw_r2, self.sample_name + ".R2.raw.fastq.gz")
        self.merge_fastq(self.clean_r1, self.sample_name + ".clean.1.fastq.gz")
        if self.clean_r2:
            self.merge_fastq(self.clean_r2, self.sample_name + ".clean.2.fastq.gz")

    def set_db(self):
        """
        结果导表
        """
        self.logger.info("开始进行结果导表")
        datasplit_api = self.api.api("datasplit.datasplit_new")
        data_json = os.path.dirname(self.work_dir) + "/data.json"
        s3_output_dir = json.loads(open(data_json).read())["output"]
        datasplit_api.sg_split_specimen_merge(self.option("merge_id"), self.output_dir, s3_output_dir)

    def run(self):
        super(SampleMergeTool, self).run()
        self.get_sample_path()
        self.run_download_files()
        self.run_merge_fastq()
        self.set_db()
        self.end()

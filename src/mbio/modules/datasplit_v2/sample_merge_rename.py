# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190125

import os
import re
import json
import time
import shutil
import subprocess
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
from biocluster.api.file.lib.transfer import MultiFileTransfer


class SampleMergeRenameModule(Module):
    """
    合并fastq文件/多样性样本改名
    """
    def __init__(self, work_id):
        super(SampleMergeRenameModule, self).__init__(work_id)
        options = [
            {"name": "sample_list", "type": "infile", "format": "datasplit.sample_merge_list"},
            {"name": "operation_type", "type": "string", "default": "merge"},  # 操作类型
            {"name": "coll_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.raw_r1, self.raw_r2, self.raw75_r1 = [], [], []
        self.clean_r1, self.clean_r2 = [], []
        self.new_paths, self.download_files = {}, {}
        self.product_type, self.file_name, self.sample_name = "", "", ""

    def check_options(self):
        if not self.option("sample_list").is_set:
            raise OptionError("缺少sample_list,请设置")
        if self.option("operation_type") not in ["merge", "rename"]:
            raise OptionError("operation_type只能是merge/rename")

    def get_sample_path(self):
        """
        获取raw和clean对应的r1、r2,若是s3开头换成要下载的对应文件
        """
        with open(self.option("sample_list").prop["path"], "rb") as f:
            lines = f.readlines()
            for i in range(1, len(lines)):
                item = lines[i].strip().split("\t")
                self.product_type = item[2]
                self.file_name = item[3]
                self.sample_name = item[4]
                if self.product_type != "mirna":  # meta合并raw和clean,mirna只合并clean,其他合并raw
                    raw_path = self.get_new_files(item[5], i, "raw")
                    self.raw_r1.append(raw_path[0])
                    if len(raw_path) == 2:
                        self.raw_r2.append(raw_path[1])
                if self.product_type == "mirna" or self.product_type == "meta":
                    clean_path = self.get_new_files(item[6], i, "clean")
                    self.clean_r1.append(clean_path[0])
                    if len(clean_path) == 2:
                        self.clean_r2.append(clean_path[1])
                    if len(item) > 7:
                        raw75_path = self.get_new_files(item[7], i, "raw75")
                        self.raw75_r1.append(raw75_path[0])

    def get_new_files(self, paths, i, type):
        new_paths_key = []
        old_paths = paths.split(";")
        for j in range(len(old_paths)):
            new_paths_key.append(type+"_"+str(i)+"_"+str(j))
            self.new_paths[type+"_"+str(i)+"_"+str(j)] = old_paths[j]
            if old_paths[j].startswith("s3"):
                self.download_files[type+"_"+str(i)+"_"+str(j)] = old_paths[j]
        return new_paths_key

    # def run_fastq_download(self):
    #     tmp_dir = os.path.join(self.work_dir, "tmp_dir")
    #     if os.path.exists(tmp_dir):
    #         shutil.rmtree(tmp_dir)
    #     os.mkdir(tmp_dir)
    #     # transfer = TransferManager()
    #     # i = 0
    #     # for key in self.download_files.keys():
    #     #     from_path = self.download_files[key]
    #     #     to_path = os.path.join(tmp_dir, os.path.basename(from_path))
    #     #     if os.path.exists(to_path):
    #     #         i += 1
    #     #         tmp_dir1 = os.path.join(tmp_dir, str(i))
    #     #         os.mkdir(tmp_dir1)
    #     #         to_path = os.path.join(tmp_dir1, os.path.basename(from_path))
    #     #     self.new_paths[key] = to_path
    #     #     try:
    #     #         transfer.add(from_uri=from_path, to_uri=to_path)
    #     #     except Exception, err:
    #     #         self.logger.info("-------------------")
    #     #         self.logger.info(err)
    #     #         if "mirna_clean/fasta" in from_path:
    #     #             from_path = from_path.replace("mirna_clean/fasta/", "")
    #     #             transfer.add(from_uri=from_path, to_uri=to_path)
    #     #     self.logger.info(from_path)
    #     #     self.logger.info(to_path)
    #     # transfer.wait()
    #     transfer = MultiFileTransfer()
    #     i = 0
    #     to_dir = tmp_dir
    #     for key in self.download_files.keys():
    #         from_path = self.download_files[key]
    #         to_path = os.path.join(tmp_dir, os.path.basename(from_path))
    #         if os.path.exists(to_path):
    #             i += 1
    #             tmp_dir1 = os.path.join(tmp_dir, str(i))
    #             to_dir = tmp_dir1
    #             os.mkdir(tmp_dir1)
    #             to_path = os.path.join(tmp_dir1, os.path.basename(from_path))
    #         self.new_paths[key] = to_path
    #         try:
    #             transfer.add_download(from_path=from_path, to_path=to_dir)
    #         except Exception, err:
    #             self.logger.info("-------------------")
    #             self.logger.info(err)
    #             if "mirna_clean/fasta" in from_path:
    #                 from_path = from_path.replace("mirna_clean/fasta/", "")
    #                 transfer.add_download(from_path=from_path, to_path=to_dir)
    #         self.logger.info(from_path)
    #         self.logger.info(to_path)
    #     # transfer.wait()
    #     time.sleep(30)
    #     if self.option("operation_type") == "merge":
    #         self.run_merge_fastq()
    #     elif self.option("operation_type") == "rename":
    #         self.run_rename_fastq()

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
            if self.option("operation_type") == "merge":
                self.run_merge_fastq()
            elif self.option("operation_type") == "rename":
                self.run_rename_fastq()
        elif len(self.tools) == 1:
            if self.option("operation_type") == "merge":
                self.tools[0].on("end", self.run_merge_fastq)
            elif self.option("operation_type") == "rename":
                self.tools[0].on("end", self.run_rename_fastq)
            self.tools[0].run()
        else:
            if self.option("operation_type") == "merge":
                self.on_rely(self.tools, self.run_merge_fastq)
            elif self.option("operation_type") == "rename":
                self.on_rely(self.tools, self.run_rename_fastq)
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
        if self.raw_r1:
            self.merge_fastq(self.raw_r1, self.file_name + ".R1.raw.fastq.gz")
        if self.raw_r2:
            self.merge_fastq(self.raw_r2, self.file_name + ".R2.raw.fastq.gz")
        if self.product_type == "meta":
            if self.clean_r2:
                self.merge_meta_fastq(self.clean_r1, self.file_name + ".clean.1.fastq", self.sample_name)
                self.merge_meta_fastq(self.clean_r2, self.file_name + ".clean.2.fastq", self.sample_name)
            else:
                self.merge_meta_fastq(self.clean_r1, self.file_name + ".clean.fastq", self.sample_name)
        elif self.product_type == "mirna":
            if self.raw75_r1:
                self.merge_fastq(self.raw75_r1, self.file_name + ".fq.gz")
            self.merge_fastq(self.clean_r1, self.file_name + ".fasta.gz")
        # else:
        #     self.merge_fastq(self.clean_r1, self.file_name + ".clean.1.fastq.gz")
        #     if self.clean_r2:
        #         self.merge_fastq(self.clean_r2, self.file_name + ".clean.2.fastq.gz")
        self.run_md5sum()
        # self.set_db()

    def run_rename_fastq(self):
        """
        多样性样本重命名
        """
        if self.product_type == "meta":
            if self.clean_r2:
                self.merge_meta_fastq(self.clean_r1, self.file_name + ".clean.1.fastq", self.sample_name)
                self.merge_meta_fastq(self.clean_r2, self.file_name + ".clean.2.fastq", self.sample_name)
            else:
                self.merge_meta_fastq(self.clean_r1, self.file_name + ".clean.fastq", self.sample_name)
        # self.set_db()
        self.run_md5sum()

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

    def merge_meta_fastq(self, key_list, file_name, sample_name):
        """
        多样性fastq文件合并和改名
        """
        path_list = []
        for key in key_list:
            path_list.append(self.new_paths[key])
        out_file = os.path.join(self.work_dir, "merge."+file_name)
        i = 0
        with open(out_file, "wb") as w:
            for f in path_list:
                new_f = f
                if f.endswith(".gz"):
                    new_f = os.path.join(self.work_dir, os.path.basename(f).split(".gz")[0])
                    os.system("gunzip -c {} > {}".format(f, new_f))
                with open(new_f) as fastq:
                    for line in fastq:
                        name = line.split("\t")  # 将序列ID重命名
                        seq_id = "@" + sample_name + "_" + str(i) + "\t" + name[1]
                        w.write('{}{}{}{}'.format(seq_id, next(fastq), next(fastq), next(fastq)))
                        i = i + 1
        out_file1 = os.path.join(self.output_dir, file_name + ".gz")
        os.system("gzip -c {} > {}".format(out_file, out_file1))
        self.logger.info("%s cat完成" % " ".join(path_list))

    def set_db(self):
        """
        结果导表
        """
        if self.option("coll_id"):
            self.logger.info("开始进行结果导表")
            datasplit_api = self.api.api("datasplit.datasplit_merge")
            data_json = os.path.dirname(self.work_dir) + "/data.json"
            s3_output_dir = json.loads(open(data_json).read())["output"]
            datasplit_api.update_sg_split_specimen_merge(self.option("coll_id"), self.output_dir, s3_output_dir, self.option("operation_type"))
        self.end()
        

    def run_md5sum(self):
        self.logger.info("开始进行md5校验")
        # self.md_tool = []
        # for f in os.listdir(self.output_dir):
            # if f in ["meta_qc", "dna_raw"]:
        options = {"fastq_dir": os.path.join(self.output_dir)}
        self.md5sum = self.add_tool("datasplit_v2.md5sum")
        self.md5sum.set_options(options)
        # self.md_tool.append(self.md5sum)
        # self.logger.info(len(self.md_tool))
        # if len(self.md_tool) > 1:
        #     self.on_rely(self.md_tool, self.set_db)
        #     for md_tool in self.md_tool:
        #         md_tool.run()
        # elif len(self.md_tool) == 1:
        self.md5sum.on("end", self.set_db)
        self.md5sum.run()
        # else:
        #     self.end()

    def run(self):
        super(SampleMergeRenameModule, self).run()
        self.get_sample_path()
        self.run_fastq_download()

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SampleMergeRenameModule, self).end()

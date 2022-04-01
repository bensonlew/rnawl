# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20210617

"""提取蛋白代谢需要上传的文件"""
import os
import re
import json
import shutil
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.api.file.lib.transfer import TransferManager


class UploadProteinMetabolomeAgent(Agent):
    def __init__(self, parent):
        super(UploadProteinMetabolomeAgent, self).__init__(parent)
        options = [
            {"name": "project_sn", "type": "string"},  # 项目的合同号
            {"name": "task_sn", "type": "string"},  # 项目的任务号
            {"name": "samples", "type": "string"},  # 样本列表，以分号分隔
            {"name": "is_upload", "type": "bool", "default": False},  # 是否上传文件到对象存储
        ]
        self.add_option(options)
        self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option("project_sn"):
            raise OptionError("必须设置项目的合同号")
        if not self.option("samples"):
            raise OptionError("必须设置项目的样本列表")

    def set_resource(self):
        self._cpu = "1"
        self._memory = "10G"

    def end(self):
        # self.add_upload_dir(self.output_dir)
        super(UploadProteinMetabolomeAgent, self).end()


class UploadProteinMetabolomeTool(Tool):
    def __init__(self, config):
        super(UploadProteinMetabolomeTool, self).__init__(config)
        self._version = 1.0
        self.python = "program/Python/bin/python"
        self.find_sample_rawdatas = self.config.PACKAGE_DIR + "/datasplit/find_sample_rawdatas.py"

    def get_sample_rawdata(self):
        """
        获取项目的rawdata path
        """
        out_table = os.path.join(self.work_dir, self.option("project_sn") + ".path.xls")
        out_dir = os.path.join(self.work_dir, self.option("project_sn") + ".output_dir.txt")
        cmd = "{} {} -MJ {} -samples {} -out {} -out_dir {}".format(self.python, self.find_sample_rawdatas,
              self.option("project_sn"), self.option("samples"), out_table, out_dir)
        if self.option("task_sn"):
            cmd += " -PM {}".format(self.option("task_sn"))
        command = self.add_command("find_sample_rawdatas", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("find_sample_rawdatas运行成功")
        else:
            self.logger.info("find_sample_rawdatas运行失败")
            self.set_db("falied", "find_sample_rawdatas运行失败")
            self.set_error("find_sample_rawdatas运行失败")

    def new_split_path(self, split_path):
        if os.path.exists(split_path):
            return split_path
        if "ilustre" in split_path:
            split_path1 = split_path.replace("ilustre", "clustre")
            if os.path.exists(split_path1):
                return split_path1
        if "sglustre" in split_path:
            split_path1 = split_path.replace("sglustre", "ilustre")
            if os.path.exists(split_path1):
                return split_path1
        return split_path

    def upload_output_dir(self):
        """
        上传文件
        """
        out_dir = os.path.join(self.work_dir, self.option("project_sn") + ".output_dir.txt")
        with open(out_dir, "rb") as f:
            for line in f:
                item = line.strip().split("\t")
                self.old_dir = item[0]
        # new_old_dir = self.old_dir.replace("ilustre", "clustre")
        new_old_dir = self.new_split_path(self.old_dir)
        if not os.path.exists(new_old_dir):
            self.set_db("falied", "文件夹:%s 没有找到，请检查" % new_old_dir)
            self.set_error("文件夹:%s 没有找到，请检查" % new_old_dir)
        data_json = os.path.dirname(self.work_dir) + "/data.json"
        s3_output_dir = json.loads(open(data_json).read())["output"]
        self.s3_dir = os.path.join(s3_output_dir, self.option("project_sn"))
        if self.option("task_sn") and new_old_dir.endswith(self.option("task_sn")):
            self.s3_dir = os.path.join(self.s3_dir, self.option("task_sn"))
        self.transfer = TransferManager()
        self.upload_file(new_old_dir, self.s3_dir)
        self.transfer.wait()

    def upload_file(self, old_dir, new_dir):
        for f in os.listdir(old_dir):
            old_path = os.path.join(old_dir, f)
            new_path = os.path.join(new_dir, f)
            if os.path.isdir(old_path):
                if f == "result" or f == "Result":  # 搜库文件夹路径
                    self.search_lib_dir[old_path] = new_path + "/"
                elif "QC" in f:
                    m = re.match(".*(QC\d+).*", f)
                    if m:
                        sample_name = m.group(1)
                        if sample_name not in self.qc_path.keys():
                            self.qc_path[sample_name] = {"path": [], "s3_path": []}
                        self.qc_path[sample_name]["path"].append(old_path)
                        self.qc_path[sample_name]["s3_path"].append(new_path + "/")
                self.upload_file(old_path, new_path)
                self.new_s3_path[old_path] = new_path + "/"
            else:
                if "QC" in f:
                    m = re.match(".*(QC\d+).*", f)
                    if m:
                        sample_name = m.group(1)
                        if sample_name not in self.qc_path.keys():
                            self.qc_path[sample_name] = {"path": [], "s3_path": []}
                        self.qc_path[sample_name]["path"].append(old_path)
                        self.qc_path[sample_name]["s3_path"].append(new_path)
                self.new_s3_path[old_path] = new_path
                self.transfer.add(from_uri=old_path, to_uri=new_path)

    def create_sample_path(self):
        """
        创建项目的样本对应的path
        """
        out_table = os.path.join(self.work_dir, self.option("project_sn") + ".path.xls")
        self.new_out_table = os.path.join(self.output_dir, self.option("project_sn") + ".sample.path.xls")
        with open(out_table, "rb") as f, open(self.new_out_table, "wb") as w:
            lines = f.readlines()
            w.write(lines[0].split("\n")[0] + "\ts3_paths\tsearch_lib_dir\ts3_search_lib_dir\n")
            search_lib_dir, s3_search_lib_dir = "", ""
            if len(self.search_lib_dir.keys()) > 0:
                search_lib_dir = self.search_lib_dir.keys()[0]
                s3_search_lib_dir = self.search_lib_dir[search_lib_dir]
            for line in lines[1:]:
                item = line.strip().split("\t")
                new_paths, s3_paths = [], []
                if len(item) > 1:
                    paths = item[1].split(";")
                    for f in paths:
                        # new_paths.append(f.replace("ilustre", "clustre"))
                        new_f = self.new_split_path(f)
                        new_paths.append(new_f)
                        if f in self.new_s3_path.keys():
                            s3_paths.append(self.new_s3_path[f])
                        else:
                            s3_paths.append("-")
                w.write(item[0] + "\t" + ";".join(new_paths) + "\t" + ";".join(s3_paths) + "\t")
                w.write(search_lib_dir + "\t" + s3_search_lib_dir + "\n")

    def create_qc_path(self):
        """
        创建qc对应的path
        """
        self.qc_table = os.path.join(self.output_dir, self.option("project_sn") + ".qc.path.xls")
        with open(self.qc_table, "wb") as w:
            w.write("project_sn\ttask_sn\tsample_name\tqc_path\ts3_qc_path\n")
            for s in self.qc_path.keys():
                w.write(self.option("project_sn") + "\t" + self.option("task_sn") + "\t" + s + "\t")
                w.write(";".join(self.qc_path[s]["path"]) + "\t" + ";".join(self.qc_path[s]["s3_path"]) + "\n")

    def set_db(self, status, desc=""):
        if self.option("is_upload"):
            api_db = datasplit_api = self.api.api("datasplit.datasplit_upload")
            if status != "end":
                api_db.update_coll_status_protein(status=status, project_sn=self.option("project_sn"), task_sn=self.option("task_sn"), samples=self.option("samples"), desc=desc)
            else:
                api_db.update_coll_path_protein(project_sn=self.option("project_sn"), task_sn=self.option("task_sn"), sample_path=self.new_out_table, work_dir=self.old_dir, s3_dir=self.s3_dir)
                api_db.update_coll_qc_path_protein(project_sn=self.option("project_sn"), task_sn=self.option("task_sn"), qc_table=self.qc_table)

    def run(self):
        super(UploadProteinMetabolomeTool, self).run()
        self.old_dir, self.s3_dir = "", ""
        self.search_lib_dir, self.qc_path = {}, {}
        self.set_db("start")
        self.get_sample_rawdata()
        self.new_s3_path = {}
        # self.link_output_dir()
        if self.option("is_upload"):
            self.upload_output_dir()
        self.create_sample_path()
        self.create_qc_path()
        self.set_db("end")
        self.end()

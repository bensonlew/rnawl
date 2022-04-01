# -*- coding: utf-8 -*-
# __author__: zengjing
# last_modify: 20190314

import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class UploadDatasplitAgent(Agent):
    """
    上传线下拆分的文库结果、原型样本结果、质控结果，将fastq上传到对象存储，统计结果上传到对应mongo数据库
    """
    def __init__(self, parent=None):
        super(UploadDatasplitAgent, self).__init__(parent)
        options = [
            {"name": "upload_file", "type": "infile", "format": "datasplit.upload_file", "required": True},  # 上传的文件
            {"name": "seq_board", "type": "string", "required": True},  # 测序板
            {"name": "upload_type": "type": "string", "required": True},  # 上传的类型,library、raw_sample、clean_sample
            {"name": "upload_id", "type": "string", "required": True},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option("upload_file").is_set:
            raise OptionError("请设置需要上传的文件")
        if self.option("upload_type") not in ["library", "raw_sample", "clean_sample"]:
            raise OptionError("上传的类型:%s只能是library/raw_sample/clesn_sample" % self.option("upload_type"))

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(UploadDatasplitAgent, self).end()


class UploadDatasplitTool(Tool):
    def __init__(self, config):
        super(UploadDatasplitTool, self).__init__(config)

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

    def library_set_output(self):
        """
        将upload_file里的文件link到output里
        """
        info = self.upload_api.check_lib_info(self.option("seq_board"), self.option("upload_type").prop["path"])

        with open(self.option("upload_type").prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                lane_name = item[2]
                lib_name = item[0]
                lib_dir = os.path.join(self.output_dir, lane_name + "." + lib_name)
                if os.path.exists(lib_dir):
                    self.set_error("文库：%s重复了，请检查" % lib_name)
                os.mkdir(lib_dir)
                for f in item[16].split(";"):
                    # self.link_file(f.replace("ilustre", "clustre"), os.path.join(lib_dir, os.path.basename(f)))
                    f_ = self.new_split_path(f)
                    self.link_file(f_, os.path.join(lib_dir, os.path.basename(f)))

    def specimen_set_output(self):
        """
        将upload_file里的文件link到output里
        """
        with open(self.option("upload_type").prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                lane_name = item[2]
                lib_name = item[0]
                lib_dir = os.path.join(self.output_dir, lane_name + "." + lib_name)
                if os.path.exists(lib_dir):
                    self.set_error("文库：%s重复了，请检查" % lib_name)
                os.mkdir(lib_dir)
                for f in item[16].split(";"):
                    self.link_file(f.replace("ilustre", "clustre"), os.path.join(lib_dir, os.path.basename(f)))

    def link_file(self, old, new):
        """
        link 文件
        """
        if os.path.isdir(old):
            if not os.path.exists(new):
                os.mkdir(new)
            for f in os.listdir(old):
                old_ = os.path.join(old, f)
                new_ = os.path.join(new, f)
                self.link_file(old_, new_)
        else:
            if os.path.exists(new):
                os.remove(new_)
            os.link(old, new)

    def run(self):
        super(UploadDatasplitTool, self).run()
        self.upload_api = self.api.api("datasplit.upload_datasplit")
        if self.option("upload_type") == "library":
            self.library_set_output()
        else:
            self.specimen_set_output()
        self.end()

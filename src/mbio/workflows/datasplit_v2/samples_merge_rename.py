# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20201105

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
# from upload_s5cmd import UploadS5cmd
from submit import Submit
import gevent
import json


class SamplesMergeRenameWorkflow(Workflow):
    """
    合并相同样本的fastq文件/样本重命名
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SamplesMergeRenameWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "sample_list", "type": "string"},
            {"name": "fx_id", "type": "string"},
            {"name": "operation_type", "type": "string", "default": "merge"}  # 操作类型
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("sample_list"):
            raise OptionError("缺少sample_list,请设置")
        if not self.option("fx_id"):
            raise OptionError("缺少fx_id,请设置")
        if self.option("operation_type") not in ["merge", "rename"]:
            raise OptionError("operation_type只能是merge/rename")

    def get_sample_info(self):
        with open(self.option("sample_list"), "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                if self.option("operation_type") == "merge":
                    self.post_sample_merge(item[1])
                elif self.option("operation_type") == "rename":
                    self.post_sample_rename(item[1])

    def post_sample_merge(self, merge_id):
        """
        启动接口，进行样本合并
        """
        merge_params = {"fx_id": self.option("fx_id"), "coll_id": merge_id, "operation_type": "merge"}
        # if self._sheet.client01:
        #     type = "tsanger"
        # else:
        #     type = "tsg"
        type = "tsanger"
        self.logger.info("开始投递合并样本的接口:%s" % merge_id)
        result = Submit(merge_params, "/s/datasplit/datasplit_sample_merge_rename", type).webapitest()
        self.logger.info(result)
        try:
            result = json.loads(result)
        except:
            self.logger.info("result不是一个json")
        datasplit_api = self.api.api("datasplit.datasplit_new")
        try:
            if result["success"]:
                datasplit_api.update_sg_split_specimen_status(merge_id, "sg_split_specimen_merge", "start", "开始运行")
                self.logger.info("样本合并:%s运行成功！" % merge_id)
            else:
                datasplit_api.update_sg_split_specimen_status(merge_id, "sg_split_specimen_merge", "failed", result["info"])
                self.logger.info("样本合并:%s运行失败！" % merge_id)
        except:
            datasplit_api.update_sg_split_specimen_status(merge_id, "sg_split_specimen_merge", "failed", "运行失败")
            self.logger.info("样本合并:%s运行失败！" % merge_id)

    def post_sample_rename(self, rename_id):
        """
        启动接口，进行样本重命名
        """
        rename_params = {"fx_id": self.option("fx_id"), "coll_id": rename_id, "operation_type": "rename"}
        # if self._sheet.client01:
        #     type = "tsanger"
        # else:
        #     type = "tsg"
        type = "tsanger"
        self.logger.info("开始投递重命名样本的接口:%s" % rename_id)
        result = Submit(rename_params, "/s/datasplit/datasplit_sample_merge_rename", type).webapitest()
        self.logger.info(result)
        try:
            result = json.loads(result)
        except:
            self.logger.info("result不是一个json")
        datasplit_api = self.api.api("datasplit.datasplit_new")
        try:
            if result["success"]:
                datasplit_api.update_sg_split_specimen_status(rename_id, "sg_split_specimen_rename", "start", "开始运行")
                self.logger.info("样本重命名:%s运行成功！" % rename_id)
            else:
                datasplit_api.update_sg_split_specimen_status(rename_id, "sg_split_specimen_rename", "failed", result["info"])
                self.logger.info("样本重命名:%s运行失败！" % rename_id)
        except:
            datasplit_api.update_sg_split_specimen_status(rename_id, "sg_split_specimen_rename", "failed", "运行失败")
            self.logger.info("样本重命名:%s运行失败！" % rename_id)

    def run(self):
        self.start_listener()
        self.fire("start")
        self.get_sample_info()
        gevent.spawn_later(5, self.end)
        super(SamplesMergeRenameWorkflow, self).run()

    def end(self):
        super(SamplesMergeRenameWorkflow, self).end()

# !usr/bin/python
# -*- coding: utf-8 -*-
from collections import OrderedDict
from bson import SON
from biocluster.config import Config
import types
import re
import datetime
import os
import json
from api_base import ApiBase


class Manhattan(ApiBase):
    """..."""
    def __init__(self, bind_object):
        super(Manhattan, self).__init__(bind_object)

    def sg_manhattan(self, task_id, project_sn, params):  # 用于测试建立主表
        """
        曼哈顿图主表
        :return:
        """
        main_id = self.add_main_table("sg_manhattan", task_id, project_sn, params, "origin_sg_manhattan",
                                      "曼哈顿图导表", "开始进行曼哈顿图导表")
        self.update_db_record("sg_manhattan", {"_id": main_id}, {"main_id": main_id, "status": "end"})
        return main_id

    def add_sg_manhattan_path(self, main_id, path):
        id = self.check_objectid(main_id)
        self.update_db_record("sg_manhattan", {"_id": id}, {"manhattan_path": path})





if __name__ == "__main__":
    a = Manhattan(None)
    task_id = "tool_lab"
    project_sn = "tool_lab_20200416"
    # a.sg_manhattan(task_id, project_sn, "manhattan")
    a.add_sg_manhattan_path("5e97eeb517b2bf6d8f8f35f8", "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/cloud_tools/qqman/output1/qqman.png")


# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 20190115

from api_base import ApiBase
from collections import defaultdict
import os
import datetime
import json


class TagDetail(ApiBase):
    """
    无参WGS导表：软件列表、样本信息、质控
    """
    def __init__(self, bind_object):
        super(TagDetail, self).__init__(bind_object)
        self._project_type = "dna_noref_wgs"

    def check_exists(self, file):
        """
        检查file文件是否存在
        """
        if not os.path.exists(file):
            raise Exception("文件：%s不存在，请检查" % file)
        return True

    def add_sg_tag_detail(self, main_id, tag_id, ref, snp):
        """
        基因结构图
        :return:
        """
        origin_id = self.check_objectid(main_id)
        data_list = []
        insert_data = {
            "origin_id": origin_id,
            "tag_id": tag_id,
            "ref": ref,
            "snp": snp
        }
        data_list.append(insert_data)
        self.col_insert_data("sg_tag_detail", data_list)


if __name__ == "__main__":
    a = TagDetail(None)
    task_id = "noref_test1"
    infile = "/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/wucan/06.genotype/test_module_data/test_output/test.tag"
    a.add_sg_structure(infile, task_id)

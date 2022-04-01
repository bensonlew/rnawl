# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from collections import OrderedDict
from bson import SON
from biocluster.config import Config
import gzip
import re
import datetime
import os
import json
from api_base import ApiBase

class VcfDistribution(ApiBase):
    def __init__(self, bind_object):
        super(VcfDistribution, self).__init__(bind_object)
        self._db_name = "sanger_tool_lab"

    def add_detail(self, main_id, picture_path):
        """
        密度分布图图片路径
        """
        insert_dict = {
            "picture":picture_path
         }
        self.update_db_record(
            "vcf_distribution",{"_id":self.check_objectid(main_id)},insert_dict
        )
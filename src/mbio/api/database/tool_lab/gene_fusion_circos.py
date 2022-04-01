# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20210714
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import datetime
import json
import os
import pandas as pd
import glob
import re
from mbio.api.database.ref_rna_v2.api_base import ApiBase

class GeneFusionCircos(ApiBase):
    def __init__(self, bind_object):
        super(GeneFusionCircos, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_s3_result(self, main_id, s3_output):
        # add main table info
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        self.update_db_record('gene_fusion_circos', main_id, status='end', main_id=main_id, s3_output = s3_output)




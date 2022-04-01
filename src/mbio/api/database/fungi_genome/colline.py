# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modify:20180622

from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId
import re


class Colline(Base):
    def __init__(self, bind_object):
        super(Colline, self).__init__(bind_object)
        self._project_type = "fungigenome"

    @report_check
    def add_detail(self, detail_file_path, main_id=None,  specimen_id=None, ref=None,data_path=None ):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        data_list = []
        collection_detail = self.db['colline_detail']
        detail_file = os.path.join(detail_file_path,"{}_block_new.xls".format(specimen_id))

        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                data = [("circos_id", main_id)]
                line = line.strip().split('\t')
                line[0]=re.sub(' ','',line[0])
                if line[6] == 'seq1':
                    line[6]=specimen_id
                else:
                    line[6]=ref

                data.extend(
                    [("block", line[0]),("specimen_id",line[6]),("strand",line[2]),("start",int(line[3])),
                     ("end",line[4]),("len",line[5]),("location",line[1][2:]),("loc_start", line[7]),('loc_end',line[8])])

                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection_detail.insert_many(data_list)
            main_collection = self.db['colline']
            main_collection.update({"_id": main_id}, {"$set": {"status": 'end',"origin_id": main_id, "cir_path":data_path}})

        except Exception as e:
            self.bind_object.logger.error("导入colline_detail %s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入colline_detail %s信息出错:%s" , variables=(detail_file, e), code="52101901")
        else:
            self.bind_object.logger.info("导入colline_detail%s信息成功!" % detail_file)


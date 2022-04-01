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


class Pseudogene(Base):
    def __init__(self, bind_object):
        super(Pseudogene, self).__init__(bind_object)
        self._project_type = "fungigenome"

    @report_check
    def add_detail(self, detail_file_path, main_id=None,  specimen_id=None):
        if not isinstance(main_id, ObjectId):
             main_id = ObjectId(main_id)
        data_list = []
        collection_detail = self.db['pseudogene_detail']
        detail_file = os.path.join(detail_file_path,"{}_pseudogene.xls".format(specimen_id))

        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                data = [("pse_id", main_id)]
                line = line.strip().split('\t')

                scf = line[2].capitalize()

                data.extend(
                    [("name", line[0]),("ref",line[1]),("location",scf),("start",line[3]),
                     ("end",line[4]),("shift",line[5]),("stop",line[6]),("specimen_id", specimen_id)])

                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection_detail.insert_many(data_list)
            main_collection = self.db['pseudogene']
            main_collection.update({"_id": main_id}, {"$set": {"status": 'end',"origin_id": main_id}})

        except Exception as e:
            self.bind_object.logger.error("导入pseudogene_detail %s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入pseudogene_detail %s信息出错:%s" , variables=(detail_file, e), code="52102501")
        else:
            self.bind_object.logger.info("导入pseudogene_detail%s信息成功!" % detail_file)

        collection_sum = self.db['pseudogene_sum']

        data = {
            "pse_id": main_id,
            "specimen_id": specimen_id,
            "num": str(len(lines)-1)
        }
        try:
            collection_sum.insert_one(data)
        except Exception:
            self.bind_object.logger.error('导入pseudogene_sum信息出错')
            self.bind_object.set_error('导入pseudogene_sum信息出错', code="52102502")
        else:
            self.bind_object.logger.info('导入pseudogene_sum信息成功')

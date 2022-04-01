# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'

from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson import ObjectId, SON
from bson.objectid import ObjectId
from types import StringTypes
# from biocluster.config import Config


class MetaSourcetracker(Base):
    """
    微生物来源组成比例分析
    """
    def __init__(self, bind_object):
        super(MetaSourcetracker, self).__init__(bind_object) #
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB
        self.task_id = ""

    @report_check
    def add_sg_sourcetracker_detail(self, id=None, file_path=None, stdev_file_path=None, name_1=None, name_2=None):
        insert_data = list()
        d_dict = dict()
        def addtwodimdict(thedict, key_a, key_b, val):
            if key_a in thedict:
                thedict[key_a].update({key_b: val})
            else:
                thedict.update({key_a: {key_b: val}})
        with open(stdev_file_path, "rb") as s:
            head_s = s.next().strip('\r\n')
            head_s = re.split('\t', head_s)
            source_name = head_s[1:]
            for line in s:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                stdev_number = line[1:]
                for i in range(0, len(source_name)):
                    addtwodimdict(d_dict, line[0], source_name[i], float(stdev_number[i]))
        print(d_dict)
        with open(file_path, 'rb') as r:
            head = r.next().strip('\r\n')  # windows换行符
            head = re.split('\t', head)
            new_head = head[1:]
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                group_num = line[1:]
                detail = dict()
                detail['sourcetracker_id'] = ObjectId(id)
                detail['file_name_1'] = name_1
                detail['file_name_2'] = name_2
                detail['sample_name'] = line[0]
                for i in range(0, len(group_num)):
                    # detail[new_head[i]] = float(group_num[i])
                    detail[new_head[i] + "_mean"] = float(group_num[i])
                    detail[new_head[i] + "_stdev"] = float(d_dict[line[0]][new_head[i]])
                    detail[new_head[i] + "_h"] = float(group_num[i]) + (d_dict[line[0]][new_head[i]])
                    detail[new_head[i] + "_l"] = float(group_num[i]) - (d_dict[line[0]][new_head[i]])
                    # detail[new_head[i] + "_stdev_plot"] = float(d_dict[line[0]][new_head[i]]*d_dict[line[0]][new_head[i]])
                    # detail[new_head[i] + "_h"] = float(group_num[i]) + (d_dict[line[0]][new_head[i]]*d_dict[line[0]][new_head[i]])
                    # detail[new_head[i] + "_l"] = float(group_num[i]) - (d_dict[line[0]][new_head[i]]*d_dict[line[0]][new_head[i]])
                insert_data.append(detail)
        try:
            collection = self.db['sg_sourcetracker_detail']
            collection.insert_many(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入sg_sourcetracker_detail表格信息出错:{}".format(e))
            self.bind_object.set_error("导入sg_sourcetracker_detail表格信息出错", code="51003801")
        else:
            self.bind_object.logger.info("导入sg_sourcetracker_detail表格成功")
# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'

from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson import ObjectId, SON
from bson.objectid import ObjectId
from types import StringTypes
# from biocluster.config import Config


class EnterotypingDb(Base):
    """
    样本菌群分型分析导表
    """
    def __init__(self, bind_object):
        super(EnterotypingDb, self).__init__(bind_object) #
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB
        self.task_id = ""

    @report_check
    def add_sg_enterotyping(self, params, from_otu_table, name = None, cluster_name = None, spe_name = None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51002601")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("sg_otu表中找不到相应记录", code="51002602")
        project_sn = result['project_sn']
        self.task_id = result['task_id']
        if not name:
            name = "enterotyping_analysis_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") #
        insert_data = {
            # "project_sn": "dd",
            # 'task_id': "ss",
            "project_sn": project_sn,
            'task_id': self.task_id,
            # 'otu_id': "ObjectId(\"" + str(from_otu_table) + "\")",
            'otu_id': from_otu_table,
            'cluster_name': cluster_name,
            'spe_name': spe_name,
            'name': name,
            "params": params,
            'status': 'end',
            'desc': 'result after Enterotyping analysis',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "show": 0,
            "type": "otu_Enterotyping"
        }
        collection = self.db["sg_enterotyping"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        #collection.update({"_id": inserted_id}, {"$set": {"main_id": inserted_id}})
        return inserted_id

    # @report_check
    # def add_sg_enterotyping_ch(self, id = None, file_path = None, table_id = None, group_id = None, from_otu_table = None, levle_id = None, Major = False):
    #     data_list = []
    #     with open(file_path, 'rb') as r:
    #         data = r.readlines()
    #         for f in data:
    #             f = f.strip().split("\t")
    #             data = [("sg_enterotyping_id",id), ("x", f[0]), ("y", f[1])]
    #             data_son = SON(data)
    #             data_list.append(data_son)
    #     try:
    #         collection = self.db["sg_enterotyping_ch"]
    #         collection.insert_many(data_list)
    #     except Exception, e:
    #         print "sg_enterotyping_ch failure{}".format(e)
    #     else:
    #         print "sg_enterotyping_ch sucess"

    @report_check
    def add_sg_enterotyping_detail(self, id = None, file_path=None, x = None, y = None, name = None, detail_name=None):
        data_list = []
        with open(file_path, 'rb') as r:
            if name == "cluster.txt":
                data = r.readlines()[1:]
            else:
                data = r.readlines()
            for f in data:
                f = f.strip().split("\t")
                if len(f) == 3:
                    data = [("enterotyping_id", ObjectId(id)), ("name_id", name), (x, f[0]), (y, f[1]), (detail_name, f[2])]
                else:
                    data = [("enterotyping_id", ObjectId(id)), ("name_id", name), (x, f[0]), (y, f[1])]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["sg_enterotyping_detail"]
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入sg_enterotyping_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入sg_enterotyping_detail表格成功")
        # except Exception, e:
        #     print "sg_enterotyping_detail failure{}".format(e)
        # else:
        #     print "sg_enterotyping_detail sucess"

    @report_check
    def add_sg_enterotyping_detail_cluster(self, id=None, file_path=None, name = None):
        insert_data = list()
        with open(file_path, 'rb') as r:
            head = r.next().strip('\r\n')  # windows换行符
            head = re.split('\t', head)
            new_head = head[1:]
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                sample_num = line[1:]
                classify_list = re.split(r"\s*;\s*", line[0])
                otu_detail = dict()
                otu_detail['enterotyping_id'] = ObjectId(id)
                otu_detail['name_id'] = name
                for cf in classify_list:
                    if cf != "":
                        otu_detail[cf[0:3].lower()] = cf
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = sample_num[i]
                insert_data.append(otu_detail)
        try:
            collection = self.db['sg_enterotyping_detail_cluster']
            collection.insert_many(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入sg_enterotyping_detail_cluster表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入sg_enterotyping_detail_cluster表格成功")

    @report_check
    def add_sg_enterotyping_detail_summary(self, id=None, file_path=None, name=None, cluster_name=None, spe_name=None):
        insert_data = list()
        with open(file_path, 'rb') as r:
            head = r.next().strip('\r\n')  # windows换行符
            head = re.split('\t', head)
            new_head = head[1:]
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                sample_num = line[1:]
                otu_detail = dict()
                otu_detail['enterotyping_group'] = line[0]
                otu_detail['enterotyping_id'] = ObjectId(id)
                otu_detail['name_id'] = name
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = sample_num[i]
                insert_data.append(otu_detail)
        try:
            collection = self.db['sg_enterotyping_detail_summary']
            collection.insert_many(insert_data)
            main_collection = self.db["sg_enterotyping"]
            self.bind_object.logger.info("开始刷新主表写物种名称和cluster.txt的名称")
            main_collection.update({"_id": ObjectId(id)},
                                   {"$set": {"cluster_name": cluster_name,
                                       "spe_name": spe_name}})
        except Exception as e:
            self.bind_object.logger.error("导入sg_enterotyping_detail_summary表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入sg_enterotyping_detail_summary表格成功")









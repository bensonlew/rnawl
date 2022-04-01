# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
#20190902

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re, os
from types import StringTypes
from biocluster.config import Config
from bson.son import SON


class Nmds(Base):
    def __init__(self, bind_object):
        super(Nmds, self).__init__(bind_object)
        self._project_type = "bac_comparative"
        #self.id = 'tsg_123'
        #self.project_sn = '188_5b5acb3018'

    @report_check
    def add_nmds(self, params=None, name=None):
        """
        pan_genomes的主表；params中记录了group_detail等字段
        :param params:
        :return:
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "nmds分析主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "Nmds_Origin",
        }
        collection = self.db["nmds"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    def add_nmds_data(self, inserted_id, nmds_data=None):
        """
        导表的cluster表
        :param inserted_id: 主表id
        :param pca_data: nmds的dir文件夹
        :return:
        """
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        if nmds_data == None:
            self.bind_object.set_error("没有生成正确的文件夹nmds_data: %s" % (nmds_data))

        site_path = os.path.join(nmds_data, 'nmds_sites.xls')
        stress_path = os.path.join(nmds_data, 'nmds_stress.xls')
        ellipse_path = os.path.join(nmds_data, 'ellipse.xls')
        try:
            self.insert_table_detail(new_inserted_id, stress_path, 'stress')
            self.insert_table_detail(new_inserted_id, site_path, 'specimen')
        except Exception, e:
            self.bind_object.set_error("更新nmds_detail出错:%s" % (e))
        try:
            self.insert_table_detail(new_inserted_id, ellipse_path, 'circle')
        except Exception, e:
            self.bind_object.logger.info("更新nmds_detail出错:%s" % (e))

    def insert_table_detail(self, inserted_id, table_path, type):
        """
        插入详情表
        :param inserted_id: pca主表
        :param table_path: 插入主表信息
        :return:
        """
        if type in ["stress"]:
            stress = ''
            with open(table_path, 'r') as f:
                lines = f.readlines()
                for line in lines[1:]:
                    line = line.strip().split("\t")
                    stress = line[0]

            main_collection = self.db["nmds"]
            try:
                main_collection.update_one({'_id': inserted_id}, {'$set': {'stress': stress}})
            except Exception, e:
                self.bind_object.set_error("更新nmds主表出错：%s" % (e))
        elif type in ["specimen"]:
            with open(table_path, "r") as f:
                lines = f.readlines()
                data_list = []
                specimen_list = "NMDS1,NMDS2"
                for line in lines[1:]:
                    line = line.strip().split("\t")
                    sample_name = line[0]
                    nmds1 = float(line[1])
                    nmds2 = float(line[2])
                    data = {
                        "nmds_id": inserted_id,
                        "specimen_id": sample_name,
                        'NMDS1': nmds1,
                        'NMDS2': nmds2,
                        "type": type
                        }
                    data_son = SON(data)
                    data_list.append(data_son)
            try:
                collection = self.db["nmds_detail"]
                collection.insert_many(data_list)
                self.bind_object.logger.info("导入nmds_detail详情表成功！")
            except Exception, e:
                self.bind_object.set_error("导入nmds_detail结果表出错:%s" % (e))
            try:
                main_collection = self.db["nmds"]
                main_collection.update_one({'_id': inserted_id}, {'$set': {'specimen': specimen_list, 'main_id':inserted_id}})
            except Exception, e:
                self.bind_object.set_error("导入nmds_detail结果表出错:%s" % (e))
        else:
            sample_dict = {}
            with open(table_path, "r") as f:
                lines = f.readlines()
                data_list = []
                header = lines[0].strip().split("\t")
                sample_list = header[1:]
                for line in lines[1:]:
                    line_dict = {}
                    line = line.strip().split("\t")
                    sample_name = line[0]
                    for sp in sample_list:
                        sp_index = sample_list.index(sp)
                        line_dict[sp] = str(line[sp_index+1])
                    sample_dict[sample_name] = line_dict

            for key in sample_dict.keys():
                data = {
                    "nmds_id": inserted_id,
                    "specimen_id": key,
                    "type": type
                    }
                function_dict = sample_dict[key]
                for sp in sample_list:
                    if sp in function_dict.keys():
                        data[sp] = function_dict[sp]
                data_son = SON(data)
                data_list.append(data_son)
            try:
                collection = self.db["nmds_detail"]
                collection.insert_many(data_list)
                self.bind_object.logger.info("导入nmds_detail详情表成功！")
                main_collection = self.db["nmds"]
                main_collection.update_one({'_id': inserted_id}, {'$set': {'circle': ','.join(sample_list), 'main_id':inserted_id}})

            except Exception, e:
                self.bind_object.set_error("导入nmds_detail结果表出错:%s" % (e))
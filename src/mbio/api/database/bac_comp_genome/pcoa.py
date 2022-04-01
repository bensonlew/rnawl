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


class Pcoa(Base):
    def __init__(self, bind_object):
        super(Pcoa, self).__init__(bind_object)
        self._project_type = "bac_comparative"
        #self.id = 'tsg_123'
        #self.project_sn = '188_5b5acb3018'

    @report_check
    def add_pcoa(self, params=None, name=None):
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
            "desc": "pcoa分析主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "Pcoa_Origin",
        }
        collection = self.db["pcoa"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    def add_pcoa_data(self, inserted_id, pcoa_data=None):
        """
        导表的cluster表
        :param inserted_id: 主表id
        :param pca_data: pca的dir文件夹
        :return:
        """
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        if pcoa_data == None:
            self.bind_object.set_error("没有生成正确的文件夹pcoa_data: %s" % (pcoa_data))

        site_path = os.path.join(pcoa_data, 'pcoa_sites.xls')
        eigenvalues_path = os.path.join(pcoa_data, 'pcoa_eigenvalues.xls')
        eigenvaluespre_path = os.path.join(pcoa_data, 'pcoa_eigenvaluespre.xls')
        ellipse_path = os.path.join(pcoa_data, 'ellipse.xls')

        try:
            self.insert_table_detail(new_inserted_id, eigenvaluespre_path, 'importance')
        except Exception, e:
            self.bind_object.set_error("更新pcoa_detail出错:%s" % (e))
        try:
            self.insert_table_detail(new_inserted_id, site_path, 'specimen')
        except Exception, e:
            self.bind_object.set_error("更新pcoa_detail出错:%s" % (e))

        try:
            self.insert_table_detail(new_inserted_id, eigenvalues_path, 'eigenvalues')
        except Exception, e:
            self.bind_object.set_error("更新pcoa_detail出错:%s" % (e))

        try:
            self.insert_table_detail(new_inserted_id, ellipse_path, 'circle')
        except Exception, e:
            self.bind_object.logger.info("更新pca_detail出错:%s" % (e))

    def insert_table_detail(self, inserted_id, table_path, type):
        """
        插入详情表
        :param inserted_id: pca主表
        :param table_path: 插入主表信息
        :return:
        """
        pcoa_list = []
        pcoa_dict = {}
        data_list = []
        if type in ["importance"]:
            with open(table_path, 'r') as f:
                lines = f.readlines()
                for line in lines[1:]:
                    line = line.strip().split("\t")
                    pcoa_name = line[0]
                    pcoa_value = float(line[1])
                    if pcoa_name not in pcoa_list:
                        pcoa_dict[pcoa_name] = pcoa_value
                        pcoa_list.append(pcoa_name)
            data = {
                "pcoa_id": inserted_id,
                "specimen_id": "Proportion of Variance",
                "type": type
            }
            for key in pcoa_list:
                if key in pcoa_dict.keys():
                    data[key] = pcoa_dict[key]
            data_son = SON(data)
            data_list.append(data_son)
            collection = self.db["pcoa_detail"]
            main_collection = self.db["pcoa"]
            try:
                collection.insert_many(data_list)
                self.bind_object.logger.info("导入%s结果表成功" % ("pcoa_detail"))
            except Exception, e:
                self.bind_object.set_error("导入%s结果表出错:%s" % ("pcoa_detail", e))
            try:
                main_collection.update_one({'_id': inserted_id}, {'$set': {'pcoa_list': ','.join(pcoa_list), 'main_id':inserted_id}})
            except Exception, e:
                self.bind_object.set_error("导入%s结果表出错:%s" % ("pcoa_detail", e))
        elif type in ["specimen", "eigenvalues"]:
            main_collection = self.db["pcoa"]
            main_table = main_collection.find_one({'_id': inserted_id})
            pca_list = main_table['pcoa_list'].split(",")
            sample_dict = {}
            sample_list = []
            with open(table_path, "r") as f:
                lines = f.readlines()
                data_list = []
                if type in ['specimen']:
                    header = lines[0].strip().split("\t")
                    sample_list = header
                    for line in lines[1:]:
                        line_dict = {}
                        line = line.strip().split("\t")
                        sample_name = line[0]
                        for sp in sample_list:
                            sp_index = sample_list.index(sp)
                            line_dict[sp] = float(line[sp_index+1])
                        sample_dict[sample_name] = line_dict
                else:
                    sample_list = []
                    sample_name = "eigenvalues"
                    line_dict = {}
                    for line in lines[1:]:
                        line = line.strip().split("\t")
                        pcoa_name = line[0]
                        pcoa_value = float(line[1])
                        if pcoa_name not in sample_list:
                            sample_list.append(pcoa_name)
                            line_dict[pcoa_name] = pcoa_value
                    sample_dict[sample_name] = line_dict
            for key in sample_dict.keys():
                data = {
                    "pcoa_id": inserted_id,
                    "specimen_id": key,
                    "type": type
                    }
                function_dict = sample_dict[key]
                for sp in pca_list:##从主表pcoa_list中求取list，为了二次检验
                    if sp in function_dict.keys():
                        data[sp] = function_dict[sp]
                data_son = SON(data)
                data_list.append(data_son)
            try:
                collection = self.db["pcoa_detail"]
                collection.insert_many(data_list)
                self.bind_object.logger.info("导入pcoa_detail详情表成功！")

            except Exception, e:
                self.bind_object.logger.info("导入pcoa_detail结果表出错:%s" % (e))
            if type in ['specimen']:
                try:
                    main_collection.update_one({'_id': inserted_id}, {'$set': {'pcoa_list': ','.join(sample_list), 'main_id':inserted_id}})
                except Exception, e:
                    self.bind_object.set_error("导入%s结果表出错:%s" % ("pcoa_detail", e))
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
                    "pcoa_id": inserted_id,
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
                collection = self.db["pcoa_detail"]
                collection.insert_many(data_list)
                self.bind_object.logger.info("导入pcoa_detail详情表成功！")
                main_collection = self.db["pcoa"]
                main_collection.update_one({'_id': inserted_id}, {'$set': {'circle': ','.join(sample_list), 'main_id':inserted_id}})

            except Exception, e:
                self.bind_object.set_error("导入pcoa_detail结果表出错:%s" % (e))
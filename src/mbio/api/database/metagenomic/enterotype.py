# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20171114
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import re
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import json
from mbio.packages.metagenomic.id_convert import name2id


class Enterotype(Base):
    def __init__(self, bind_object):
        super(Enterotype, self).__init__(bind_object)
        self._project_type = "metagenomic"

    @report_check
    def add_enterotype(self, main_id=None, main=False, task_id=None, anno_id=None, params=None, anno_type=None,
                       cluster_name=None, spe_name=None, name=None):
        self._tables = []  # 记录存入了哪些表格
        _main_collection = self.db['enterotype']
        if main:
            insert_mongo_json = {
                'project_sn': self.bind_object.sheet.project_sn,
                'task_id': task_id,
                'name': name if name else 'Enterotypes_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end',
                'desc': '',
                'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'anno_type': anno_type
            }
            if anno_id in ['', 'None', None, 'none']:
                pass
            else:
                insert_mongo_json['anno_id'] = ObjectId(anno_id)
            main_id = _main_collection.insert_one(insert_mongo_json).inserted_id
        else:
            if not main_id:
                self.bind_object.set_error('不写入主表时，需要提供主表ID', code="52800801")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        result = _main_collection.find_one({'_id': main_id})
        if not result:
            self.bind_object.set_error('找不到主表id对应的表', code="52800802")
        _main_collection.update_one({'_id': main_id}, {'$set': {'cluster_name': cluster_name}})
        _main_collection.update_one({'_id': main_id}, {'$set': {'spe_name': spe_name}})
        self.bind_object.logger.info("主表导入成功")

    @report_check
    def add_enterotype_detail(self, file_path, table_type, update_id, coll_name='enterotype_detail',
                              main_coll='enterotype', update_column=True, db=None):
        """
        """
        self._tables.append(table_type)
        if not update_id:
            self.bind_object.set_error('需要提供主表ID', code="52800803")
        if not isinstance(update_id, ObjectId):
            update_id = ObjectId(update_id)
        if not db:
            db = self.db
        collection = db[coll_name]
        main_collection = self.db[main_coll]
        main_info = main_collection.find_one({'_id': update_id})
        task_name2id = main_info["task_id"]
        self.sample_2_id = name2id(task_name2id, type="task")
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            columns = all_lines[0].rstrip().split('\t')[1:]
            data_temp = []
            for line in all_lines[1:]:
                values = line.rstrip().split('\t')
                name = values[0]
                if table_type in ['BCA_point','pcoa_point']:
                    name = self.sample_2_id[name]
                insert_data = {
                    'enterotype_id': update_id,
                    'table_type': table_type,
                    'name': name
                }
                if table_type in ['cluster','summary','pcoa_point','BCA_circle','pcoa_circle']:
                    values_dict = dict(zip(columns,values[1:]))
                else:
                    values_dict = dict(zip(columns, map(lambda x: round(float(x), 4), values[1:])))
                data_temp.append(dict(insert_data, **values_dict))
            if data_temp:
                collection.insert_many(data_temp)
            else:
                return None
        if update_column:
            default_column = ['ch', 'BCA_circle', 'BCA_point','pcoa_circle','pcoa_point', 'cluster', 'summary']
            if table_type in default_column:
                main_collection.update_one({'_id': update_id},
                                           {'$set': {table_type: ','.join(columns)}}
                                           )
            else:
                self.bind_object.set_error('错误的表格类型：%s不能在主表中插入相应表头', variables=(table_type),
                                           code="52800804")

    @report_check
    def add_enterotype_detail_cluster(self, id=None, file_path=None, name=None, update_column=True):
        insert_data = list()
        main_collection = self.db['enterotype']
        main_info = main_collection.find_one({'_id': ObjectId(id)})
        task_name2id = main_info["task_id"]
        self.sample_2_id = name2id(task_name2id, type="task")
        with open(file_path, 'rb') as r:
            # head = r.next().strip('\r\n')  # windows换行符
            # head = re.split('\t', head)
            all_lines =  r.readlines()
            head = all_lines[0].rstrip('\r\n')
            head = re.split('\t', head)
            new_head = head[1:-1]  # 第一行后多了的\t，所以为-1
            sample = new_head[:-2]
            for i in range(0, len(sample)):
                sample[i] = self.sample_2_id[sample[i]]
            for line in all_lines[1:101]:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                # classify_list = re.split(r"\s*;\s*", line[0])
                otu_detail = dict()
                otu_detail['tax'] = line[0]
                otu_detail['enterotype_id'] = ObjectId(id)
                otu_detail['cluster_id'] = name
                # for cf in classify_list:
                #     if cf != "":
                #         otu_detail[cf[0:3].lower()] = cf
                for i in range(0, len(sample)):
                    otu_detail[sample[i]] = line[i+1]
                otu_detail[new_head[-2]] = line[-2]
                otu_detail[new_head[-1]] = line[-1]
                insert_data.append(otu_detail)
        try:
            collection = self.db['enterotype_detail_cluster']
            collection.insert_many(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入enterotype_detail_cluster表格信息出错:{}".format(e))
            self.bind_object.set_error("导入enterotype_detail_cluster表信息出错", code="52800805")
        else:
            self.bind_object.logger.info("导入enterotype_detail_cluster表格成功")
        if update_column:
            main_collection.update_one({'_id': ObjectId(id)}, {'$set': {str(name): ','.join(sample)}})

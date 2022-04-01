# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20180310
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId


class AnnoPhi(Base):
    def __init__(self, bind_object):
        super(AnnoPhi, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_anno_phi(self, anno_phi_path, main=False, task_id=None, project_sn=None, main_id=None,
                     params=None, name=None, update_id=None, specimen_id=None):
        self.specimen_id = specimen_id
        if main:
            if task_id is None:
                task_id = self.bind_object.sheet.id
            project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': 'PHI注释',
                'created_ts': created_ts,
                'name': name if name else 'AnnoPhi_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end',
                'version': '3.0',
                'settled_params' : json.dumps({"version": "phi_v4.9"})   # by zhaozhigang 20200929
            }
            collection = self.db['anno_phi']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
            collection.update({'_id': main_id}, {'$set': {'main_id': main_id}})  #guanqing 20180813
        else:
            if main_id is None:
                raise Exception("main为False时需提供main_id!")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        self.add_anno_phi_type(anno_phi_path + '/' + specimen_id + '_phi_stat.xls', main_id)
        self.add_anno_phi_detail(anno_phi_path + '/' + specimen_id + '_phi_anno.xls', main_id, update_id)
        self.add_anno_phi_pie(anno_phi_path + '/' + specimen_id + '_phi_phenotype.xls', main_id)
        self.add_anno_phi_stat(main_id, anno_phi_path + '/sample_stat.xls')
        return main_id

    def add_anno_phi_type(self, type_file, main_id=None):
        data_list = []
        collection_type = self.db['anno_phi_type']
        with open(type_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                data = [("anno_phi_id", main_id)]
                line = line.strip().split('\t')
                data.extend(
                    [("specimen_id", self.specimen_id), ("phi", line[1]), ("des", line[2]), ("anno_num", int(line[3]))])
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection_type.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入anno_phi_type%s信息出错:%s" % (type_file, e))
        else:
            self.bind_object.logger.info("导入anno_phi_type%s信息成功!" % type_file)

    @report_check
    def add_anno_phi_stat(self, main_id, sample_stat):
        data_list = []
        with open(sample_stat, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "anno_phi_id": main_id,
                    "specimen_id": line[0],
                    "gene_num": int(line[1]),
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["anno_phi_stat"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (sample_stat, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % sample_stat)

    def add_anno_phi_pie(self, pie_file, main_id=None):
        data_list = []
        collection_type = self.db['anno_phi_pie']
        with open(pie_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                data = [("anno_phi_id", main_id)]
                line = line.strip().split('\t')
                data.extend(
                    [("specimen_id", self.specimen_id), ("type", line[1]), ("anno_num", int(line[2]))])
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection_type.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入anno_phi_type%s信息出错:%s" % (pie_file, e))
        else:
            self.bind_object.logger.info("导入anno_phi_type%s信息成功!" % pie_file)

    def add_anno_phi_detail(self, detail_file, main_id=None, update_id=None):
        data_list = []
        collection_detail = self.db['anno_phi_detail']
        summary_collection = self.db['anno_summary_detail']
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                data = [("anno_phi_id", main_id)]
                line = line.strip().split('\t')
                data.extend(
                    [("gene_id", line[0]), ("location", line[1]), ("specimen_id", self.specimen_id), ("phi", line[3]),
                     ("protein", line[4]), ("gene", line[5]), ("tax", line[6]), ("pathogen", line[7]),
                     ("type", line[8]), ("host_des", line[9]), ("host_id", line[10]), ("host_species", line[11]),
                     ("function", line[12]),
                     ("identity", float(line[13])), ("evalue", float(line[14])), ("disease", line[15]),
                     ("coverage", float(line[16])), ("score", float(line[17]))])
                data_son = SON(data)
                data_list.append(data_son)
                if update_id:
                    pass
                    # if not isinstance(update_id, ObjectId):
                    #     update_id = ObjectId(update_id)
                    # summary_collection.update_one(
                    #     {'summary_id': update_id, 'gene_id': line[0], 'specimen_id': self.specimen_id},
                    #     {'$set': {'phi_function': line[12]}},
                    #     upsert=False)
                    # self.bind_object.logger.info(update_id)
                    # self.bind_object.logger.info(line[0])
                    # self.bind_object.logger.info(self.specimen_id)
                    # self.bind_object.logger.info(summary_collection.find_one({'summary_id': update_id, 'gene_id': line[0], 'specimen_id': self.specimen_id}))
        try:
            collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入anno_phi_detail%s信息出错:%s" % (detail_file, e))
        else:
            self.bind_object.logger.info("导入anno_phi_detail%s信息成功!" % detail_file)

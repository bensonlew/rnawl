# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20180308
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId


class AnnoVfdb(Base):
    def __init__(self, bind_object):
        super(AnnoVfdb, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_anno_vfdb(self, anno_vfdb_path, main=False, task_id=None, project_sn=None, main_id=None,
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
                'desc': 'VFDB注释',
                'created_ts': created_ts,
                'name': name if name else 'AnnoVfdb_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end',
                "version": "3.0",
                "settled_params": json.dumps({"version": "vfdb_v20200703"})
            }
            collection = self.db['anno_vfdb']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
            collection.update({'_id': main_id},{'$set':{'main_id':main_id}}) #guanqing 20180813
        else:
            if main_id is None:
                raise Exception("main为False时需提供main_id!")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        self.add_anno_vfdb_type(anno_vfdb_path + '/' + specimen_id + '_vfdb_level.xls', main_id)
        self.add_anno_vfdb_detail(anno_vfdb_path + '/' + specimen_id + '_vfdb_anno.xls', main_id, update_id=update_id)
        return main_id

    def add_anno_vfdb_type(self, type_file, main_id=None):
        data_list = []
        collection_type = self.db['anno_vfdb_type']
        with open(type_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                data = [("anno_vfdb_id", main_id)]
                line = line.strip().split('\t')
                data.extend(
                    [("specimen_id", self.specimen_id), ("level1", line[1]), ("level2", line[2]),
                     ("anno_num", int(line[3]))])
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection_type.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入anno_vfdb_type%s信息出错:%s" % (type_file, e))
        else:
            self.bind_object.logger.info("导入anno_vfdb_type%s信息成功!" % type_file)

    def add_anno_vfdb_detail(self, detail_file, main_id=None, update_id=None):
        data_list = []
        collection_detail = self.db['anno_vfdb_detail']
        summary_collection = self.db['anno_summary_detail']
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                data = [("anno_vfdb_id", main_id)]
                line = line.strip().split('\t')
                data.extend(
                    [("gene_id", line[0]), ("location", line[1]), ("specimen_id", self.specimen_id), ("vfdb", line[3]),
                     ("vfs", line[4]), ("species", line[5]), ("des", line[6]), ("level1", line[7]), ("level2", line[8]),
                     ("identity", float(line[9])), ("evalue", float(line[10])), ("score", float(line[11])), ("coverage",float(line[12]))])
                data_son = SON(data)
                data_list.append(data_son)
                if update_id:
                    pass
                    # if not isinstance(update_id, ObjectId):
                    #     update_id = ObjectId(update_id)
                    # summary_collection.update_one(
                    #     {'summary_id': update_id, 'gene_id': line[0], 'specimen_id': self.specimen_id},
                    #     {'$set': {'vfdb_id': line[3]}},
                    #     upsert=False)
        try:
            collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入anno_vfdb_detail%s信息出错:%s" % (detail_file, e))
        else:
            self.bind_object.logger.info("导入anno_vfdb_detail%s信息成功!" % detail_file)
        collection_stat = self.db['anno_vfdb_stat']
        insert_data = {"anno_vfdb_id": main_id, "specimen_id": self.specimen_id, "anno_num": len(lines[1:])}
        try:
            collection_stat.insert_one(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入anno_vfdb_stat%s信息出错:%s" % (detail_file, e))
        else:
            self.bind_object.logger.info("导入anno_vfdb_stat%s信息成功!" % detail_file)

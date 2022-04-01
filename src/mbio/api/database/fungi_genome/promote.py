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


class Promote(Base):
    def __init__(self, bind_object):
        super(Promote, self).__init__(bind_object)
        self._project_type = "fungigenome"

    @report_check
    def add_promote(self, promote_path, main=False, task_id=None, project_sn=None, main_id=None,
                    params=None, name=None, update_id=None, specimen_id=None):
        summary_collection = self.db['anno_summary_detail']
        if main:
            if task_id is None:
                task_id = self.bind_object.sheet.id
            project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '启动子预测',
                'created_ts': created_ts,
                'name': name if name else 'Promote_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end'
            }
            collection = self.db['promote']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
        else:
            if main_id is None:
                #raise Exception("main为False时需提供main_id!")
                self.bind_object.set_error("main为False时需提供main_id!", code="52102401")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        collection_detail = self.db["promote_detail"]
        data_list = []
        with open(promote_path + '/' + specimen_id + '_promoter_result.xls', 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                data = [("promote_id", main_id)]
                line = line.strip().split('\t')
                data.extend(
                    [("gene_id", line[0]), ("location", line[1].capitalize()), ("specimen_id", specimen_id),
                     ("upstream", int(line[3])),
                     ("len", int(line[4])), ("seq", line[5]), ("lsp", int(line[6])), ("lspe", float(line[7])),
                     ("dmaxp", int(line[8])), ("dmx", float(line[9])), ("dave", float(line[10]))])
                data_son = SON(data)
                data_list.append(data_son)
                if update_id:
                    if not isinstance(update_id, ObjectId):
                        update_id = ObjectId(update_id)
                    summary_collection.update_one(
                        {'summary_id': update_id, 'gene_id': line[0], 'specimen_id': specimen_id},
                        {'$set': {'promotor': "Yes"}},
                        upsert=False)
        try:
            collection_detail.insert_many(data_list)
            main_collection = self.db['promote']
            main_collection.update({"_id": main_id}, {"$set": {"status": 'end',"origin_id": main_id}})

        except Exception as e:
            self.bind_object.logger.error("导入promote%s信息出错:%s" % (promote_path, e))
            self.bind_object.set_error("导入promote%s信息出错:%s" , variables=(promote_path, e), code="52102402")
        else:
            self.bind_object.logger.info("导入promote%s信息成功!" % promote_path)
        return main_id

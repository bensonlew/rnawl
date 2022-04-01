# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20180312
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId


class AnnoSignalp(Base):
    def __init__(self, bind_object):
        super(AnnoSignalp, self).__init__(bind_object)
        self._project_type = "fungigenome"

    @report_check
    def add_anno_signalp(self, anno_signalp_path, main=False, task_id=None, project_sn=None, main_id=None,
                         params=None, name=None, specimen_id=None, update_id=None):
        if main:
            if task_id is None:
                task_id = self.bind_object.sheet.id
            project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '分泌蛋白分析',
                'created_ts': created_ts,
                'name': name if name else 'AnnoSignalp_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end'
            }
            collection = self.db['anno_signalp']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
        else:
            if main_id is None:
                #raise Exception("main为False时需提供main_id!")
                self.bind_object.set_error("main为False时需提供main_id!", code="52101101")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        self.add_anno_signalp_detail(anno_signalp_path + '/' + specimen_id + '_SignalP.xls', update_id=update_id, main_id=main_id ,specimen_id=specimen_id )

        return main_id

    @report_check
    def add_anno_signalp_detail(self, detail_file, main_id=None, update_id=None, specimen_id=None):
        data_list = []
        self.specimen_id = specimen_id
        collection_detail = self.db['anno_signalp_detail']
        summary_collection = self.db['anno_summary_detail']
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                data = [("anno_signalp_id", main_id)]
                line = line.strip().split('\t')
                data.extend(
                    [("gene_id", line[0]), ("location", line[1]),("specimen_id",line[2]) ,("cmax", float(line[3])), ("cpos", int(line[4])),
                     ("ymax", float(line[5])), ("ypos", int(line[6])), ("smax", float(line[7])), ("spos", int(line[8])),
                     ("smean", float(line[9])), ("network", line[10])])
                data_son = SON(data)
                data_list.append(data_son)
                if update_id:
                    if not isinstance(update_id, ObjectId):
                        update_id = ObjectId(update_id)
                    summary_collection.update_one(
                        {'summary_id': update_id, 'gene_id': line[0], 'specimen_id': self.specimen_id},
                        {'$set': {'signalp_id': "YES"}},
                        upsert=False)
        try:
            collection_detail.insert_many(data_list)
            main_collection = self.db['anno_signalp']
            main_collection.update({"_id": main_id}, {"$set": {"status": 'end',"origin_id": main_id}})

        except Exception as e:
            self.bind_object.logger.error("导入anno_signalp_detail%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入anno_signalp_detail%s信息出错:%s" , variables=(detail_file, e), code="52101102")
        else:
            self.bind_object.logger.info("导入anno_signalp_detail%s信息成功!" % detail_file)

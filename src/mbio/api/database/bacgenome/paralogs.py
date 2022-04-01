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


class Paralogs(Base):
    def __init__(self, bind_object):
        super(Paralogs, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_paralogs(self, paralogs_path, len, main=False, task_id=None, project_sn=None, main_id=None,
                     params=None, name=None, specimen_id=None):
        self.specimen_id = specimen_id
        #self.search_id = search_id
        self.len = len
        if main:
            if task_id is None:
                task_id = self.bind_object.sheet.id
            project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '旁系同源分析',
                'created_ts': created_ts,
                'name': name if name else 'Paralogs_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end'
            }
            collection = self.db['paralogs']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
            collection.update({'_id': main_id},{'$set':{'main_id':main_id}}) #guanqing 20180813
        else:
            if main_id is None:
                #raise Exception("main为False时需提供main_id!")
                self.bind_object.set_error("main为False时需提供main_id!", code="51403101")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        self.add_paralogs_detail(paralogs_path, main_id, task_id)
        return main_id

    def add_paralogs_detail(self, detail_file, main_id=None, task_id=None):
        data_list = []
        collection_detail = self.db['paralogs_detail']
        collection = self.db['paralogs']
        anno = self.db['anno_summary']
        #summary_id = anno.find_one({"task_id": task_id})['_id']
        summary_id = anno.find_one({"task_id": task_id})['main_id']
        summary_detail = self.db['anno_summary_detail']
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            if len(lines) >1:
                for line in lines[1:]:
                    data = [("paralogs_id", main_id)]
                    line = line.strip().split('\t')
                    search_id = line[0]
                    if line[1] != search_id:
                        if float(line[2]) < 50:
                            continue
                        detail = summary_detail.find_one({"summary_id": summary_id, "specimen_id": self.specimen_id, "gene_id": line[1]})
                        coverage = round(abs(float(line[7])-float(line[6]))/float(line[12]),4)*100
                        #coverage = float(line[14])*100
                        if detail:
                            data.extend(
                                [("gene_id", line[1]), ("location", detail['location']), ("len", float(line[13])),
                                ("des", detail['gene_des']), ("identity", float(line[2])), ("evalue", float(line[10])),
                                ("score", float(line[11])), ("search_id",search_id), ('coverage',float(coverage))])   #zouguanqing 20190401 增加coverage和search_id
                        else:
                            self.bind_object.logger.info({"summary_id": summary_id, "specimen_id": self.specimen_id, "gene_id": line[1]})
                            continue
                        data_son = SON(data)
                        data_list.append(data_son)
            else:
                p_collection = self.db['paralogs']
                p_collection.update({'_id': main_id}, {'$set': {'main_id': main_id}})
                self.bind_object.logger.info("无旁系同源基因")
        if data_list:
            try:
                collection_detail.insert_many(data_list)
                p_collection = self.db['paralogs']
                p_collection.update({'_id': main_id}, {'$set': {'main_id': main_id}})
            except Exception as e:
                self.bind_object.logger.error("导入paralogs_detail%s信息出错:%s" % (detail_file, e))
                self.bind_object.set_error("导入paralogs_detail%s信息出错:%s" , variables=(detail_file, e), code="51403102")
            else:
                self.bind_object.logger.info("导入paralogs_detail%s信息成功!" % detail_file)
        #detail = summary_detail.find_one(
        #    {"summary_id": summary_id, "specimen_id": self.specimen_id, "gene_id": self.search_id})
        # collection.update_one({'_id': main_id}, {
        #     '$set': {'specimen_id': self.specimen_id, 'search_id': self.search_id, 'location': detail['location'],
        #              'len': self.len, 'des': detail['gene_des']}},
        #                       upsert=False)
        collection.update_one({'_id': main_id}, {'$set': {'main_id': main_id}})


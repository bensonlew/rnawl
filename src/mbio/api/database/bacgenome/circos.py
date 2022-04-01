# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last modify: 2018.04.9

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from collections import Counter
from bson.objectid import ObjectId
import shutil


class Circos(Base):
    def __init__(self, bind_object):
        super(Circos, self).__init__(bind_object)
        self._project_type = "bacgenome"

    def add_circos(self, file, params=None, project_sn=None, task_id=None, specimen_id=None, location="Scaffold",
                   main_id=None, main=False, update=True, name=None, link_path=None):
        task_id = task_id if task_id else self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        collection = self.db["circos"]
        if main:
            insert_data = {
                "project_sn": project_sn,
                "task_id": task_id,
                "status": "end",
                "name": name if name else "Circos_Orgin",
                "desc": "circos分析",
                "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
                "created_ts": created_ts,
                "version": "3.0",
            }

            self.bind_object.logger.info(insert_data)
            main_id = collection.insert_one(insert_data).inserted_id
            collection.update({'_id': main_id}, {'$set': {'main_id': main_id}})
        else:
            if not main_id:
                #raise Exception('不写入主表时，需要提供主表ID')
                self.bind_object.set_error('不写入主表时，需要提供主表ID', code="51402201")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
            ##guanqing.zou 20180903
            lab = collection.find_one({"_id":main_id})['lab']
            if 'ncrna' in lab:
                rrna_str = self.rrna_type(task_id,specimen_id)
                rna_str = 'tRNA,'+ rrna_str
                new_lab = lab.replace('ncrna',rna_str)
                collection.update_one({"_id":main_id},{'$set': {'lab':new_lab}})
        collection.update_one({"_id":main_id},{'$set': {'img_path': link_path + '','main_id': main_id}})
        if update:
            self.update_table(link_path + 'pre_file/', task_id, location, specimen_id)
        return main_id

    @report_check
    def update_table(self, link_path, task_id, location, specimen_id):
        circos_table = self.db['circos_table']
        table_id = circos_table.find_one({"task_id": task_id})['_id']
        detail_table = self.db['circos_table_detail']
        if location.startswith('Scaffold'):
            detail_table.update_one({"circos_table_id": table_id, "specimen_id": specimen_id},
                                    {'$set': {'file_path': link_path}}, upsert=True)
        else:
            detail_table.update_one({"circos_table_id": table_id, "location": location, "specimen_id": specimen_id},
                                    {'$set': {'file_path': link_path}}, upsert=True)
    ##guanqing.zou 20180903
    def rrna_type(self,task_id, specimen_id):
        rrna_main = self.db['rrna_predict']
        r_main_id = rrna_main.find_one({"task_id":task_id})['_id']
        rrna_detail = self.db['rrna_predict_detail']
        f_result = rrna_detail.find({"predict_id":r_main_id, 'specimen_id':specimen_id})
        r_type = []
        for i in f_result:
            if i['type'] not in r_type:
                r_type.append(i['type'])
        tmp = len('S_rRNA')
        s_type = sorted(r_type,key=lambda e:int(e[:-tmp]))
        return ','.join(s_type)


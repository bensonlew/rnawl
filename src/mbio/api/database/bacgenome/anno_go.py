# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
from bson.objectid import ObjectId
import pandas as pd


class AnnoGo(Base):
    def __init__(self, bind_object):
        super(AnnoGo, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_anno_go(self, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "GO注释",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "category_list": { 'molecular_function': [], 'cellular_component': [], 'biological_process':[]},
            "name": "AnnoGO_Origin",
            "version" : "3.0"
        }
        collection = self.db["anno_go"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_anno_go_detail(self, inserted_id, specimen_id, anno):
        data_list = []
        go_des_table = self.db["anno_go_specimen"]
        ann = pd.read_table(anno, sep='\t', header=0)
        if len(ann) < 1:
            return
        for i in range(len(ann)):
            data = {
                "go_id":  ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id":  ann["Gene ID"][i],
                "location": ann["Location"][i],
                "gene_des": ann["Gene Description"][i],
                "go_list": ann["GO"][i],
            }
            go_des = []
            for i in  ann["GO"][i].split(';'):
                f_go = go_des_table.find_one({"go_id":ObjectId(inserted_id),"specimen_id":specimen_id,"go":i})
                if f_go:
                    go_des.append(f_go["func_des"])
            data['go_des'] = ';'.join(go_des)
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_go_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)

    @report_check
    def add_anno_go_specimen(self, inserted_id, specimen_id, anno):
        data_list = []
        ann = pd.read_table(anno, sep='\t', header=0)
        if 'GO ID (Lev2)' in ann.columns:
            go_id_head = 'GO ID (Lev2)'
        else:
            go_id_head = 'GO ID'
        if 'GO Term (Lev2)' in ann.columns:
            go_term = 'GO Term (Lev2)'
        else:
            go_term = 'GO Term'

        if len(ann) < 1:
            return
        category = list(set(ann["GO (Lev1)"]))
        for i in range(len(ann)):
            data = {
                "go_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "category": ann["GO (Lev1)"][i],
                "func_des": ann[go_term][i],
                "gene_num": ann["Seq Number"][i],
                "go": ann[go_id_head][i],
                "gene_list": ann["Seq List"][i],
            }
            if 'Ls_Draw' in ann.columns:  #zouguanqing 20190403
                data['draw'] = ann['Ls_Draw'][i]
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_go_specimen"]
            collection.insert_many(data_list)
            main_collection = self.db["anno_go"]
            category_list = main_collection.find_one({"_id": ObjectId(inserted_id)})['category_list']
            #self.bind_object.logger.info(category)
            #self.bind_object.logger.info(">>>>>>>>>>>>>")
            #self.bind_object.logger.info(category_list)
            for cate in category:
                if cate in category_list.keys():
                    category_list[cate].append(specimen_id)
                else:
                    self.bind_object.logger.info("GO主表中不存在该GO (Lev1)，请核实！")
            main_collection.update({"_id": ObjectId(inserted_id)}, {"$set": {"category_list": category_list}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)

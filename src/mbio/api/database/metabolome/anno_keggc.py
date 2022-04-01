# -*- coding: utf-8 -*-
from biocluster.config import Config
from biocluster.api.database.base import Base,report_check
from bson import ObjectId
import json
import types
import datetime
from bson.son import SON


class AnnoKeggc(Base):
    def __init__(self, bind_object):
        super(AnnoKeggc, self).__init__(bind_object)
        self._project_type = 'metabolome'

    def add_anno_keggc_main(self, name, params):
        import datetime
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "desc": "kegg代谢物注释表",
            "name": name,
            "created_ts": created_ts,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "status": "end",
            "version": "2.0"
        }
        collection = self.db['anno_keggc']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({'_id': main_id}, {'$set': {'main_id': main_id}})
        return main_id

    def add_level_detail(self, main_id, level_path):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54701001")
        data_list = list()
        with open(level_path, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")
                data = [
                    ('kegg_id', main_id),
                    ('name', line[3]),
                    ('first_category', line[1]),
                    ('second_category', line[2]),
                    ('compound_id', line[0]),
                    ('metab_id', line[4]),
                    ('hyperlink', line[5])
                ]
                data = SON(data)
                data_list.append(data)

        level_collection = self.db['anno_keggc_level']
        try:
            if data_list:
                level_collection.insert_many(data_list)
            else:
                self.bind_object.logger.info("level表为空")
        except Exception, e:
            self.bind_object.set_error("导入%s信息出错：%s" , variables=(level_path, e), code="54701002")
        else:
            self.bind_object.logger.info("导入level信息成功")

    def add_stat_detail(self, main_id, stat_path):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54701003")
        data_list = list()
        with open(stat_path, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")
                data = [
                    ("kegg_id", main_id),
                    ("first_category", line[0]),
                    ("second_category", line[1]),
                    ("hyperlink", line[2]),
                    ("metab_id", line[3]),
                    ("count", int(line[4]))
                ]
                data = SON(data)
                data_list.append(data)
        stat_coll = self.db['anno_keggc_stat']
        try:
            if data_list:
                stat_coll.insert_many(data_list)
            else:
                self.bind_object.logger.info("stat表为空")
        except Exception, e:
            self.bind_object.set_error("导入%s信息出错： %s" , variables=(stat_path,e), code="54701004")
        else:
            self.bind_object.logger.info("导入stat信息成功")

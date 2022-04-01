# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'
# 20181130

from biocluster.api.database.base import Base, report_check
import datetime,json
from bson.objectid import ObjectId
import pandas as pd

class CommonApi(Base):
    def __init__(self, bind_object):
        super(CommonApi, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_main(self,table, name=None, params=None, others=None,desc=None, task_id=None):
        if task_id == None:
            task_id = self.bind_object.sheet.id
        else:
            task_id = task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": desc if desc else "Job has been finished",
            "name": name,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }

        if others:
            for k in others.keys():      #others {mongo字段名: mongo字段值}
                insert_data[k] = others[k]
        collection = self.db[table]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id}, {'$set': {"main_id": ObjectId(inserted_id)}})
        return inserted_id

    @report_check
    def add_main_detail(self, infile, detail_table, main_id, mongo_key, has_head =True, main_name='main_id',
                        main_table=None, update_dic=None, convert_dic=None):
        """
        :param infile:  需要导表的数据文件
        :param detail_table: detail表名称
        :param main_id: detail表对应的主表_id
        :param mongo_key: detail表中的字段名称,逗号分割
        :param has_head: 数据文件是否有表头，有则删除表头
        :param main_name: detail表中对应主表id的字段名称
        :param main_table: detail表对应的主表名称
        :param update_dic: 需要更新主表中的内容
        :param convert_dic: 更新detail表时数据格式的转换方法，需自定义函数{“字段1”： 转换函数1，“字段2”： 转换函数2}
                            不需要转换则不用写进参数变量中
        :return:
        """
        data_list = []
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        key_list = mongo_key.split(',')
        num = len(key_list)
        with open(infile, 'rb') as r:
            if has_head:
                r.readline()
            for line in r:
                insert_data = {main_name:main_id}
                spline = line.strip("\n\t").split("\t")
                if len(spline) == num:
                    for i in range(num):
                        if key_list[i] != '':
                            insert_data[key_list[i]] = self.convert_format(spline[i], key_list[i], convert_dic)
                else:
                    self.bind_object.set_error("data incomplete: %s ", variables=(line), code="52803302")

                data_list.append(insert_data)
        try:
            collection = self.db[detail_table]
            self.bind_object.logger.info("开始导入%s数据")
            collection.insert_many(data_list)
            if update_dic:
                main_table = self.db[main_table]
                main_table.update({'_id':main_id},{'$set':update_dic})
        except Exception, e:
            self.bind_object.logger.info("导入CommonApi数据出错:%s" % e)
            self.bind_object.set_error("import CommonApi data ERROR", code="52803301")
        else:
            self.bind_object.logger.info("%s导入CommonApi数据成功" % infile)

    def convert_format(self, value, key_name, convert_dic):
        if convert_dic and convert_dic.has_key(key_name):
            try:
                value = convert_dic[key_name](value)
                return value
            except:
                self.bind_object.logger.error("字段%s,数据值%s,用方法%s转换不成功" % (key_name, value, convert_dic[key_name].__name__))
                return value
        else:
            return value

    def update_genome_taxon(self, genome_id, detail_table, task_id=None):
        if not task_id:
            task_id = "_".join(self.bind_object.sheet.id.split("_")[:2])
        result = self.db["genome"].find_one({"task_id": task_id, "genome_id": genome_id})
        data = pd.read_table(detail_table)
        data.sort_values("Identity(%)", ascending=False, inplace=True)
        taxon = data["Taxonomy"][0]
        self.db["genome"].update({"_id": result["_id"]}, {"$set": {"taxon": taxon}})

    def add_sg_status(self, db_name, desc=None, submit_location=None, params=None, table_id=None, table_name=None, genome_id=None, type_name=None):
        task_id ='_'.join(self.bind_object.sheet.id.split("_")[0:2])
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": desc if desc else "Job has been finished",
            "submit_location": submit_location,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "table_id": ObjectId(table_id),
            "table_name": table_name,
            "genome_id": genome_id,
            "type_name": type_name
        }
        collection = self.db[db_name]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def update_status_genome_id(self, table_id, genome_id):
        if not isinstance(table_id, ObjectId):
            table_id = ObjectId(table_id)
        collection = self.db["sg_status"]
        collection.update({"table_id": table_id}, {"$set": {"genome_id": genome_id}})

    def update_mongo(self,db_name,search, change):
        db = self.db[db_name]
        ret = db.find_one(search)
        if ret:
            db.update({"_id":ret["_id"]},{"$set":change})
            return str(ret["_id"])
        else:
            return 'not find'


# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'
# 20181130

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import pandas as pd
import json



class CommonApi(Base):
    def __init__(self, bind_object):
        super(CommonApi, self).__init__(bind_object)
        self._project_type = "metagenomic"

    @report_check
    def add_main(self,table, name=None, params=None, others=None, software_ver={}):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "Job has been finished",
            "name": name,
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }

        if others:
            for k in others.keys():      #others {mongo字段名: mongo字段值}
                insert_data[k] = others[k]
        if software_ver:
            insert_data.update(software_ver)
        collection = self.db[table]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_main_detail(self, infile,detail_table, main_id, mongo_key, has_head =True, main_name='main_id',
                        main_table=None, update_dic=None, convert_dic={}):
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
                spline = line.strip("\n").split("\t")
                if len(spline) == num:
                    for i in range(num):
                        if key_list[i] !='':
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
            self.bind_object.logger.info("导入CommonApi数据成功")


    def convert_format(self, value, key_name, convert_dic):
        if key_name in convert_dic:
            try:
                value = convert_dic[key_name](value)
                return value
            except:
                self.bind_object.logger.error("字段%s,数据值%s,用方法%s转换不成功" % (key_name, value, convert_dic[key_name].__name__))
                return value
        else:
            return value

    def add_detail(self, file_path, detail_col, main_id, main_id_name, main_table=None,
                   has_head=True, mongo_keys=[], key_map={}, update_main={}):
        """
        file_path: 导表文件路径，tab分隔文件
        detail_table: 详情表名称
        main_id: 主表id
        has_head: 是否有表头
        mongo_key: 需要导入的详情表的列名
        key_map: 表头名称修改的字典 {old_name: new_name}
        """
        main_id = ObjectId(main_id)
        header = (None, 0)[int(has_head)]
        detail = self.db[detail_col]
        data = pd.read_csv(file_path, sep='\t', header=header, chunksize=50000)
        for chunk in data:
            chunk = chunk.rename(columns=key_map)
            header = mongo_keys or chunk.columns
            chunk = chunk[header]
            chunk[main_id_name] = main_id
            insert_data = chunk.to_dict('records')
            detail.insert_many(insert_data)
        if main_table and update_main:
            self.db[main_table].update_one({'_id': main_id},
                                           {'$set': update_main})

    def find_one(self, col_name, finder):
        col = self.db[col_name]
        return col.find_one(finder)

    def update_soft_db_info(self, mail_col, main_id, info):
        """
        main_col: 主表 collection 名称
        main_id； 主表id
        info: dict 更新的数据库软件信息 {name: "daimond", params: ["blastp", "evalue < 1e-5"], db: "nr"}
        """
        if isinstance(info, str):
            try:
                info = json.loads(str)
            except Exception as e:
                self.bind_object.set_error("wrong format for soft_db info {}: {}".format(info, e))
        if not isinstance(dict, info):
            self.bind_object.set_error("wrong format for soft_db info {}".format(info))
        info["soft"] = self.get_soft_db(info["soft"])
        info["db"] = self.get_soft_db(info["db"])
        soft_db_info = {"soft_db_info": info}
        self.db[mail_col].update_one({"_id": ObjectId(main_id)},
                                     {"$set", soft_db_info})

    def insert_soft_db_task(self, task_id, soft, database):
        if not isinstance(soft, list) or not isinstance(database, list):
            self.bind_object.set_error("soft and database must be list type")
        soft = [self.get_soft_db(n) for n in soft]
        database = [self.get_soft_db(n) for n in database]
        col = self.db["soft_db_task"]
        if col.find_one({"task_id": task_id}):
            col.delete_one({"task_id": task_id})
        col.insert_one({"task_id": task_id,
                        "soft": soft,
                        "database": database})

    def get_soft_db(self, name, version=None):
        soft_db_file = self._config.SOFTWARE_DIR + "/database/metagenome/soft_db.list"
        df = pd.read_csv(soft_db_file, sep='\t')
        if version:
            df = df[(df["s_name"] == name) & (df["version"] == version)]
        else:
            df = df[df["s_name"] == name]
        df.sort_values("acc_id", axis=0, inplace=True)
        return list(df["acc_id"])[-1]

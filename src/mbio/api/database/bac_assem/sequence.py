# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
# last_modify:20190417
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId
import pandas as pd


class Sequence(Base):
    def __init__(self, bind_object):
        super(Sequence, self).__init__(bind_object)
        self._project_type = "bac_assem"

    @report_check
    def add_sample_info(self, params, task_id=None, project_sn=None, name=None, sample_info={}):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': 'sample_info',
            'created_ts': created_ts,
            'name': name if name else 'sample_info',
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end'
        }
        insert_data.update(sample_info)
        collection = self.db['sample_info']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': main_id},{'$set':{'main_id':main_id}})
        self.db["sg_task"].update({'task_id': task_id}, {'$set': {'comp_pe': 1}})
        return main_id

    @report_check
    def add_raw_stat_detail(self, main_id, stat_file, fastx_dir):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！')
        data = pd.read_table(stat_file, header=0, names=["samp", "lib_type", "lib_name", "insert_size", "read_len", "raw_num", "raw_base", "raw_q20", "raw_q30", "map_sample_detail"])
        data.apply(self.parse_draft, axis=1, main_id=main_id, db_name="raw_stat_detail", fastx_dir=fastx_dir)

    @report_check
    def add_qc_stat_detail(self, main_id, stat_file, fastx_dir, fq_dir):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！')
        data = pd.read_table(stat_file, header=0, names=["samp", "lib_type", "lib_name", "insert_size", "read_len", "clean_num", "clean_base", "clean_q20", "clean_q30", "map_sample_detail"])
        list_data = pd.read_table(os.path.join(fq_dir, "list.txt"))
        data.apply(self.parse_draft, axis=1, main_id=main_id, db_name="qc_stat_detail", fastx_dir=fastx_dir, list_data=list_data)

    @report_check
    def add_comp_stat_detail(self, main_id, stat_dir):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！')
        stat_file = os.path.join(stat_dir, "statistics.xls")
        data = pd.read_table(stat_file, header=0, names=["samp", "lib_type", "raw_num", "raw_base", "largest", "aver"])
        self.bind_object.logger.info(data)
        data.apply(self.parse_comp, axis=1, main_id=main_id, db_name="raw_stat_detail", stat_dir=stat_dir)

    @report_check
    def add_comp_clean_stat_detail(self, main_id, stat_dir, sample, lib_type):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！')
        stat_file2 = os.path.join(stat_dir, sample + ".PacBio_statistics.xls")
        data2 = pd.read_table(stat_file2, header=0, names=["clean_num", "clean_base", "largest", "aver"])
        data2["samp"] = sample
        data2["lib_type"] =lib_type
        data2.apply(self.parse_comp, axis=1, main_id=main_id, db_name="qc_stat_detail", stat_dir=stat_dir)

    @report_check
    def parse_draft(self, df, main_id, db_name, fastx_dir, list_data=None):
        if db_name == "raw_stat_detail":
            insert_data = df[["samp", "lib_type", "lib_name", "insert_size", "read_len", "raw_num", "raw_base", "raw_q20", "raw_q30"]].to_dict()
            insert_data["info_id"] = main_id
            fastx_l_file = os.path.join(fastx_dir, df["map_sample_detail"] + "_l.raw_fastxstat")
            fastx_r_file = os.path.join(fastx_dir, df["map_sample_detail"] + "_r.raw_fastxstat")
        elif db_name == "qc_stat_detail":
            insert_data = df[["samp", "lib_type", "lib_name", "insert_size", "read_len", "clean_num", "clean_base", "clean_q20", "clean_q30"]].to_dict()
            insert_data["info_id"] = main_id
            if df["lib_type"] == "PE":
                file_path = list_data.loc[(list_data["sample"]==df["samp"]) & (list_data["lib"]=="PE") & (list_data["insert"]==df["insert_size"]), "file"].iloc[0]
                file_list = file_path.split(",")
                insert_data["fq1"] = os.path.join(self.bind_object.sheet.output, "data_qc/cleandata", file_list[0])
                insert_data["fq2"] = os.path.join(self.bind_object.sheet.output, "data_qc/cleandata", file_list[1])
            fastx_l_file = os.path.join(fastx_dir, df["map_sample_detail"] + "_l.clean_fastxstat")
            fastx_r_file = os.path.join(fastx_dir, df["map_sample_detail"] + "_r.clean_fastxstat")
        collection = self.db[db_name]
        try:
            seq_id = collection.insert_one(insert_data).inserted_id
            self.add_graphic_detail(seq_id, fastx_l_file, type="left")
            self.add_graphic_detail(seq_id, fastx_r_file, type="right")
        except Exception,e:
            self.bind_object.set_error("%s export error: %s" % (db_name, e))
        else:
            self.bind_object.logger.info("%s export success" % db_name)

    def add_graphic_detail(self, main_id, fastx_file, type):
        data = pd.read_table(fastx_file, header=0, names=["column", "count", "min", "max", "sum", "mean", "q1", "median", "q3", "IQR", "lw", "rw", "a", "c", "g", "t", "n", "max_count"])
        data["data_stat_id"] = main_id
        data["seq_type"] = "raw" if "raw" in fastx_file else "clean"
        data["type"] = type
        data["err_qual"] = 10 ** (data["mean"] / (-10)) * 100
        data_list = data.reindex(columns=["data_stat_id", "seq_type", "type", "column", "min", "max", "err_qual", "q1", "median", "q3", "a", "c", "g", "t", "n"]).to_dict(orient="records")
        collection = self.db["graphic_detail"]
        try:
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("graphic_detail export %s error: %s" % (fastx_file, e))
        else:
            self.bind_object.logger.info("导入graphic_detail信息成功")

    def parse_comp(self, df, main_id, db_name, stat_dir):
        len_file = ''
        type = ''
        if os.path.exists(os.path.join(stat_dir, str(df["samp"]) + ".len.xls")):
            len_file = os.path.join(stat_dir, str(df["samp"]) + ".len.xls")
            type = "raw"
            insert_data = df[["samp", "lib_type", "raw_num", "raw_base", "largest", "aver"]].to_dict()
        elif os.path.exists(os.path.join(stat_dir, str(df["samp"]) + ".clean.len.xls")):
            len_file = os.path.join(stat_dir, str(df["samp"]) + ".clean.len.xls")
            type = "clean"
            insert_data = df[["samp", "lib_type", "clean_num", "clean_base", "largest", "aver"]].to_dict()
        insert_data["info_id"] = main_id
        insert_data["lib_name"] = "-"
        db = self.db[db_name]
        try:
            seq_id = db.insert_one(insert_data).inserted_id
            self.add_clrlen_detail(seq_id, len_file, type)
        except Exception,e:
            self.bind_object.set_error("%s export error: %s" % (db_name, e))
        else:
            self.bind_object.logger.info("%s export success" % db_name)

    def add_clrlen_detail(self, main_id, len_file, type):
        if type in ['raw']:
            data = pd.read_table(len_file, header=0, names=["categories", "data", "line_data"])
            data["data_stat_id"] = main_id
            data["type"] = type
            insert_data = data.to_dict(orient="records")
        elif type in ['clean']:
            data = pd.read_table(len_file, header=0, names=["categories", "data"])
            data["data_stat_id"] = main_id
            data["type"] = type
            insert_data = data.to_dict(orient="records")
        try:
            self.db["clrlen_detail"].insert_many(insert_data)
        except Exception,e:
            self.bind_object.set_error("clrlen_detail export %s error: %s" % (len_file, e))
        else:
            self.bind_object.logger.info("导入clrlen_detail信息成功")


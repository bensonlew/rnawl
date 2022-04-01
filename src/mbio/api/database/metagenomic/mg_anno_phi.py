# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base,report_check
from bson import ObjectId
from bson.son import SON
from mbio.packages.metagenomic.id_convert import name2id
import datetime
import json
import types
import os

class MgAnnoPhi(Base):
    def __init__(self, bind_object):
        super(MgAnnoPhi, self).__init__(bind_object)
        self._project_type = "metagenomic"
        self.sample_2_id = ""  # name2id(self.bind_object.sheet.id, type="task")

    @report_check
    def add_anno_phi(self, params, geneset_id, anno_file, name="PHI_Origin"):
        """
        params根据接口设计调整
        anno_file路径同其他主表
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.sample_2_id = name2id(task_id, type="task")
        specimen = ",".join(self.sample_2_id.values())
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "desc": "phi注释",
            "created_ts": created_ts,
            "status": "end",
            "name": name,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "geneset_id": geneset_id,
            "specimen": specimen,
            "anno_file": anno_file,
            "lowest_level": "protein",
            "is_origin": 1,
            "settled_params": json.dumps({"version": "phi_v4.9"})  # by zhaozhigang 20200929
        }
        collection = self.db["anno_phi"]
        main_id = collection.insert_one(insert_data).inserted_id
        return main_id

    @report_check
    def add_anno_phi_host(self, main_id, table):
        """
        宿主个数统计表
        :param main_id:
        :param table:
        :return:
        """
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串！", code="52803801")
        data_list = list()
        result = self.db["anno_phi"].find_one({'_id': main_id})
        self.sample_2_id = name2id(result["task_id"], type="task")
        with open(table, "r") as f1:
            file = f1.readlines()
            head_list = file[0].split("\t")[1:-2]
            for line in file[1:]:
                line = line.strip().split("\t")
                data = [
                    ("host", line[0]),
                    ("host_desc", line[-1]),
                    ("Total", float(line[-2])),
                    ("phi_id", main_id)
                ]
                for index,i in enumerate(line[1:-2]):
                    data.append((self.sample_2_id[head_list[index]], float(i)))
                data = SON(data)
                data_list.append(data)
        detail_coll = self.db["anno_phi_host"]
        try:
            detail_coll.insert_many(data_list)
        except Exception,e:
            self.bind_object.set_error("anno_phi_host error: %s" , variables=( e), code="52803802")
        else:
            self.bind_object.logger.info("导入anno_phi_host成功")

    @report_check
    def add_anno_phi_pathogen(self, main_id, table):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串！", code="52803803")
        data_list = list()
        result = self.db["anno_phi"].find_one({'_id': main_id})
        self.sample_2_id = name2id(result["task_id"], type="task")
        with open(table, "r") as f1:
            file = f1.readlines()
            head_list = file[0].split("\t")[1:-1]
            for line in file[1:]:
                line = line.strip().split("\t")
                data = [
                    ("pathogen", line[0]),
                    ("Total", float(line[-1])),
                    ("phi_id", main_id)
                ]
                for index,i in enumerate(line[1:-1]):
                    data.append((self.sample_2_id[head_list[index]], float(i)))
                data = SON(data)
                data_list.append(data)
        detail_coll = self.db["anno_phi_pathogen"]
        try:
            detail_coll.insert_many(data_list)
        except Exception,e:
            self.bind_object.set_error("anno_phi_pathogen error: %s" , variables=( e), code="52803804")
        else:
            self.bind_object.logger.info("导入anno_phi_pathogen成功")

    @report_check
    def add_anno_phi_phenotype(self, main_id, table):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串！", code="52803805")
        data_list = list()
        result = self.db["anno_phi"].find_one({'_id': main_id})
        self.sample_2_id = name2id(result["task_id"], type="task")
        with open(table, "r") as f1:
            file = f1.readlines()
            head_list = file[0].split("\t")[1:-1]
            for line in file[1:]:
                line = line.strip().split("\t")
                data = [
                    ("phenotype", line[0]),
                    ("Total", float(line[-1])),
                    ("phi_id", main_id)
                ]
                for index,i in enumerate(line[1:-1]):
                    data.append((self.sample_2_id[head_list[index]], float(i)))
                data = SON(data)
                data_list.append(data)
        detail_coll = self.db["anno_phi_phenotype"]
        try:
            detail_coll.insert_many(data_list)
        except Exception,e:
            self.bind_object.set_error("anno_phi_phenotype error: %s" , variables=( e), code="52803806")
        else:
            self.bind_object.logger.info("导入anno_phi_phenotype成功")

    @report_check
    def add_anno_phi_detail(self, main_id, table):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串！", code="52803807")
        data_list = list()
        # result = self.db["anno_phi"].find_one({'_id': main_id})
        with open(table, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")
                data = [
                    ("phi_id", main_id),
                    ("gene_id", line[0]),
                    ("phi_accession", line[1]),
                    ("protein_id", line[2]),
                    ("gene_name", line[3]),
                    ("tax", line[4]),
                    ("pathogen", line[5]),
                    ("phenotype", line[6]),
                    ("host_desc", line[7]),
                    ("host_id", line[8]),
                    ("host", line[9]),
                    ("function", line[10]),
                ]
                if len(line) == 14:
                    data += [
                        ("disease", line[11]),
                        ("identity", float(line[12])),
                        ("align_length", float(line[13])),
                    ]
                else:
                    data += [
                        ("identity", float(line[11])),
                        ("align_length", float(line[12])),
                    ]
                data = SON(data)
                data_list.append(data)
        detail_coll = self.db["anno_phi_detail"]
        try:
            detail_coll.insert_many(data_list)
            self.db['anno_phi'].update_one({'_id': main_id}, {'$set':{'main_id':main_id,"settled_params": json.dumps({"version": "phi_v4.9"})}})
        except Exception,e:
            self.bind_object.set_error("anno_phi_detail error: %s" , variables=( e), code="52803808")
        else:
            self.bind_object.logger.info("导入anno_phi_detail成功")

    @report_check
    def update_anno_file(self, main_id, table):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串！", code="52803809")
        collection = self.db["anno_phi"]
        collection.update({'_id': main_id}, {'$set':{"anno_file": table}})

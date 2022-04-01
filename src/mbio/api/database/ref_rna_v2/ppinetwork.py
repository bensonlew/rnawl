# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last_modify:20170425
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import datetime
import json

class Ppinetwork(Base):
    def __init__(self, bind_object):
        super(Ppinetwork, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + '_ref_rna'

    @report_check
    def add_ppi_main_id(self, geneset_id, combine_score, gene_type, species):
        params_dict = dict()
        params_dict["combine_score"] = combine_score,
        params_dict["gene_type"] = gene_type,
        params_dict["geneset_id"] = geneset_id,
        params_dict["species"] = species,
        params_dict["submit_location"] = "ppinetwork",
        params_dict["task_type"] = "workflow"
        data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "status": "end",
            "name": 'PPINETWORK_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3],
            "geneset_id": geneset_id,
            "desc": "",
            "created_ts": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            "params": json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
        }
        col = self.db["sg_ppinetwork"]
        main_id = col.insert_one(data).inserted_id
        self.bind_object.logger.info("主表创建成功")
        return main_id

    # @report_check
    def add_network_attributes(self, file1_path, file2_path, table_id=None, major=False):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                raise Exception("table_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        data_list1 = []
        with open(file1_path, "rb") as r, open(file2_path, "rb") as w:
            data1 = r.readlines()
            data2 = w.readlines()
            for line1 in data1:
                temp1 = line1.rstrip().split("\t")
            for line2 in data2:
                temp2 = line2.rstrip().split("\t")
                data_list1.append(temp2[1])
            data = [("ppi_id", table_id), ("num_of_nodes", eval(data_list1[0])), ("num_of_edges", eval(data_list1[1])),
                    ("average_node_degree", eval(data_list1[2])), ("average_path_length", eval(data_list1[3])),
                    ("average_cluster_coefficient", "None" if data_list1[4] == "NA" else eval(data_list1[4])),
                    ("transitivity", eval(temp1[1]))]
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["sg_ppinetwork_structure_attributes"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入%s信息出错:%s" % (file1_path, e))
            self.bind_object.set_error("导入%s信息出错:%s" % (file2_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file1_path)
            self.bind_object.logger.info("导入%s信息成功!" % file2_path)
        return data_list, table_id

    # @report_check
    def add_network_cluster_degree(self, file1_path, file2_path, table_id=None, major=False):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                raise Exception("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        with open(file1_path, 'rb') as r, open(file2_path, 'rb') as w:
            data1 = r.readlines()[1:]
            data2 = w.readlines()[1:]
            for line2 in data2:
                temp2 = line2.rstrip().split("\t")
                for line1 in data1:
                    temp1 = line1.rstrip().split("\t")
                    if temp1[1] == temp2[1]:
                        data = [("ppi_id", table_id), ("node_id", eval(temp1[0])), ("node_name", temp1[1]),
                                ("degree", eval(temp1[2])), ("clustering", eval(temp2[2]))]
                        data_son = SON(data)
                        data_list.append(data_son)
        try:
            collection = self.db["sg_ppinetwork_structure_node"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入%s信息出错:%s" % (file1_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file1_path)
        return data_list

    # @report_check
    def add_network_centrality(self, file_path, table_id=None, major=False):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                raise Exception("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        with open(file_path, 'rb') as r:
            i = 0
            for line in r:
                if i == 0:
                    i = 1
                else:
                    line = line.strip('\n')
                    line_data = line.split('\t')
                    data = [("ppi_id", table_id), ("node_id", eval(line_data[0])),
                            ("node_name", line_data[1]), ("degree_centrality", eval(line_data[2])),
                            ("closeness_centrality", eval(line_data[3])), ("betweenness_centrality", eval(line_data[4]))]
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["sg_ppinetwork_centrality_node"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    # @report_check
    def add_node_table(self, file_path, table_id=None, major=False):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                raise Exception("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        with open(file_path, 'rb') as r:
            data_line = r.readlines()[1:]
            for line in data_line:
                line_data = line.strip().split('\t')
                data = [("ppi_id", table_id), ("node_name", line_data[0]), ("degree", eval(line_data[1])),
                        ("gene_id", line_data[2]), ("string_id", line_data[3])]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["sg_ppinetwork_node_table"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    #@report_check
    def add_edge_table(self, file_path, table_id=None, major=False):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                raise Exception("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        with open(file_path, 'rb') as r:
            data1 = r.readlines()[1:]
            for line in data1:
                line_data = line.strip().split('\t')
                data = [("ppi_id", table_id), ("from", line_data[0]), ("to", line_data[1]),
                        ("combined_score", eval(line_data[2]))]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["sg_ppinetwork_structure_link"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    @report_check
    def add_degree_distribution(self, file_path, table_id=None, major=False):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                raise Exception("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        with open(file_path, 'rb') as r:
            data_line = r.readlines()[1:]
            for line in data_line:
                line_data = line.strip().split('\t')
                data = [("ppi_id", table_id), ("degree", eval(line_data[0])), ("node_num", eval(line_data[1]))]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["sg_ppinetwork_distribution_node"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

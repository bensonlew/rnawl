# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
from biocluster.api.database.base import Base, report_check
import re
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
import gridfs
import datetime
import os
# from biocluster.config import Config


class CorrNetwork(Base):
    def __init__(self, bind_object):
        super(CorrNetwork, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB

    # @report_check
    def add_network_attributes(self, file_path, table_id = None, group_id = None, from_otu_table = None, level_id = None, major = False, params=None):
        if major:
            table_id = self.create_corrnetwork(self, params, group_id, from_otu_table, level_id)
        else:
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
                else:
                    self.bind_object.set_error("table_id必须为ObjectId对象或者其对应的字符串！", code="51000701")
        data_list = []
        with open(file_path, "rb") as r:
            data = r.readlines()[1:]
            for line in data:
                line = line.strip().split('\t')
            data = [("corr_network_id", table_id), ("transitivity", eval(line[0])),
                    ("diameter", str(line[1])), ("average_shortest_path_length", str(line[2]))]
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["sg_corr_network_structure_attributes"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_corr_network"]
            #main_collection.update({"_id": table_id}, {"$set": {"main_id": table_id}})
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list, table_id

    # @report_check
    def add_network_degree_distribution(self, file_path, table_id = None, group_id = None, from_otu_table = None, level_id = None, major = False,params=None):
        if major:
            table_id = self.create_corrnetwork(self, params, group_id, from_otu_table, level_id)
        else:
            if table_id is None:
                self.bind_object.set_error("major为False时需提供table_id!", code="51000702")
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
                else:
                    self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51000703")
        data_list = []
        with open(file_path, 'rb') as r:
            data = r.readlines()[1:]
            for line in data:
                line = line.strip().split("\t")
                data = [("corr_network_id", table_id), ("degree", line[0]), ("num", line[1])]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["sg_corr_network_distribution_node"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    # @report_check
    def add_network_centrality(self, file_path, table_id = None, group_id = None, from_otu_table = None, level_id = None, major = False,params=None):
        if major:
            table_id = self.create_corrnetwork(self, params, group_id, from_otu_table, level_id)
        else:
            if table_id is None:
                self.bind_object.set_error("major为False时需提供table_id!", code="51000704")
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
                else:
                    self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51000705")
        data_list = []
        with open(file_path, 'rb') as r:
            i = 0
            for line in r:
                if i == 0:
                    i = 1
                else:
                    line = line.strip('\n')
                    line_data = line.split('\t')
                    data = [("corr_network_id", table_id), ("node_id", eval(line_data[0])),
                            ("node_name", line_data[1]),("degree_centrality", eval(line_data[2])),
                            ("closeness_centrality", eval(line_data[3])), ("betweenness_centrality", eval(line_data[4]))]
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["sg_corr_network_centrality_node"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    # @report_check
    def add_network_cluster_degree(self,  file1_path, file2_path, params=None, group_id=None, from_otu_table=None, table_id=None, major=False,level_id=None):
        #file1_path:节点的degree表，file2_path:节点的cluster表
        if major:
            table_id = self.create_corrnetwork(self, params, group_id, from_otu_table, level_id)
        else:
            if table_id is None:
                self.bind_object.set_error("major为False时需提供table_id!", code="51000706")
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
                else:
                    self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51000707")
        data_list = []
        with open(file1_path, 'rb') as r, open(file2_path, 'rb') as w:
            data1 = r.readlines()[1:]
            data2 = w.readlines()[1:]
            for line2 in data2:
                temp2 = line2.strip().split("\t")
                for line1 in data1:
                    temp1 = line1.strip().split("\t")
                    if temp1[1] == temp2[1]:
                        data = [("corr_network_id", table_id), ("node_id", eval(temp1[0])), ("node_name", temp1[1]), ("degree", eval(temp1[2])),
                                ("clustering", eval(temp2[2]))]
                        data_son = SON(data)
                        data_list.append(data_son)
        try:
            collection = self.db["sg_corr_network_structure_node"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file1_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file1_path)
        return data_list

    #@report_check
    def add_network_links_table(self, file_path, node_id_file, table_id=None, group_id=None, from_otu_table=None, level_id=None,
                               major=False,params=None, ):
        if major:
            table_id = self.create_corrnetwork(self, params, group_id, from_otu_table, level_id)
        else:
            if table_id is None:
                self.bind_object.set_error("major为False时需提供table_id!", code="51000708")
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
                else:
                    self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51000709")
        id_name = {}
        with open(node_id_file) as fr:
            fr.readline()
            for line in fr:
                line = line.strip()
                spline = line.split("\t")
                if spline[1] not in id_name:
                    id_name[spline[1]] = spline[0]

        data_list = []
        with open(file_path, 'rb') as r:
            data1 = r.readlines()[3:]
            for line in data1:
                line_data = line.strip().split('\t')
                data = [("corr_network_id", table_id), ("source", int(id_name[line_data[0]])), ("target", int(id_name[line_data[1]])),
                        ("value", eval(line_data[2]))]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["sg_corr_network_structure_link"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    # @report_check
    def add_network_abundance_table(self, file_path, table_id=None, group_id=None, from_otu_table=None, level_id=None,
                               major=False,params=None):
        if major:
            table_id = self.create_corrnetwork(self, params, group_id, from_otu_table, level_id)
        else:
            if table_id is None:
                self.bind_object.set_error("major为False时需提供table_id!", code="51000710")
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
                else:
                    self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51000711")
        data_list = []
        with open(file_path, 'rb') as r:
            data = r.readlines()[1:]
            for line in data:
                line = line.strip().split("\t")
                data = [("corr_network_id", table_id), ("node_name", line[0]), ("abundance", eval(line[1])), ("phylum", line[2])]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["sg_corr_network_structure_abundance"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.error("导入%s信息成功!" % file_path)
        return data_list

    #@report_check
    def create_corrnetwork(self, params, group_id=0, from_otu_table=0, name=None, level_id=0):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51000712")
        if group_id not in ["all","All","ALL"]:
            if group_id != 0 and not isinstance(group_id, ObjectId):
                if isinstance(group_id, StringTypes):
                    group_id = ObjectId(group_id)
                else:
                    self.bind_object.set_error("group_detail必须为ObjectId对象或其对应的字符串!", code="51000713")
        if level_id not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level_id)
            self.bind_object.set_error("level参数不在允许的范围内", code="51000714")

        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "corrnetwork分析"
        if name:
            table_name = name
        else:
            table_name = self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else "corrnetwork_" + \
            datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "otu_id": from_otu_table,
            "group_id": group_id,
            "name": table_name,
            "params": params,
            "level_id": level_id,
            "desc": desc,
            "status": "end",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_corr_network"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        #collection.update({"_id": inserted_id}, {"$set": {"main_id": inserted_id}})
        return inserted_id


    @report_check
    def check_id(self, object_id):
        if not isinstance(object_id, ObjectId):
            object_id = ObjectId(object_id)
        return object_id

    ##heatmap
    @report_check
    def add_heatmap_corr_detail(self, main_id, corr_file, p_file, tree_file=None):
        main_id = self.check_id(main_id)
        if not os.path.exists(corr_file):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(corr_file), code="51000715")
        data_list = []
        result = self.db['sg_corr_network'].find_one({'_id': main_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(main_id), code="51000716")
        else:
            task_id = result['task_id']

        with open(corr_file, 'rb') as f:
            lines = f.readlines()

        head = lines[0]
        sams = head.rstrip().split("\t")
        sample_map = {}
        spe_list = []
        num = 1
        for  s in sams:
            if s == '':
                continue
            tmp = '{}:vs{}'.format(s,str(num))
            spe_list.append(tmp)
            sample_map[s] = 'vs%s'% num
            num += 1
        new_spe_list = '|'.join(spe_list)

        insert_list = []
        for line in lines[1:]:
            line = line.strip().split('\t')
            name = line[0]
            insert_data = {
                'corr_network_id': main_id,
                'spe': name,
                'type': 'corr'
            }
            for i in range(1, len(line)):
                sam_corr = float(line[i])
                k = sample_map[sams[i]]
                insert_data[k] = sam_corr
            insert_list.append(insert_data)

        try:
            collection = self.db['sg_corr_network_heatmap']
            collection.insert_many(insert_list)
        except Exception as e:
            self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(corr_file, e), code="51000717")
        else:
            pass
        self.bind_object.logger.info("导入表格%s信息成功!" % corr_file)

        data_list = []
        with open(p_file, 'rb') as f:
            head = f.next()
            sams = head.rstrip().split("\t")
            for line in f:
                line = line.strip().split('\t')
                name = line[0]
                insert_data = {
                    'corr_network_id': main_id,
                    'spe': name,
                    'type': 'pvalue'
                }
                for i in range(1, len(line)):
                    sam_corr = float(line[i])
                    k = sample_map[sams[i]]
                    insert_data[k] = sam_corr
                data_list.append(insert_data)
        try:
            collection = self.db['sg_corr_network_heatmap']
            collection.insert_many(data_list)

        except Exception as e:
            self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(corr_file, e), code="51000718")
        else:
            pass
        self.bind_object.logger.info("导入表格%s信息成功!" % corr_file)


        self.bind_object.logger.info("更新主表sg_corr_network")
        main_collection = self.db['sg_corr_network']
        if tree_file:
            tree = open(tree_file).readline()
            main_collection.update({'_id':main_id},{"$set":{"spe_list": new_spe_list, "spe_tree": tree}})
        else:
            main_collection.update({'_id':main_id},{"$set":{"spe_list": new_spe_list}})

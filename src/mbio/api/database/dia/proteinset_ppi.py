# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last_modify: liubinxu 201804027
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import pandas as pd
import datetime
import json
import unittest
import os
from mbio.api.database.dia.proteinset import Proteinset

class ProteinsetPpi(Proteinset):
    def __init__(self, bind_object):
        super(ProteinsetPpi, self).__init__(bind_object)
        self.node2id = dict()
        self.node2acc = dict()
        self.acc2id = dict()

    @report_check
    def add_ppi_main_id(self, geneset_id, combine_score, gene_type, species):
        params_dict = dict()
        params_dict["combine_score"] = combine_score
        params_dict["gene_type"] = gene_type
        params_dict["geneset_id"] = geneset_id
        params_dict["species"] = species
        params_dict["submit_location"] = "ppinetwork"
        params_dict["task_type"] = "workflow"
        params_dict["task_id"] = self.bind_object.sheet.id
        data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "status": "end",
            "name": 'ProteinsetPpi_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3],
            "geneset_id": geneset_id,
            "desc": "",
            "created_ts": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            "params": json.dumps(params_dict, sort_keys=True, separators=(',', ':')),
            "gene_type": gene_type
        }
        col = self.db["sg_ppinetwork"]
        main_id = col.insert_one(data).inserted_id
        self.bind_object.logger.info("主表创建成功")
        return main_id

    @report_check
    def add_network_attributes(self, file2_path, table_id=None, major=False):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                raise Exception("table_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        data_list1 = []
        with open(file2_path, "rb") as w:
            data2 = w.readlines()
            for line2 in data2:
                temp2 = line2.rstrip().split("\t")
                data_list1.append(temp2[1])
            data = [("ppi_id", table_id), ("num_of_nodes", eval(data_list1[0])), ("num_of_edges", eval(data_list1[1])),
                    ("average_node_degree", eval(data_list1[2])), ("average_path_length", eval(data_list1[3])),
                    ("average_cluster_coefficient", "None" if data_list1[4] == "NA" else eval(data_list1[4])),
                    ("transitivity", "None" if data_list1[5] == "NA" else eval(data_list1[5]))]
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["sg_proteinset_ppi_stat"]
            collection.insert_many(data_list)
        except Exception as e:
            # self.bind_object.set_error("导入%s信息出错:%s" % (file1_path, e))
            self.bind_object.set_error("导入%s信息出错:%s" % (file2_path, e))
        else:
            # self.bind_object.logger.info("导入%s信息成功!" % file1_path)
            self.bind_object.logger.info("导入%s信息成功!" % file2_path)
        return data_list, table_id

    @report_check
    def add_network_cluster_degree(self, file1_path, file2_path, file3_path, table_id=None, major=False):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                raise Exception("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        with open(file1_path, 'rb') as r, open(file2_path, 'rb') as w, open(file3_path, "rb") as m:
            data1 = r.readlines()[1:]
            data2 = w.readlines()[1:]
            data3 = m.readlines()[1:]
            for line2 in data2:
                temp2 = line2.rstrip().split("\t")
                for line1 in data1:
                    temp1 = line1.rstrip().split("\t")
                    for line3 in data3:
                        temp3 = line3.rstrip().split("\t")
                        if temp1[1] == temp2[1] and temp1[1] == temp3[0]:
                            node_id = self.acc2id[temp3[2]]
                            data = [("ppi_id", table_id), ("node_id", node_id), ("node_name", temp1[1]),
                                    ("degree", eval(temp1[2])), ("clustering", eval(temp2[2])), ("accession_id", temp3[2])]
                            data_son = SON(data)
                            data_list.append(data_son)
        try:
            collection = self.db["sg_proteinset_ppi_nodes_cluster"]
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入%s信息出错:%s" % (file1_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file1_path)
        return data_list

    @report_check
    def add_network_centrality(self, file_path, file1_path, table_id=None, major=False):
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
                    with open(file1_path, 'rb') as m:
                        data1 = m.readlines()[1:]
                        for line1 in data1:
                            temp = line1.rstrip().split("\t")
                            if temp[1] == line_data[1]:
                                node_id = line_data[0]
                                data = [("ppi_id", table_id), ("node_id", node_id),
                                        ("node_name", line_data[1]), ("degree_centrality", eval(line_data[2])),
                                        ("closeness_centrality", eval(line_data[3])),
                                        ("betweenness_centrality", eval(line_data[4])), ("accession_id", temp[3])]
                                data_son = SON(data)
                                data_list.append(data_son)
        try:
            collection = self.db["sg_proteinset_ppi_nodes_centrality"]
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    @report_check
    def add_node_table(self, file_path, group_id=None, table_id=None, major=False):
        self.bind_object.logger.info(file_path)
        update_status = False
        if not table_id:
            update_status = True
            pset_name = os.path.dirname(file_path).split('/')[-2]
            proteinset_id = self.get_proteinset_id(pset_name).__str__()
            with open(file_path) as nr:
                _ = nr.readline()
                line = nr.readline()
                if line:
                    specie = nr.readline().split('.')[0]
                else:
                    specie = "9096"
            params = {
                "proteinset_id": proteinset_id,
                "combine_score":"300",
                "species": specie,
                "submit_location": "ppinetwork",
                "task_id": self.task_id,
                "task_type": "2",
                "taxon":"All",
                # "version":'v3'
            }
            name = 'DiffPpi_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3] + '_' + pset_name
            table_id = self.add_main_table('sg_proteinset_ppi', params, name)
        else:
            table_id = ObjectId(table_id)
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                raise Exception("table_id必须为ObjectId对象或其对应的字符串!")
        nodes_table = pd.read_table(file_path, header=0)
        nodes_table.sort_values(['degree'],  ascending = [False], inplace=True)
        nodes_table['node_id'] = range(0, len(nodes_table))
        nodes_table.rename(columns={'node':'node_name','STRING_id':'string_id'}, inplace=True)
        nodes_table['ppi_id'] = table_id
        self.node2acc =dict(zip(nodes_table['string_id'],nodes_table['accession_id']))
        self.node2id =dict(zip(nodes_table['string_id'],nodes_table['node_id']))
        self.acc2id =dict(zip(nodes_table['accession_id'],nodes_table['node_id']))
        nodes_table_list = nodes_table.to_dict('records')
        self.create_db_table('sg_proteinset_ppi_node', nodes_table_list)
        if update_status:
            if group_id:
                params.update({'group_id': str(group_id)})
            self.add_sg_status(submit_location='ppinetwork', params=params, table_id=ObjectId(table_id),
                               table_name=name, type_name='sg_proteinset_ppi')
        self.update_db_record('sg_proteinset_ppi', ObjectId(table_id), main_id=ObjectId(table_id))
        return table_id

    @report_check
    def add_edge_table(self, file_path, table_id=None, major=False):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                raise Exception("table_id必须为ObjectId对象或其对应的字符串!")
        edge_table = pd.read_table(file_path, header=0)
        edge_table['from_accession'] = edge_table['from'].map(lambda x:self.node2acc[x])
        edge_table['to_accession'] = edge_table['to'].map(lambda x:self.node2acc[x])
        edge_table['from_id'] = edge_table['from'].map(lambda x:self.node2id[x])
        edge_table['to_id'] = edge_table['to'].map(lambda x:self.node2id[x])
        transfer_list = []
        for method in ['database', 'coexpression', 'experiments', 'neighborhood', 'textmining']:
            edge_table[method] = edge_table[method + '_transferred'].map(lambda x:float(x))
            transfer_list.append(method + '_transferred')
        for method in ['cooccurence', 'homology', 'combined_score', 'fusion']:
            edge_table[method] = edge_table[method].map(lambda x:float(x))
        edge_table.drop(transfer_list, axis=1, inplace=True)

        edge_table.rename(columns={'preferred_name.x':'from_name','preferred_name.y':'to_name'}, inplace=True)
        edge_table['ppi_id'] = table_id
        edge_table_list = edge_table.to_dict('records')
        self.create_db_table('sg_proteinset_ppi_edge', edge_table_list)

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
            collection = self.db["sg_proteinset_ppi_nodes_degreestat"]
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

class TestFunction(unittest.TestCase):
    """
    测试导表函数

    """
    def test_mongo(test):
        from mbio.workflows.itraq_and_tmt.itraq_test_api import ItraqApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            "id": "itraq_tmt",
            #+ str(random.randint(1,10000)),
            #"id": "denovo_rna_v2",
            "project_sn": "itraq_tmt",
            #+ str(random.randint(1,10000)),
            "type": "workflow",
            "name": "itraq_and_tmt.protein_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = ItraqApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("itraq_and_tmt.proteinset_ppi")

        test_dir = '/mnt/ilustre/users/sanger-dev/workspace/20180409/ProteinsetPpi_itraq_tmt_4321_9222/PpinetworkAnalysis/output'
        all_nodes_path = test_dir + '/ppinetwork_predict/all_nodes.txt'   # 画图节点属性文件
        interaction_path = test_dir + '/ppinetwork_predict/interaction_detail.txt'  # 画图的边文件
        network_stats_path = test_dir + '/ppinetwork_predict/network_stats.txt'  # 网络全局属性统计
        network_centrality_path = test_dir + '/ppinetwork_topology/protein_interaction_network_centrality.txt'
        network_clustering_path = test_dir + '/ppinetwork_topology/protein_interaction_network_clustering.txt'
        network_transitivity_path = test_dir + '/ppinetwork_topology/protein_interaction_network_transitivity.txt'
        degree_distribution_path = test_dir + '/ppinetwork_topology/protein_interaction_network_degree_distribution.txt'
        network_node_degree_path = test_dir + '/ppinetwork_topology/protein_interaction_network_node_degree.txt'
        if not os.path.isfile(all_nodes_path):
            raise Exception("找不到报告文件:{}".format(all_nodes_path))
        if not os.path.isfile(interaction_path):
            raise Exception("找不到报告文件:{}".format(interaction_path))
        if not os.path.isfile(network_stats_path):
            raise Exception("找不到报告文件:{}".format(network_stats_path))
        if not os.path.isfile(network_clustering_path):
            raise Exception("找不到报告文件:{}".format(network_clustering_path))
        if not os.path.isfile(network_centrality_path):
            raise Exception("找不到报告文件:{}".format(network_centrality_path))
        if not os.path.isfile(network_transitivity_path):
            raise Exception("找不到报告文件:{}".format(network_transitivity_path))
        if not os.path.isfile(degree_distribution_path):
            raise Exception("找不到报告文件:{}".format(degree_distribution_path))
        if not os.path.isfile(network_node_degree_path):
            raise Exception("找不到报告文件:{}".format(network_node_degree_path))

        print('start insert')
        ppi_id = "5acac971a4e1af6b2a6b054c"
        wf.test_api.add_node_table(file_path=all_nodes_path, table_id=ppi_id)   # 节点的属性文件（画网络图用）
        wf.test_api.add_edge_table(file_path=interaction_path, table_id=ppi_id)  # 边信息
        wf.test_api.add_network_attributes(file1_path=network_transitivity_path, file2_path=network_stats_path,
                                              table_id=ppi_id)  # 网络全局属性
        wf.test_api.add_network_cluster_degree(file1_path=network_node_degree_path,
                                                  file2_path=network_clustering_path, file3_path=all_nodes_path,
                                                  table_id=ppi_id)  # 节点的聚类与degree，画折线图
        wf.test_api.add_network_centrality(file_path=network_centrality_path, file1_path=all_nodes_path,
                                              table_id=ppi_id)  # 中心信息
        wf.test_api.add_degree_distribution(file_path=degree_distribution_path, table_id=ppi_id)
        print('end insert')

if __name__ == '__main__':
    unittest.main()

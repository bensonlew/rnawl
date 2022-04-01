# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metagenomic.id_convert import name2id
from mbio.packages.metagenomic.id_convert import id2name
import json
import pandas as pd


class NetworkCorfd(Base):
    def __init__(self, bind_object):
        super(NetworkCorfd, self).__init__(bind_object)
        self._project_type = "metaasv"

    @report_check
    def add_network_corfd(self, attributes_file, geneset_id=None, main=True, main_table_id=None, params=None,
                          name = None):
        if not os.path.exists(attributes_file):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('attributes_file所指定的路径不存在，请检查！')
        if main:
            if not isinstance(geneset_id, ObjectId):
                if isinstance(geneset_id, types.StringTypes):
                    geneset_id = ObjectId(geneset_id)
                else:  # 如果是其他类型，则报错
                    self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！')
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            with open(attributes_file, "rb") as f:
                for line in f:
                    line = line.strip().split("\t")
                    transitivity = line[0]
                    diameter = line[1]
                    average_shortest_path_length = line[2]
                    if diameter == "none":
                        diameter = 0
                    if average_shortest_path_length == "none":
                        average_shortest_path_length = 0
                    if transitivity == "none":
                        transitivity = 0
            if params != None:
                params = params
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else "NetworkCorfd_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                'transitivity': float(transitivity),
                'diameter': float(diameter),
                'min_dist': float(average_shortest_path_length),
            }
            try:
                collection = self.db['two_corr_network']
                network_corfd_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入network_corfd主表异常:{}'.format(e))
                self.bind_object.set_error("导入network_corfd主表异常")
        else:
            if not os.path.exists(attributes_file):
                self.bind_object.logger.error('attributes_file所指定的路径{}不存在，请检查！'.format(attributes_file))
                self.bind_object.set_error("attributes_file路径不存在")
            if main_table_id is None:
                self.bind_object.set_error("main为False时需提供main_table_id!")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！')
            data_list = []
            with open(attributes_file, "rb") as f:
                lines = f.readlines()
                for line in lines[1:]:
                    line = line.strip().split("\t")
                    transitivity = line[0]
                    diameter = line[1]
                    average_shortest_path_length = line[2]
                    if diameter == "none":
                        diameter = 0
                    if average_shortest_path_length == "none":
                        average_shortest_path_length = 0
                    if transitivity == "none":
                        transitivity = 0
                    insert_data = {
                        'corfd_id': main_table_id,
                        'transitivity': float(transitivity),
                        'diameter': float(diameter),
                        'min_dist': float(average_shortest_path_length),
                        "type": "single"
                    }
                    data_list.append(insert_data)
            try:
                collection = self.db['two_corr_network_node']
                collection.insert_many(data_list)
                network_corfd_id = main_table_id
            except Exception as e:
                self.bind_object.logger.error('更新network_corfd主表attributes_file出错:{}'.format(e))
                self.bind_object.set_error("更新主表attributes_file出错")
        return network_corfd_id

    @report_check
    def add_network_corfd_link(self, network_corfd_id, network_corfd_dir, anno_type="none"):
        if not isinstance(network_corfd_id, ObjectId):  # 检查传入的NetworkCor_id是否符合ObjectId类型
            if isinstance(network_corfd_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                network_corfd_id = ObjectId(network_corfd_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('network_corfd_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(network_corfd_dir):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('network_corfd_dir所指定的路径不存在，请检查！')
        data_list = []
        result = self.db['two_corr_network'].find_one({'_id': network_corfd_id})
        if not result:
            self.bind_object.set_error('找不到corr_network_by_cut.txt对应的主表id')
        self.id_names = {}
        if not os.path.exists(network_corfd_dir + "/corr_network_centrality.txt"):
            self.bind_object.set_error('corr_network_centrality所指定的路径不存在，请检查！')
        with open(network_corfd_dir + "/corr_network_centrality.txt", 'rb') as f1:
            f1.next()
            for line in f1:
                line = line.strip().split('\t')
                node_id = line[0]
                node_name = line[1]
                node_name = node_name.split("; ")[-1]
                self.id_names[node_name] = node_id
        with open(network_corfd_dir + "/corr_network_by_cut.txt", 'rb') as f:
            head1 = f.next()
            head2 = f.next()
            head = f.next()
            if "Node1_Name" in head:
                heads = head.strip().split("\t")
            else:
                self.bind_object.set_error('corr_network_by_cut.txt文件错误！')
            for line in f:
                line = line.strip().split('\t')
                source1 = line[0]
                target1 = line[1]
                source1 = source1.split("; ")[-1]
                target1 = target1.split("; ")[-1]
                source = int(self.id_names[source1]) - 1
                target = int(self.id_names[target1]) - 1
                # weight = line[2]
                Coefficient = line[2]
                insert_data = {
                    'corfd_id': network_corfd_id,
                    'source': int(source),
                    'target': int(target),
                    'coefficient': float(Coefficient),
                    "p_value" : float(line[3]),
                }
                data_list.append(insert_data)
            try:
                collection = self.db['two_corr_network_link']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入link表格%s信息出错:%s" % (network_corfd_dir, e))
                self.bind_object.set_error("导入network_corfd_dir信息出错")
            else:
                self.bind_object.logger.info("导入link表格%s信息成功!" % network_corfd_dir)

    def add_network_corfd_node(self, network_corfd_id, network_corfd_dir, profile1, profile2,color_level_id):
        if not isinstance(network_corfd_id, ObjectId):  # 检查传入的network_corfd_id是否符合ObjectId类型
            if isinstance(network_corfd_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                network_corfd_id = ObjectId(network_corfd_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('network_corfd_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(network_corfd_dir):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('network_corfd_dir所指定的路径不存在，请检查！')
        data_list = []
        result = self.db['two_corr_network'].find_one({'_id': network_corfd_id})
        if not result:
            self.bind_object.set_error('找不到real_node_table.txt对应的主表id')
        node_ids = {}
        node_de = {}
        node_clo = {}
        node_be = {}
        with open(network_corfd_dir + "/corr_network_centrality.txt", 'rb') as f1:
            head = f1.next()
            for line in f1:
                line = line.strip().split('\t')
                node_id = line[0]
                node_name = line[1]
                node_name = node_name.split("|")[0].strip()
                Degree_Centrality = line[2]
                Closeness_Centrality = line[3]
                Betweenness_Centrality = line[4]
                node_ids[node_name] = node_id
                node_de[node_name] = float(Degree_Centrality)
                node_clo[node_name] = float(Closeness_Centrality)
                node_be[node_name] = float(Betweenness_Centrality)
        name_abu1 = {}
        name_abu2 = {}
        nr_names = {}
        table1 = pd.read_table(profile1,sep="\t",header=0,index_col=0)
        table1['sum'] = table1.apply(lambda x: x.sum(), axis=1)
        table1['tmp'] = table1.index
        sel1 = pd.DataFrame(table1['tmp'].astype('str').str.split("|",expand=True)).iloc[:,0]
        table1["name"] = sel1
        need_table = table1[["name","sum"]]
        need_table2 = table1[["name"]]
        name_abu1 = need_table.set_index('name').T.to_dict('records')[0]
        self.bind_object.logger.info("-------------------")
        self.bind_object.logger.info(profile2)
        level2 = "env"
        if level2 == "env":
            table2 = pd.read_table(profile2,sep="\t",header=0,index_col=0)
            table2 = table2.T
        else:
             table2 = pd.read_table(profile2,sep="\t",header=0,index_col=0)
        table2['sum'] = table2.apply(lambda x: x.sum(), axis=1)
        table2['tmp'] = table2.index
        sel2 = pd.DataFrame(table2['tmp'].astype('str').str.split("|",expand=True)).iloc[:,0]
        table2["name"] = sel2
        need_table_2 = table2[["name","sum"]]
        need_table2_2 = table2[["name"]]
        all_names = pd.concat([need_table2,need_table2_2])
        all_names.index = all_names["name"]
        name_abu2 = need_table_2.set_index('name').T.to_dict('records')[0]
        nr_names =  all_names.set_index(all_names.index).T.to_dict('records')[0]

        name_clustering = {}
        with open(network_corfd_dir + "/corr_network_clustering.txt", 'rb') as f3:
            head = f3.next()
            for line in f3:
                line = line.strip().split('\t')
                node_name = line[1]
                node_name = node_name.split("|")[0].strip()
                Clustering = line[2]
                name_clustering[node_name] = Clustering

        with open(network_corfd_dir + "/corr_network_node_degree.txt", 'rb') as f:
            head = f.next()
            for line in f:
                line = line.strip().split('\t')
                node_ID = line[0]
                node_name_ori = line[1]
                #node_name = node_name_ori.split(";")[-1].strip()
                node_name = node_name_ori.split("|")[0].strip()
                #if len(node_name_ori.split("|")) > 1:
                sp_names = node_name_ori.split(";")
                if len(sp_names) > 1:
                    color_level = sp_names[color_level_id-1].strip()
                else:
                    color_level = sp_names[0]

                degree = line[2]
                degree_centrality = node_de[node_name]
                closeness_centrality = node_clo[node_name]
                betweenness_centrality = node_be[node_name]
                if name_abu1.has_key(node_name):
                    abu = name_abu1[node_name]

                elif name_abu2.has_key(node_name):
                    abu = name_abu2[node_name]
                else:
                    self.bind_object.logger.info("-------------------")
                    self.bind_object.logger.info(node_name)
                    self.bind_object.logger.info(name_abu1)
                    self.bind_object.logger.info(name_abu2)
                self.bind_object.logger.info(name_clustering)
                self.bind_object.logger.info(node_name)
                self.bind_object.logger.info(network_corfd_dir + "/corr_network_node_degree.txt")
                clus = name_clustering[node_name]
                node_data = {
                    "node": node_name.split('; ')[-1].strip(),
                    'corfd_id': network_corfd_id,
                    "type": "node",
                    'degree': int(degree),
                    'clustering': float(clus),
                    'group': color_level,
                    "abu": float(abu),
                    "rank": int(self.id_names[node_name.split('; ')[-1].strip()]) - 1,
                }
                insert_data = {
                    'corfd_id': network_corfd_id,
                    'node': node_name.split('; ')[-1].strip(),
                    'node_id': node_ID,
                    'degree_c': float(degree_centrality),
                    'closeness_c': float(closeness_centrality),
                    'betweenness_c': float(betweenness_centrality),
                    "type": "centrality"
                }
                data_list.append(node_data)
                data_list.append(insert_data)
            try:
                collection = self.db['two_corr_network_node']
                collection.insert_many(data_list)
                main_collection = self.db['two_corr_network']
                settled_params = {'software': "R-3.3.1 (stat),python-2.7 (stat)"}
                settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
                single_table_data = {"table_data": ["transitivity", "diameter", "min_dist"], "condition":{"type": "single"}}
                single_table_data_json = json.dumps(single_table_data, sort_keys=True, separators=(',', ':'))
                node_table_data = {"table_data": ["node", "degree", "clustering"], "condition":{"type": "node"}}
                node_table_data_json = json.dumps(node_table_data, sort_keys=True, separators=(',', ':'))
                centrality_table_data = {"table_data": ['node_id', 'node', "degree_c", 'closeness_c','betweenness_c'], "condition":{"type": "centrality"}}
                centrality_table_data_json = json.dumps(centrality_table_data, sort_keys=True, separators=(',', ':'))
                centrality_line_data = {"line_data":{"name":"specimen","data": {"x":"node_id", "y":["degree_c", "closeness_c", "betweenness_c"]} ,"condition":{"type": "centrality"}}}
                centrality_line_data_json = json.dumps(centrality_line_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id": ObjectId(network_corfd_id)},
                                                {"$set": {"status": 'end',
                                                    "main_id": network_corfd_id,
                                                    "settled_params": settled_params_json,
                                                    "single_table_data": single_table_data_json,
                                                    "node_table_data": node_table_data_json,
                                                    "centrality_line_data": centrality_line_data_json,
                                                    "centrality_table_data":centrality_table_data_json}})

            except Exception as e:
                self.bind_object.logger.error("导入node表格%s信息出错:%s" % (network_corfd_dir, e))
                self.bind_object.set_error("导入network_corfd_dir信息出错")
            else:
                self.bind_object.logger.info("导入node表格%s信息成功!" % network_corfd_dir)

    def add_network_corfd_degree(self, network_corfd_id, network_corfd_dir):
        if not isinstance(network_corfd_id, ObjectId):  # 检查传入的network_corfd_id是否符合ObjectId类型
            if isinstance(network_corfd_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                network_corfd_id = ObjectId(network_corfd_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('network_corfd_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(network_corfd_dir):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('network_corfd_dir所指定的路径不存在，请检查！')
        print network_corfd_id
        data_list = []
        result = self.db['two_corr_network'].find_one({'_id': network_corfd_id})
        if not result:
            self.bind_object.set_error('找不到network_corfd_nog对应的主表id')
        else:
            task_id = result['task_id']
            #samples_dic = id2name(task_id, type="task")
        with open(network_corfd_dir + "/corr_network_degree_distribution.txt", 'rb') as f:
            head = f.next()
            if "Degree" in head:
                heads = head.strip().split("\t")
            else:
                self.bind_object.set_error('corr_network_degree_distribution.txt文件错误！')
            for line in f:
                line = line.strip().split('\t')
                degree = line[0]
                count = line[1]
                insert_data = {
                    'corfd_id': network_corfd_id,
                    'x': int(degree),
                    'y': int(count),
                    "name": ""
                    #'type':2
                }
                data_list.append(insert_data)
            try:
                collection = self.db['two_corr_network_degree']
                collection.insert_many(data_list)
                stat_node_table_data = {"table_data": ["x", "y",]}
                stat_node_table_data_json = json.dumps(stat_node_table_data, sort_keys=True, separators=(',', ':'))
                stat_line_data = {"line_data": {"name":"name"}}
                stat_line_datajson = json.dumps(stat_line_data, sort_keys=True, separators=(',', ':'))
                main_collection = self.db['two_corr_network']
                main_collection.update({"_id":network_corfd_id},{"$set":{"main_id": network_corfd_id,
                                                                         "stat_line_data": stat_line_datajson,
                                                                        "stat_node_table_data": stat_node_table_data_json }})

            except Exception as e:
                self.bind_object.logger.error("导入degree表格%s信息出错:%s" % (network_corfd_dir, e))
                self.bind_object.set_error("导入network_corfd_dir信息出错")
            else:
                self.bind_object.logger.info("导入degree表格%s信息成功!" % network_corfd_dir)



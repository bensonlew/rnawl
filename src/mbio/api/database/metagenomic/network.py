# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20180330
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metagenomic.id_convert import name2id
from mbio.packages.metagenomic.id_convert import id2name
import json
import pandas as pd

class Network(Base):
    def __init__(self, bind_object):
        super(Network, self).__init__(bind_object)
        # self._db_name = Config().MONGODB + '_metagenomic'
        self._project_type = "metagenomic"

    @report_check
    def add_network(self, anno_type, attributes_file, geneset_id=None, group_detail=None, main=True,
                    main_table_id=None, params=None, name = None):
        if not os.path.exists(attributes_file):
            self.bind_object.logger.error('attributes_file所指定的路径{}不存在，请检查！'.format(attributes_file))
            self.bind_object.set_error('attributes_file路径不存在，请检查！', code="52802501")
        if main:
            if not isinstance(geneset_id, ObjectId):
                if isinstance(geneset_id, types.StringTypes):
                    geneset_id = ObjectId(geneset_id)
                else:
                   self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52802502")
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            with open(attributes_file, "rb") as f:
                for line in f:
                    line = line.strip().split("\t")
                    if line[0] == "Transitivity":
                        transitivity = line[1]
                    if line[0] == "Diameter":
                        diameter = line[1]
                    if line[0] == "Average_shortest_path_length":
                        average_shortest_path_length = line[1]
            if params != None:
                params = params
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else "Network_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                'anno_type': anno_type,
                'transitivity': float(transitivity),
                'diameter': float(diameter),
                'average_shortest_path_length': float(average_shortest_path_length),
            }
            try:
                collection = self.db['network']
                network_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入network主表异常:{}'.format(e))
                self.bind_object.set_error("导入network主表异常", code="52802503")
        else:
            with open(attributes_file, "rb") as f:
                for line in f:
                    line = line.strip().split("\t")
                    if line[0] == "Transitivity":
                        transitivity = line[1]
                    if line[0] == "Diameter":
                        diameter = line[1]
                    if line[0] == "Average_shortest_path_length":
                        average_shortest_path_length = line[1]
            if main_table_id is None:
                self.bind_object.set_error("main为False时需提供main_table_id!", code="52802504")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52802505")
            try:
                self.db['network'].update_one({'_id': main_table_id}, {
                    '$set': {'transitivity': float(transitivity), 'diameter': float(diameter),
                             'average_shortest_path_length': float(average_shortest_path_length)}})
                network_id = main_table_id
            except Exception as e:
                self.bind_object.logger.error('更新network主表attributes_file出错:{}'.format(e))
                self.bind_object.set_error("更新network主表attributes_file出错", code="52802506")
        return network_id

    @report_check
    def add_network_link(self, network_id, network_dir, anno_type):
        if not isinstance(network_id, ObjectId):
            if isinstance(network_id, types.StringTypes):
                network_id = ObjectId(network_id)
            else:
                self.bind_object.set_error('network_id必须为ObjectId对象或其对应的字符串！', code="52802507")
        file1 = network_dir + "/network_centrality.txt"
        file2 = network_dir + "/network_edge_table.txt"
        if not os.path.exists(file1):
            self.bind_object.logger.error('network_centrality所指定的路径{}不存在，请检查！'.format(file1))
            self.bind_object.set_error("network_centrality路径不存在，请检查！", code="52802508")
        if not os.path.exists(file2):
            self.bind_object.logger.error('network_edge_table所指定的路径{}不存在，请检查！'.format(file2))
            self.bind_object.set_error("network_edge_table路径不存在，请检查！", code="52802509")
        data_list = []
        result = self.db['network'].find_one({'_id': network_id})
        if not result:
            self.bind_object.set_error('找不到network_link对应的主表id', code="52802510")
        else:
            task_id = result['task_id']
            # samples_dic samples_dic = name2id(task_id, type="task")
        id_names = {}
        with open(file1, 'rb') as f1:
            f1.next()
            for line in f1:
                line = line.strip().split('\t')
                node_id = line[0]
                node_name = line[1]
                if anno_type == "nr":
                    node_name = node_name.split(";")[-1]
                id_names[node_name] = node_id
        with open(file2, 'rb') as f:
            head = f.next()
            if "from" in head:
                heads = head.strip().split("\t")
            else:
                self.bind_object.set_error('network_edge_table.txt文件错误！', code="52802511")
            for line in f:
                line = line.strip().split('\t')
                source1 = line[0]
                target1 = line[1]
                eweight = line[2]
                if anno_type == "nr":
                    source1 = source1.split(";")[-1]
                    target1 = target1.split(";")[-1]
                source = int(id_names[source1]) - 1
                target = int(id_names[target1]) - 1
                insert_data = {
                    'network_id': network_id,
                    'source': int(source),
                    'target': int(target),
                    'eweight': float(eweight)
                }
                data_list.append(insert_data)
            try:
                collection = self.db['network_link']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入link表格%s信息出错:%s" % (network_dir, e))
                self.bind_object.set_error("导入link表格出错", code="52802512")
            else:
                self.bind_object.logger.info("导入link表格%s信息成功!" % network_dir)

    @report_check
    def add_network_node(self, network_id, group_file, network_dir, anno_type, level, abu_file):
        if not isinstance(network_id, ObjectId):
            if isinstance(network_id, types.StringTypes):
                network_id = ObjectId(network_id)
            else:
                self.bind_object.set_error('network_id必须为ObjectId对象或其对应的字符串！', code="52802507")
        if not os.path.exists(network_dir):
            self.bind_object.set_error('network_dir所指定的路径不存在，请检查！', code="52802513")
        data_list = []
        result = self.db['network'].find_one({'_id': network_id})
        if not result:
            self.bind_object.set_error('找不到network_nodes_table.txt对应的主表id', code="52802514")
        else:
            task_id = result['task_id']
            samples_dic = name2id(task_id, type="task")
            group_id = json.loads(result["params"])["group_id"]
        file1 = network_dir + "/network_centrality.txt"
        file2 = network_dir + "/network_nodes_table.txt"
        if not os.path.exists(file1):
            self.bind_object.logger.error('network_centrality所指定的路径{}不存在，请检查！'.format(file1))
            self.bind_object.set_error('network_centrality路径不存在，请检查！', code="52802515")
        if not os.path.exists(file2):
            self.bind_object.logger.error('network_nodes_table所指定的路径{}不存在，请检查！'.format(file2))
            self.bind_object.set_error('network_nodes_table路径不存在，请检查！', code="52802516")
        if not os.path.exists(group_file):
            self.bind_object.logger.error('group_file所指定的路径{}不存在，请检查！'.format(group_file))
            self.bind_object.set_error('group_file路径不存在，请检查！', code="52802517")
        if not os.path.exists(abu_file):
            self.bind_object.logger.error('abu_file所指定的路径{}不存在，请检查！'.format(abu_file))
            self.bind_object.set_error('abu_file路径不存在，请检查！', code="52802518")
        abu_dict = self.abu_save(abu_file)
        node_ids = {}
        node_de = {}
        node_clo = {}
        node_be = {}
        self.bind_object.logger.info("save Degree_Centrality、Closeness_Centrality、Betweenness_Centrality")
        with open(file1, 'rb') as f1:
            head = f1.next()
            for line in f1:
                line = line.strip().split('\t')
                node_id = line[0]
                node_name = line[1]
                Degree_Centrality = line[2]
                Closeness_Centrality = line[3]
                Betweenness_Centrality = line[4]
                node_ids[node_name] = node_id
                node_de[node_name] = Degree_Centrality
                node_clo[node_name] = Closeness_Centrality
                node_be[node_name] = Betweenness_Centrality
        self.bind_object.logger.info("save group！")
        groups = {}
        same = []
        self.bind_object.logger.info("group path %s!" % group_file)
        with open(group_file, 'rb') as f2:
            head = f2.next()
            for line in f2:
                line = line.strip().split('\t')
                name = line[0]
                group = line[1]
                groups[name] = group
                if name == group:
                    same.append(name)
        # if len(same) == len(groups.keys()):
        if group_id == "all" or group_id == "All" or group_id == "ALL":
            for each in groups.keys():
                groups[each] = "all"
        self.bind_object.logger.info("save node!")
        with open(file2, 'rb') as f:
            head = f.next()
            for line in f:
                line = line.strip().split('\t')
                node_name = line[0]
                node_type = line[2]
                if anno_type == "nr":
                    node_name = node_name.split(";")[-1]
                if node_type == "tax_fun_node":
                    type = 1
                    group = level
                elif node_type == "sample_node":
                    type = 2
                    if groups.has_key(node_name):
                        group = groups[node_name]
                    else:
                        group = node_name
                if abu_dict.has_key(node_name):
                    abu = abu_dict[node_name]
                else:
                    abu = max(abu_dict.values())*5
                degree = line[3]
                weighted_degree = line[4]
                node_ID = node_ids[node_name]
                degree_centrality = node_clo[node_name]
                closeness_centrality = node_de[node_name]
                betweenness_centrality = node_be[node_name]
                if samples_dic.has_key(node_name):
                    node_name = samples_dic[node_name]
                insert_data = {
                    'network_id': network_id,
                    'node_name': node_name,
                    'degree': float(degree),
                    'group': group,
                    'node_ID': int(node_ID),
                    "abu": float(abu),
                    'degree_centrality': float(degree_centrality),
                    'closeness_centrality': float(closeness_centrality),
                    'betweenness_centrality': float(betweenness_centrality),
                    'weighted_degree': float(weighted_degree),
                    "type": type
                }
                data_list.append(insert_data)
            try:
                collection = self.db['network_node']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入node表格%s信息出错:%s" % (network_dir, e))
                self.bind_object.set_error("导入node表格信息出错", code="52802519")
            else:
                self.bind_object.logger.info("导入node表格%s信息成功!" % network_dir)

    @report_check
    def add_network_degree(self, network_id, network_dir, anno_type):
        if not isinstance(network_id, ObjectId):
            if isinstance(network_id, types.StringTypes):
                network_id = ObjectId(network_id)
            else:
                self.bind_object.set_error('network_id必须为ObjectId对象或其对应的字符串！', code="52802507")
        if not os.path.exists(network_dir):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('network_dir所指定的路径不存在，请检查！', code="52802513")
        file1 = network_dir + "/network_sample_degree.txt"
        file2 = network_dir + "/network_tax_fun_degree.txt"
        file3 = network_dir + "/network_nodes_degree.txt"
        if not os.path.exists(file1):
            self.bind_object.logger.error('network_sample_degree所指定的路径{}不存在，请检查！'.format(file1))
            self.bind_object.set_error('network_sample_degree路径不存在，请检查！', code="52802520")
        if not os.path.exists(file2):
            self.bind_object.logger.error('network_tax_fun_degree所指定的路径{}不存在，请检查！'.format(file2))
            self.bind_object.set_error('network_tax_fun_degree路径不存在，请检查！', code="52802521")
        if not os.path.exists(file3):
            self.bind_object.logger.error('network_nodes_degree所指定的路径{}不存在，请检查！'.format(file3))
            self.bind_object.set_error('network_tax_fun_degree路径不存在，请检查！', code="52802522")
        data_list = []
        result = self.db['network'].find_one({'_id': network_id})
        if not result:
            self.bind_object.set_error('找不到network_nog对应的主表id', code="52802523")
        else:
            task_id = result['task_id']
            #samples_dic = id2name(task_id, type="task")
        self.bind_object.logger.info("save network_sample_degree")
        with open(file1, 'rb') as f:
            head1 = f.next()
            head = f.next()
            if "Degree" in head:
                heads = head.strip().split("\t")
            else:
                self.bind_object.set_error('network_sample_degree.txt文件错误！', code="52802524")
            for line in f:
                line = line.strip().split('\t')
                degree = line[0]
                count = line[1]
                insert_data = {
                    'network_id': network_id,
                    'degree': int(degree),
                    'count': int(count),
                    'type': 2
                }
                data_list.append(insert_data)
        self.bind_object.logger.info("network_tax_fun_degree")
        with open(file2, 'rb') as f2:
            head1 = f2.next()
            head = f2.next()
            if "Degree" in head:
                heads = head.strip().split("\t")
            else:
                self.bind_object.set_error('network_tax_fun_degree.txt文件错误！', code="52802525")
            for line in f2:
                line = line.strip().split('\t')
                degree = line[0]
                count = line[1]
                insert_data = {
                    'network_id': network_id,
                    'degree': int(degree),
                    'count': int(count),
                    'type': 1
                }
                data_list.append(insert_data)
        self.bind_object.logger.info("network_nodes_degree")
        with open(file3, 'rb') as f3:
            head1 = f3.next()
            head = f3.next()
            if "Degree" in head:
                heads = head.strip().split("\t")
            else:
                self.bind_object.set_error('network_nodes_degree.txt文件错误！', code="52802526")
            for line in f3:
                line = line.strip().split('\t')
                degree = line[0]
                count = line[1]
                insert_data = {
                    'network_id': network_id,
                    'degree': int(degree),
                    'count': int(count),
                    'type': 3
                }
                data_list.append(insert_data)
            try:
                collection = self.db['network_degree']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入all_degree表格%s信息出错:%s" % (e))
                self.bind_object.set_error("导入all_degree出错", code="52802527")
            else:
                self.bind_object.logger.info("导入all_degree表格%s信息成功!" % network_dir)

    @report_check
    def abu_save(self, abu_file):
        abu_dict ={}
        abu_table = pd.read_table(abu_file, sep='\t', header=0)
        sample_list = abu_table.columns[1:len(abu_table.columns)]
        abu_table['Total'] = abu_table.loc[:, sample_list].apply(lambda x: x.sum(),axis=1)
        data = abu_table['Total']
        data.index = abu_table.iloc[:,0]
        abu_dict = data.to_dict()
        return abu_dict
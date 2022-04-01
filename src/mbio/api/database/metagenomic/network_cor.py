# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20171114
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


class NetworkCor(Base):
    def __init__(self, bind_object):
        super(NetworkCor, self).__init__(bind_object)
        # self._db_name = Config().MONGODB + '_metagenomic'
        self._project_type = "metagenomic"

    @report_check
    def add_network_cor(self, anno_type, attributes_file, geneset_id=None, main=True, main_table_id=None, params=None,
                        name = None):
        if not os.path.exists(attributes_file):  # 检查要上传的数据表路径是否存在
           self.bind_object.set_error('attributes_file所指定的路径不存在，请检查！', code="52802601")
        if main:
            if not isinstance(geneset_id, ObjectId):
                if isinstance(geneset_id, types.StringTypes):
                    geneset_id = ObjectId(geneset_id)
                else:  # 如果是其他类型，则报错
                    self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52802602")
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
                'name': name if name else "NetworkCor_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                #'specimen': specimen,
                'anno_type': anno_type,
                'transitivity': float(transitivity),
                'diameter': float(diameter),
                'average_shortest_path_length': float(average_shortest_path_length),
            }
            try:
                collection = self.db['network_cor']
                network_cor_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入network_cor主表异常:{}'.format(e))
                self.bind_object.set_error("导入network_cor主表异常", code="52802603")
        else:
            if not os.path.exists(attributes_file):
                self.bind_object.logger.error('attributes_file所指定的路径{}不存在，请检查！'.format(attributes_file))
                self.bind_object.set_error("attributes_file路径不存在", code="52802604")
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
            if main_table_id is None:
                self.bind_object.set_error("main为False时需提供main_table_id!", code="52802605")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52802606")
            try:
                self.db['network_cor'].update_one({'_id': ObjectId(main_table_id)}, {
                    '$set': {'transitivity': float(transitivity), 'diameter': float(diameter),
                             'average_shortest_path_length': float(average_shortest_path_length)}})
                network_cor_id = main_table_id
            except Exception as e:
                self.bind_object.logger.error('更新network_cor主表attributes_file出错:{}'.format(e))
                self.bind_object.set_error("更新主表attributes_file出错", code="52802607")
        return network_cor_id

    @report_check
    def add_network_cor_link(self, network_cor_id, network_cor_dir, anno_type):
        if not isinstance(network_cor_id, ObjectId):  # 检查传入的NetworkCor_id是否符合ObjectId类型
            if isinstance(network_cor_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                network_cor_id = ObjectId(network_cor_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('network_cor_id必须为ObjectId对象或其对应的字符串！', code="52802608")
        if not os.path.exists(network_cor_dir):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('network_cor_dir所指定的路径不存在，请检查！', code="52802609")
        data_list = []
        result = self.db['network_cor'].find_one({'_id': network_cor_id})
        if not result:
            self.bind_object.set_error('找不到corr_network_by_cut.txt对应的主表id', code="52802610")
        else:
            task_id = result['task_id']
            samples_dic = name2id(task_id, type="task")
        id_names = {}
        if not os.path.exists(network_cor_dir + "/corr_network_centrality.txt"):
            self.bind_object.set_error('corr_network_centrality所指定的路径不存在，请检查！', code="52802611")
        with open(network_cor_dir + "/corr_network_centrality.txt", 'rb') as f1:
            f1.next()
            for line in f1:
                line = line.strip().split('\t')
                node_id = line[0]
                node_name = line[1]
                if anno_type == "nr":
                    node_name = node_name.split(";")[-1]
                else:
                    node_name = node_name.split("|")[0]
                id_names[node_name] = node_id
        with open(network_cor_dir + "/corr_network_by_cut.txt", 'rb') as f:
            head1 = f.next()
            head2 = f.next()
            head = f.next()
            if "Node1_Name" in head:
                heads = head.strip().split("\t")
            else:
                self.bind_object.set_error('corr_network_by_cut.txt.txt文件错误！', code="52802612")
            for line in f:
                line = line.strip().split('\t')
                source1 = line[0]
                target1 = line[1]
                if anno_type == "nr":
                    source1 = source1.split(";")[-1]
                    target1 = target1.split(";")[-1]
                else:
                    source1 = source1.split("|")[0]
                    target1 = target1.split("|")[0]
                source = int(id_names[source1]) - 1
                target = int(id_names[target1]) - 1
                # weight = line[2]
                Coefficient = line[2]
                insert_data = {
                    'network_cor_id': network_cor_id,
                    'source': int(source),
                    'target': int(target),
                    'coefficient': float(Coefficient)
                }
                data_list.append(insert_data)
            try:
                collection = self.db['network_cor_link']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入link表格%s信息出错:%s" % (network_cor_dir, e))
                self.bind_object.set_error("导入network_cor_dir信息出错", code="52802613")
            else:
                self.bind_object.logger.info("导入link表格%s信息成功!" % network_cor_dir)

    def add_network_cor_node(self, network_cor_id, profile,level ,network_cor_dir,color_level=None, anno_type = None):
        if not isinstance(network_cor_id, ObjectId):  # 检查传入的network_cor_id是否符合ObjectId类型
            if isinstance(network_cor_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                network_cor_id = ObjectId(network_cor_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('network_cor_id必须为ObjectId对象或其对应的字符串！', code="52802608")
        if not os.path.exists(network_cor_dir):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('network_cor_dir所指定的路径不存在，请检查！', code="52802609")
        data_list = []
        result = self.db['network_cor'].find_one({'_id': network_cor_id})
        if not result:
            self.bind_object.set_error('找不到real_node_table.txt对应的主表id', code="52802614")
        else:
            task_id = result['task_id']
            samples_dic = name2id(task_id, type="task")
        node_ids = {}
        node_de = {}
        node_clo = {}
        node_be = {}
        with open(network_cor_dir + "/corr_network_centrality.txt", 'rb') as f1:
            head = f1.next()
            for line in f1:
                line = line.strip().split('\t')
                node_id = line[0]
                node_name = line[1]
                if anno_type == "nr":
                    node_name = node_name.split(";")[-1]
                else:
                    node_name = node_name.split("|")[0]
                Degree_Centrality = line[2]
                Closeness_Centrality = line[3]
                Betweenness_Centrality = line[4]
                node_ids[node_name] = node_id
                node_de[node_name] = float(Degree_Centrality)
                node_clo[node_name] = float(Closeness_Centrality)
                node_be[node_name] = float(Betweenness_Centrality)
        name_abu = {}
        table = pd.read_table(profile,sep="\t",header=0,index_col=0)
        table['sum'] = table.apply(lambda x: x.sum(), axis=1)
        table['tmp'] = table.index
        if anno_type == "nr":
            sel = pd.DataFrame(table['tmp'].astype('str').str.split(";",expand=True)).iloc[:,-1]
        else:
            sel = pd.DataFrame(table['tmp'].astype('str').str.split("|",expand=True)).iloc[:,0]
        table["name"] = sel
        need_table = table[["name","sum"]]
        name_abu = need_table.set_index('name').T.to_dict('records')[0]
        self.bind_object.logger.info("test why >>>>>>>>>>>>>>>>>>>>>>>>>>..")
        self.bind_object.logger.info(name_abu)
        '''
        with open(profile, 'rb') as f2:
            head = f2.next()
            for line in f2:
                line = line.strip().split('\t')
                abu = line[-1]
                name = line[0]
                if anno_type == "nr":
                    name = name.split(";")[-1]
                else:
                    name = name.split("|")[0]
                name_abu[name] = abu
        '''
        name_clustering = {}
        with open(network_cor_dir + "/corr_network_clustering.txt", 'rb') as f3:
            head = f3.next()
            for line in f3:
                line = line.strip().split('\t')
                node_name = line[1]
                if anno_type == "nr":
                    node_name = node_name.split(";")[-1]
                else:
                    node_name = node_name.split("|")[0]
                Clustering = line[2]
                name_clustering[node_name] = Clustering
        self.bind_object.logger.info("NR color")
        all_levels = ["Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"]
        self.bind_object.logger.info(anno_type)
        if anno_type == "nr":
            nr_level_index = all_levels.index(color_level)
        with open(network_cor_dir + "/corr_network_node_degree.txt", 'rb') as f:
            head = f.next()
            for line in f:
                line = line.strip().split('\t')
                node_ID = line[0]
                node_name_ori = line[1]
                if anno_type == "nr":
                    group = node_name_ori.split(";")[nr_level_index]
                    node_name = node_name_ori.split(";")[-1]
                    group = group.split("__", 1)[1]
                else:
                    node_name = node_name_ori.split("|")[0]
                    if color_level and len(node_name_ori.split("|")) > 1:
                        group = node_name_ori.split("|")[1]
                    else:
                        group = level
                degree = line[2]
                degree_centrality = node_de[node_name]
                closeness_centrality = node_clo[node_name]
                betweenness_centrality = node_be[node_name]
                abu = name_abu[node_name]
                clus = name_clustering[node_name]
                insert_data = {
                    'network_cor_id': network_cor_id,
                    'node_name': node_name,
                    'degree': int(degree),
                    'group': group,
                    'node_ID': node_ID,
                    "abu": float(abu),
                    'degree_centrality': float(degree_centrality),
                    'closeness_centrality': float(closeness_centrality),
                    'betweenness_centrality': float(betweenness_centrality),
                    'clustering': float(clus)
                }
                data_list.append(insert_data)
            try:
                collection = self.db['network_cor_node']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入node表格%s信息出错:%s" % (network_cor_dir, e))
                self.bind_object.set_error("导入network_cor_dir信息出错", code="52802615")
            else:
                self.bind_object.logger.info("导入node表格%s信息成功!" % network_cor_dir)

    def add_network_cor_degree(self, network_cor_id, network_cor_dir, anno_type):
        if not isinstance(network_cor_id, ObjectId):  # 检查传入的network_cor_id是否符合ObjectId类型
            if isinstance(network_cor_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                network_cor_id = ObjectId(network_cor_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('network_cor_id必须为ObjectId对象或其对应的字符串！', code="52802608")
        if not os.path.exists(network_cor_dir):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('network_cor_dir所指定的路径不存在，请检查！', code="52802609")
        print network_cor_id
        data_list = []
        result = self.db['network_cor'].find_one({'_id': network_cor_id})
        if not result:
            self.bind_object.set_error('找不到network_cor_nog对应的主表id', code="52802616")
        else:
            task_id = result['task_id']
            #samples_dic = id2name(task_id, type="task")
        with open(network_cor_dir + "/corr_network_degree_distribution.txt", 'rb') as f:
            head = f.next()
            if "Degree" in head:
                heads = head.strip().split("\t")
            else:
                self.bind_object.set_error('corr_network_degree_distribution.txt文件错误！', code="52802617")
            for line in f:
                line = line.strip().split('\t')
                degree = line[0]
                count = line[1]
                insert_data = {
                    'network_cor_id': network_cor_id,
                    'degree': int(degree),
                    'count': int(count),
                    #'type':2
                }
                data_list.append(insert_data)
            try:
                collection = self.db['network_cor_degree']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入degree表格%s信息出错:%s" % (network_cor_dir, e))
                self.bind_object.set_error("导入network_cor_dir信息出错", code="52802618")
            else:
                self.bind_object.logger.info("导入degree表格%s信息成功!" % network_cor_dir)


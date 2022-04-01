# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
# last_modify: 20200701
import datetime
from types import StringTypes

from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from bson.son import SON
import pandas as pd
from api_base import ApiBase
import json


class CorrNetwork(ApiBase):
    def __init__(self, bind_object):
        super(CorrNetwork, self).__init__(bind_object)
        self._task_id = self.bind_object._sheet.task_id
        self._name = self.bind_object.name
        self._work_dir = self.bind_object.work_dir
        self._output_dir = self.bind_object.output_dir

        self._main_table_id = None

    def add(self, main_id, main_table, file_list):
        '''
        :param file_list: 固定的文件和顺序
                          [corr_network_centrality.txt, cabundance_rank.txt, orr_network_node_degree.txt, corr_network_clustering.txt, corr_network_by_cut.txt, shared.0.03.pearson.corr]
        '''
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        name2id = self.add_node(file_list[1:4], main_id)
        self.add_link(main_id, file_list[4:], name2id)
        self.add_centrality(main_id, file_list[0])

        self.update_main('corr_network', main_id)

    def add_node(self, files, main_id):
        degree = pd.read_csv(files[1], header=0, sep='\t')
        col_d = {c: c.lower() for c in degree.columns}
        degree = degree.rename(columns=col_d)
        
        clustering = pd.read_csv(files[2], header=0, sep='\t')
        col_d = {c: c.lower() for c in clustering.columns}
        clustering = clustering.rename(columns=col_d)
        
        node = pd.merge(degree, clustering, on=['node_name', 'node_id'])

        rank = pd.read_csv(files[0], header=None, sep='\t', names=['node_name', 'abundance_rank'])
        node = pd.merge(node, rank, on='node_name')
        node['type'] = 'node'
        node['network_id'] = main_id
        node_data = self.df_to_mongo(node)
        self._run_add('corr_network_node', node_data)

        name2id = dict(zip(node['node_name'], node['node_id']))
        return name2id

    def add_link(self, main_id, files, name2id):
        link = pd.read_csv(files[0], header=None, sep='\t', skiprows=3, names=['source', 'target', 'value'])
        pc = pd.read_csv(files[1], header=None, sep='\t')
        pvalues = {'_'.join([r[1], r[2]]): r[-1] for r in pc.itertuples()}
        link['p_value'] = link.apply(lambda x: pvalues['_'.join([x[0], x[1]])] if '_'.join([x[0], x[1]]) in pvalues else pvalues['_'.join([x[1], x[0]])], axis=1)
        link['name'] = link['source'] + '_' + link['target']
        link['type'] = 'link'
        link['network_id'] = main_id
        link['source'] = link['source'].map(lambda x: name2id[x])
        link['target'] = link['target'].map(lambda x: name2id[x])
        link_data = self.df_to_mongo(link)
        self._run_add('corr_network_detail', link_data)

        #link.to_csv('corr_network_links.txt', sep='\t', index=False)

    def add_centrality(self, main_id, flpath):
        centrality = pd.read_csv(flpath, header=0, sep='\t')
        centrality['network_id'] = main_id
        centrality_data = self.df_to_mongo(centrality)
        self._run_add('corr_network_centrality', centrality_data)


    def df_to_mongo(self, df):
        keys = map(lambda x: x.lower(), df.columns)
        mongo_data = []
        df.apply(
            lambda x: mongo_data.append(SON(dict(zip(keys, x)))), axis=1
            )
        return mongo_data

    def update_main(self, table, table_id):
        tb = self.db[table]
        
        link_data = {
                'name': 'name', 'value': 'value',
                'condition': {'type': 'link'}
                }
        node_data = {
                'name': 'node_name',
                'condition': {'type': 'node'} 
                }
        node_column = [
                {'field': 'node_name', 'title': 'node_name', 'filter': 'true', 'sort': 'true', 'type': 'string'},
                {'field': 'degree', 'title': 'degree', 'filter': 'true', 'sort': 'true', 'type': 'int'},
                {'field': 'clustering', 'title': 'clustering', 'filter': 'true', 'sort': 'true', 'type': 'float'},
                {'field': 'node_id', 'title': 'node_id', 'filter': 'true', 'sort': 'true', 'type': 'string'},
                ]
        #node_table = {'column': json.dumps(node_column), 'condition': {}}
        node_table = {'column': node_column, 'condition': {}}

        centrality_column = [
                {'field': 'node_id', 'title': 'Node_ID','filter':'true', 'sort': 'true', 'type': 'string'},
                {'field': 'node_name', 'title': 'Node_Name','filter':'true', 'sort': 'true', 'type': 'string'},
                {'field': 'degree_centrality', 'title': 'Degree_Centrality','filter':'true', 'sort': 'true', 'type': 'float'},
                {'field': 'closeness_centrality', 'title': 'Closeness_Centrality','filter':'true', 'sort': 'true', 'type': 'float'},
                {'field': 'betweenness_centrality', 'title': 'Betweenness_Centrality','filter':'true', 'sort': 'true', 'type': 'float'},
                ]
        #centrality_table = {'column': json.dumps(centrality_column), 'condition': {}}
        centrality_table = {'column': centrality_column, 'condition': {}}

        updata_info = {
                "main_id": table_id,
                "link_data": json.dumps(link_data),
                "node_data": json.dumps(node_data),
                "node_table_data": json.dumps(node_table),
                "centrality_table_data": json.dumps(centrality_table),
                }

        tb.update({'_id': table_id},
                  {'$set': updata_info})
        self.bind_object.logger.info('更新主表成功！')

    def _run_add(self, name, data, main=False):
        try:
            collection = self.db[name]
            collection.insert_many(data)
        except Exception, e:
            self.bind_object.set_error('导入表%s出错：%s') % (name, e)
        else:
            self.bind_object.logger.info('导入表%s成功！' % name)

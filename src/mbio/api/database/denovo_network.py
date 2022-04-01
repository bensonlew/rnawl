# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# last_modify:20161027
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import bson.binary
from cStringIO import StringIO
import json
import re
from gridfs import GridFS
class DenovoNetwork(Base):
    def __init__(self, bind_object):
        super(DenovoNetwork, self).__init__(bind_object)
        self._db_name = Config().MONGODB + '_rna'

    @report_check
    def add_network(self, params, softpower_png, module_png, module=None, name=None):
        """
        if not isinstance(express_id, ObjectId):
            if isinstance(express_id, types.StringTypes):
                express_id = ObjectId(express_id)
            else:
                raise Exception('express_matrix_id必须为ObjectId对象或其对应的字符串！')
        """
        if not os.path.exists(softpower_png):
            raise Exception('softpower_png所指定的路径:{}不存在，请检查！'.format(softpower_png))
        if not os.path.exists(module_png):
            raise Exception('module_png所指定的路径:{}不存在，请检查！'.format(module_png))
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        collection = self.db['sg_denovo_network']
        softpower_fs = GridFS(self._db_name)
        softpower_id = softpower_fs.put(open(softpower_png, 'r'))
        module_fs = GridFS(self._db_name_)
        module_id = module_fs.put(open(module_png, 'r'))
        #params['diff_fpkm'] = str(express_id)
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'network_table_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'desc': '差异基因网络分析主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'status': 'end',
            'softpower_png': ObjectId(softpower_id),
            'module_png': ObjectId(module_png),
            'module': module
        }
        network_id = collection.insert_one(insert_data).inserted_id
        return network_id

    @report_check
    def add_network_detail(self, network_id, node_path, edge_path):
        if not isinstance(network_id, ObjectId):
            if isinstance(network_id, types.StringTypes):
                network_id = ObjectId(network_id)
        if not os.path.exists(node_path):
            raise Exception('node_path所指定的路径:{}不存在，请检查！'.format(node_path))
        if not os.path.exists(edge_path):
            raise Exception('edge_path所指定的路径:{}不存在，请检查！'.format(edge_path))
        data_list = []
        gene_color = {}
        with open(node_path, 'rb') as n, open(edge_path, 'rb') as f:
            n.readline()
            for line in n:
                line = line.strip().split('\t')
                gene_color[line[0]] = line[2] #提取all_nodes.txt文件第1、3列
            f.readline()
            for line in f:
                line = line.strip().split('\t')
                weight_value = line[2]
                if re.search(r'e',weight_value):
                     float_num = re.search(r'e-0(\d)',weight_value).group(1)
                     new_weight_value = round(float(weight_value),int(float_num) + 3)
                else:
                     new_weight_value = round(float(weight_value),3)
                data = [
                    ('network_id', network_id),
                    ('gene_id1', {'name': line[0], 'color': gene_color[line[0]]}),
                    ('gene_id2', {'name': line[1], 'color': gene_color[line[1]]}),
                    ('weight', new_weight_value), #提取all_edges.txt文件前三列
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_denovo_network_detail"]
            collection.insert_many(data_list)
            print collection
        except Exception, e:
            self.bind_object.logger.error("导入网络表达统计表：%s，%s信息出错:%s" % (node_path, edge_path, e))
        else:
            self.bind_object.logger.info("导入网络表达统计表:%s， %s信息成功!" % (node_path, edge_path))
        
    @report_check
    def add_network_module(self, network_id, module_path, module_color):
        if not isinstance(network_id, ObjectId):
            if isinstance(network_id, types.StringTypes):
                network_id = ObjectId(network_id)
        if not os.path.exists(module_path):
            raise Exception('module_path所指定的路径:{}不存在，请检查！'.format(module_path))
        data_list = []
        with open(module_path, 'rb') as f:
            f.readline()
            for line in f:
                line = line.strip().split('\t')
                weight_value = line[2]
                if re.search(r'e',weight_value):
                     float_num = re.search(r'e-0(\d)',weight_value).group(1)
                     new_weight_value = round(float(weight_value),int(float_num) + 3)
                else:
                     new_weight_value = round(float(weight_value),3) #保留三位有效的数字比如4.234345345e-04变为4.234e-04
                data = [
                    ('network_id', network_id),
                    ('gene_id1', line[0]),
                    ('gene_id2', line[1]),
                    ('weight', new_weight_value),
                    ('module_color', module_color),
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_denovo_network_module"]
            collection.insert_many(data_list)
            print collection
        except Exception, e:
            self.bind_object.logger.error("导入网络表达统计表：%s信息出错:%s" % (module_path, e))
        else:
            self.bind_object.logger.info("导入网络表达统计表:%s信息成功!" % module_path)

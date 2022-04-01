# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'
# last_modify:20211124
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
from bson.son import SON
from bson.objectid import ObjectId
import pandas as pd
from mbio.packages.metabolome.common import check_metab


class RelationCorrNetwork(Base):
    def __init__(self, bind_object):
        super(RelationCorrNetwork, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_relation_corr_network(self, name=None, main_id=None, params =None):
        if not main_id:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'AssoCorrNetwork_Origin',
                'params': params if params else '',
                'status': 'end',
                'main_id': ''
            }
            try:
                collection = self.db['relation_corr_network']
                association_corr_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", association_corr_id, association_corr_id)
            except Exception, e:
                self.bind_object.set_error('导入relation_corr_network主表异常:%s', variables=(e), code="54700101")
        else:
            self.update_table("main_id", main_id, main_id)
            association_corr_id = main_id
        return association_corr_id

    @report_check
    def add_relation_corr_network_detail(self, relation_corrnetwork_id, result_file=None, gene_id2name=None):
        if gene_id2name:
            df_geneid2name = pd.read_table(gene_id2name,"\t")
            geneid2name = df_geneid2name.fillna("-")
            dict_gene2name = geneid2name.set_index(["gene_id"])["gene_name"].to_dict()
        data_list = []
        with open(result_file, 'r') as result1:
            list1 = result1.readlines()
            for i in list1[1:]:
                line = i.strip().split('\t')
                insert_data = {
                    'corr_id': relation_corrnetwork_id,
                    'metabolite':line[0],
                    'gene_id':line[1],
                    'pvalue':float(line[3]),
                    'qvalue':float(line[4]),
                    'corr':float(line[2])
                }
                if gene_id2name:
                    insert_data["gene_name"] = dict_gene2name[line[1]]
                else:
                    insert_data["gene_name"] = ""
                data_son = SON(insert_data)
                data_list.append(data_son)
        try:
            collection = self.db['relation_corr_network_detail']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(result_file, e), code="54700107")
        else:
            pass
        self.bind_object.logger.info("导入表格%s信息成功!" % result_file)

    @report_check
    def update_table(self, str, name, main_table_id):
        try:
            self.db['relation_corr_network'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {str: name}})
        except Exception as e:
            self.bind_object.set_error('relation_corr_network%s字段出错:%s', variables=(str,e), code="54700109")

    @report_check
    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type), code="54700110")
        return object_id

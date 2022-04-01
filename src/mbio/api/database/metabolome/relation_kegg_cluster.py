# -*- coding: utf-8 -*-

from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import pandas as pd
from bson.objectid import ObjectId
from mbio.packages.metabolome.common import check_metab


class RelationKeggCluster(Base):
    def __init__(self, bind_object):
        super(RelationKeggCluster, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_relation_kegg_heatmap(self, name=None, main_id=None, params=None, pathway_tree=None, set_tree=None):
        if not main_id:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'RelationKeggCluster_Origin',
                'params': params if params else '',
                'status': 'end',
                'pathway_tree': '',
                'set_tree': '',
                "setname": '',
                'main_id': ''
            }
            try:
                collection = self.db['relation_kegg_heatmap']
                relation_heatmap_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", relation_heatmap_id, relation_heatmap_id)
            except Exception, e:
                self.bind_object.set_error('导入relation_kegg_heatmap主表异常:%s', variables=(e), code="54700101")
        else:
            self.bind_object.logger.info("main_id1为{}".format(main_id))
            self.update_table("main_id", main_id, main_id)
            relation_heatmap_id = main_id
            self.bind_object.logger.info("relation_heatmap_id0为{}".format(relation_heatmap_id))
        if pathway_tree:
            if not os.path.exists(pathway_tree):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pathway_tree), code="54700102")
            with open(pathway_tree, "r") as f:
                pathwaytree = f.readline().strip()
            self.update_table("pathway_tree", pathwaytree, relation_heatmap_id)
        if set_tree:
            if not os.path.exists(set_tree):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(set_tree), code="54700103")
            with open(set_tree, "r") as f2:
                settree = f2.readline().strip()
            self.update_table("set_tree", settree, relation_heatmap_id)
        return relation_heatmap_id

    def add_relation_kegg_heatmap_detail(self, relation_heatmap_id, pathway_table=None, exp_table=None, pvalue_table=None,id2db=None,id2term=None,select=None,species=None):

        self.bind_object.logger.info("relation_heatmap_id1为{}".format(relation_heatmap_id))
        relation_heatmap_id = self.check_id(relation_heatmap_id, "heatmap_id")
        self.bind_object.logger.info("relation_heatmap_id2为{}".format(relation_heatmap_id))
        result = self.db['relation_kegg_heatmap'].find_one({'_id': relation_heatmap_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(relation_heatmap_id), code="54700806")
        else:
            task_id = result['task_id']
        result = self.db['relation_kegg_heatmap'].find_one({'_id': relation_heatmap_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(relation_heatmap_id), code="54700106")
        else:
            task_id = result['task_id']
        df_pathway_cluster = pd.read_table(pathway_table, '\t', header=None)
        # 获取pathway和聚类层级对应关系
        cluster_list = []
        pathway_list = []
        for i in df_pathway_cluster.index:
            pathway = df_pathway_cluster.loc[i, 1].split(';')
            for j in pathway:
                pathway_list.append(j)
                cluster_list.append(df_pathway_cluster.loc[i, 0])
        dict_p2cluster = dict(zip(pathway_list,cluster_list))
        self.bind_object.logger.info("字典信息为{}".format(dict_p2cluster))
        df_id2db = pd.read_table(id2db, '\t')
        dict_id2db = df_id2db.set_index(["ID"])["Database"].to_dict()
        df_id2term = pd.read_table(id2term, '\t')
        dict_id2term = df_id2term.set_index(["ID"])["Term"].to_dict()
        df_exp = pd.read_table(exp_table, '\t')
        df_pvalue = pd.read_table(pvalue_table, '\t')
        sample_list = df_exp.columns.tolist()[1:]
        set_name = []
        num = 0
        for name in sample_list:
            num += 1
            set_name.append((name, "set"+str(num)))
        set_dict = dict(set_name)

        data_list = []
        for i in df_exp.index:
            pathway_id = df_exp.loc[i, "ID"]
            insert_data = {
                "heatmap_id": relation_heatmap_id,
                "pathway_id": pathway_id,
                "database": dict_id2db[df_exp.loc[i, "ID"]],
                "pathway_desc": dict_id2term[df_exp.loc[i, "ID"]],
                "ncluster": dict_p2cluster[df_exp.loc[i, "ID"]],
                "type": select
            }
            for j in sample_list:
                insert_data[set_dict[j]+"_value"] = df_exp.loc[i, j]
                insert_data[set_dict[j]+"_pvalue"] = df_pvalue.loc[i, j]
            data_list.append(insert_data)
        try:
            collection = self.db['relation_kegg_heatmap_detail']
            collection.insert_many(data_list)
            self.bind_object.logger.info("relation_heatmap_id3为{}".format(relation_heatmap_id))
            self.update_table('setname', set_dict, relation_heatmap_id)
        except Exception as e:
            self.bind_object.set_error("导入表格relation_kegg_heatmap信息出错:%s", variables=(e), code="54700107")
        else:
            pass
        self.bind_object.logger.info("导入表格信息成功!")

    @report_check
    def update_table(self, str, name, main_table_id):
        try:
            self.db['relation_kegg_heatmap'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {str: name}})
        except Exception as e:
            self.bind_object.set_error('relation_kegg_heatmap%s字段出错:%s', variables=(str, e), code="54700109")

    @report_check
    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type), code="54700110")
        return object_id

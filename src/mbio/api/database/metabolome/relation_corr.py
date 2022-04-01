# -*- coding: utf-8 -*-

from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import pandas as pd
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metabolome.common import check_metab


class RelationCorr(Base):
    def __init__(self, bind_object):
        super(RelationCorr, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_relation_corr(self,  name=None, main_id=None, params =None, trans_tree=None, metab_tree=None, gene_id2name=None, trans_select_table=None):
        if not main_id:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'RelationCorr_Origin',
                'params': params if params else '',
                'status': 'end',
                'gene_id': '',
                'gene_name': '',
                'trans_tree': '',
                'metab_tree': '',
                'main_id': ''
            }
            try:
                collection = self.db['relation_corr']
                relation_corr_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", relation_corr_id, relation_corr_id)
            except Exception, e:
                self.bind_object.set_error('导入association_corr主表异常:%s', variables=(e), code="54700101")
        else:
            self.update_table("main_id", main_id, main_id)
            relation_corr_id = main_id
        if trans_tree:
            if not os.path.exists(trans_tree):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(trans_tree), code="54700102")
            with open(trans_tree,"r") as f:
                trans_tree = f.readline().strip()
            self.update_table("trans_tree", trans_tree, relation_corr_id)
        if metab_tree:
            if not os.path.exists(metab_tree):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(metab_tree), code="54700103")
            with open(metab_tree,"r") as f2:
                metabtree = f2.readline().strip()
                metabtree = check_metab(metabtree)
            self.update_table("metab_tree", metabtree, relation_corr_id)
        
        trans_select_table = pd.read_table(trans_select_table,'\t')
        gene_id_list = trans_select_table["gene_id"].tolist()
        self.update_table("gene_id", gene_id_list, relation_corr_id)
        if gene_id2name:
            if not os.path.exists(gene_id2name):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(gene_id2name), code="54700104")
            gene_id2name = pd.read_table(gene_id2name,'\t')
            gene_id2name = gene_id2name.fillna("-")
            id2name_dict = gene_id2name.set_index(["gene_id"])["gene_name"].to_dict()
            trans_name_list = []
            for i in gene_id_list:
                trans_name_list.append(id2name_dict[i])
            self.update_table("gene_name", trans_name_list, relation_corr_id)
        else:
            empty_list = []
            self.update_table("gene_name", empty_list, relation_corr_id)
        return relation_corr_id

    @report_check
    def add_relation_corr_detail(self, relation_corr_id, corr_file, p_file, metab_anno=None, metab_hmdb_anno=None):
        relation_corr_id = self.check_id(relation_corr_id, "relation_corr_id")
        if not os.path.exists(corr_file):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(corr_file), code="54700105")
        result = self.db['relation_corr'].find_one({'_id': relation_corr_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(relation_corr_id), code="54700106")
        else:
            task_id = result['task_id']
        # 提取代谢物对应的注释信息
        metab_anno_table = pd.read_table(metab_anno,'\t')
        if metab_hmdb_anno:
            metab_hmdb_anno_table = pd.read_table(metab_hmdb_anno,'\t')
            metab_list = metab_hmdb_anno_table["metab"].tolist()
        data_list = []
        with open(corr_file, 'rb') as f:
            head = f.next()
            sams = head.strip().split("\t")
            self.bind_object.logger.info("sams为{}".format(sams))
            number1 = 1
            dict1 = {}
            for s in sams[1:]:
                dict1[s] = "key" + str(number1)
                number1 += 1
            self.bind_object.logger.info("dict1为{}".format(dict1))
            for line in f:
                line = line.strip().split('\t')
                name = check_metab(line[0])
                insert_data = {
                        'corr_id': relation_corr_id,
                        'metab':name,
                        'type': 'corr'
                        }
                for x in metab_anno_table.index:
                    if metab_anno_table.loc[x, "metab"] == name:
                        insert_data["kegg_category"] = metab_anno_table.loc[x, "compound_first_category"]
                        insert_data["metab_id"] = metab_anno_table.loc[x, "metab_id"]
                        
                if metab_hmdb_anno:
                    for y in metab_hmdb_anno_table.index:
                        if metab_hmdb_anno_table.loc[y, "metab"] == name:
                            insert_data["hmdb_superclass"] = metab_hmdb_anno_table.loc[y, "superclass"]
                            insert_data["hmdb_class"] = metab_hmdb_anno_table.loc[y, "class"]
                            insert_data["hmdb_subclass"] = metab_hmdb_anno_table.loc[y, "subclass"]
                    if name not in metab_list:
                        insert_data["hmdb_superclass"] = "-"
                        insert_data["hmdb_class"] = "-"
                        insert_data["hmdb_subclass"] = "-"
                else:
                    insert_data["hmdb_superclass"] = "-"
                    insert_data["hmdb_class"] = "-"
                    insert_data["hmdb_subclass"] = "-"
                for i in range(1, len(sams)):
                    sam_corr = float(line[i])
                    insert_data[dict1[sams[i]]] = sam_corr
                data_son1 = SON(insert_data)
                data_list.append(data_son1)

        with open(p_file, 'rb') as f:
            head = f.next()
            sams = head.strip().split("\t")
            number2 = 1
            dict2 = {}
            for s in sams[1:]:
                dict2[s] = "key" + str(number2)
                number2 += 1
            self.bind_object.logger.info("dict2为{}".format(dict2))
            for line in f:
                line = line.strip().split('\t')
                name = check_metab(line[0])
                insert_data1 = {
                    'corr_id': relation_corr_id,
                    'metab':name,
                    'type': 'pvalue'
                }
                for x in metab_anno_table.index:
                    if metab_anno_table.loc[x, "metab"] == name:
                        insert_data1["kegg_category"] = metab_anno_table.loc[x, "compound_first_category"]
                        insert_data1["metab_id"] = metab_anno_table.loc[x, "metab_id"]
                if metab_hmdb_anno:
                    for y in metab_hmdb_anno_table.index:
                        if metab_hmdb_anno_table.loc[y, "metab"] == name:
                            insert_data1["hmdb_superclass"] = metab_hmdb_anno_table.loc[y, "superclass"]
                            insert_data1["hmdb_class"] = metab_hmdb_anno_table.loc[y, "class"]
                            insert_data1["hmdb_subclass"] = metab_hmdb_anno_table.loc[y, "subclass"]
                    if name not in metab_list:
                        insert_data1["hmdb_superclass"] = "-"
                        insert_data1["hmdb_class"] = "-"
                        insert_data1["hmdb_subclass"] = "-"
                else:
                    insert_data1["hmdb_superclass"] = "-"
                    insert_data1["hmdb_class"] = "-"
                    insert_data1["hmdb_subclass"] = "-"                
                for i in range(1, len(sams)):
                    sam_corr = float(line[i])
                    insert_data1[dict2[sams[i]]] = sam_corr
                data_son2 = SON(insert_data1)
                data_list.append(data_son2)
        self.bind_object.logger.info(data_list)
        collection = self.db['relation_corr_detail']
        collection.insert_many(data_list)
        self.update_table('relate_geneid', dict1, relation_corr_id)
        self.bind_object.logger.info("导入表格relation_corr_detail信息成功!")

    @report_check
    def update_table(self, str, name, main_table_id):
        try:
            self.db['relation_corr'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {str: name}})
        except Exception as e:
            self.bind_object.set_error('relation_corr%s字段出错:%s', variables=(str,e), code="54700109")

    @report_check
    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type), code="54700110")
        return object_id

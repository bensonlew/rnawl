# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'
# last_modify:20211124
import pandas as pd
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
from bson.son import SON
from bson.objectid import ObjectId
from biocluster.config import Config


class RelationKeggpVenn(Base):
    def __init__(self, bind_object):
        super(RelationKeggpVenn, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_relation_keggp_venn(self, name=None, main_id=None, params =None):
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
                relation_keggpvenn_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", relation_keggpvenn_id, relation_keggpvenn_id)
            except Exception, e:
                self.bind_object.set_error('导入relation_corr_network主表异常:%s', variables=(e), code="54700101")
        else:
            self.update_table("main_id", main_id, main_id)
            relation_keggpvenn_id = main_id
        return relation_keggpvenn_id

    def add_relation_keggp_venn_detail(self, main_id, metab_desc, metab_keggp_table, trans_keggp_table, species_name, gene_name2id=None, version=None):
        if version in ['v94.2', '94.2']:
            ref_db = self.ref_db["kegg_v94.2_organisms"]
        elif version in ['kegg', '']:
            ref_db = self.ref_db["kegg_organisms"]
        else:
            ref_db = self.ref_db["kegg_v2021.09.18_organisms"]
        df_metab_desc = pd.read_table(metab_desc, '\t')
        df_metab_desc = df_metab_desc.fillna("-")
        dict_metab = df_metab_desc.set_index(["metab_id"])["Metabolite"].to_dict()
        df_metab_keggp_table = pd.read_table(metab_keggp_table, '\t')
        metab_pathway_list = df_metab_keggp_table["pathway_id"].tolist()
        set_metab = set(metab_pathway_list)
        if gene_name2id:
            df_gene_name2id = pd.read_table(gene_name2id, '\t')
            df_trans_desc = df_gene_name2id.fillna("-")
            dict_trans = df_trans_desc.set_index(["gene_id"])["gene_name"].to_dict()
        
        df_trans_keggp_table = pd.read_table(trans_keggp_table, '\t')
        trans_pathway_list = df_trans_keggp_table["pathway"].tolist()
        self.bind_object.logger.info("trans_pathway_list为{}".format(len(trans_pathway_list)))
        trans_newpathway_list = []
        for pathway in trans_pathway_list:
            if not pathway.startswith("map"):
                species = pathway[0:-5]
                rename = pathway.replace(species, "map")
                trans_newpathway_list.append(rename)
            else:
                trans_newpathway_list.append(pathway)
                species = "map"
        self.bind_object.logger.info("trans_newpathway_list为{}".format(len(trans_newpathway_list)))
        sp_species = species_name.split(";")
        self.bind_object.logger.info("sp_species为{}".format(sp_species))
        if len(sp_species) == 1:
            pass
        else:
            res = ref_db.find_one({"second_category": sp_species[1]})
            species = res["name"]
        self.bind_object.logger.info("species为{}".format(species))
        set_trans = set(trans_newpathway_list)
        set1 = set_metab - set_trans
        set2 = set_trans - set_metab
        set3 = set_metab & set_trans
        self.bind_object.logger.info("set1为{}".format(set1))
        self.bind_object.logger.info("set1为{}".format(len(set1)))
        self.bind_object.logger.info("set2为{}".format(set2))
        self.bind_object.logger.info("set2为{}".format(len(set2)))
        self.bind_object.logger.info("set3为{}".format(set3))
        self.bind_object.logger.info("set3为{}".format(len(set3)))
        list1 = []
        list2 = []
        list3 = []
        for p in set1:
            list1.append(p.replace(p[0:-5], species))
        for p in set2:
            list2.append(p.replace(p[0:-5], species))
        for p in set3:
            list3.append(p.replace(p[0:-5], species))
        insert_data1 = {
            'keggp_id': main_id,
            'pathway_list': list1,
            "pathway_num": len(list1),
            'type': "metabolome"
        }
        insert_data2 = {
            'keggp_id': main_id,
            'pathway_list': list2,
            "pathway_num": len(list2),
            'type': "transcript"
        }
        insert_data3 = {
            'keggp_id': main_id,
            'pathway_list': list3,
            "pathway_num": len(list3),
            'type': "transcript&metabolome"
        }
        try:
            collection = self.db['relation_keggp_venn']
            collection.insert(insert_data1, check_keys=False)
            collection.insert(insert_data2, check_keys=False)
            collection.insert(insert_data3, check_keys=False)
        except Exception as e:
            self.bind_object.set_error("导入表格relation_keggp_venn信息出错:%s",e)
        else:
            self.bind_object.logger.info("导入表格relation_keggp_venn信息成功!")

        data_list = []
        empty_list = []
        for j in df_metab_keggp_table.index:
            metab_id = df_metab_keggp_table.loc[j, "metab_id"].split(';')
            metab_list = []
            for x in metab_id:
                metab_list.append(dict_metab[x])
            if df_metab_keggp_table.loc[j, "pathway_id"] in set1:
                pathway_id = df_metab_keggp_table.loc[j, "pathway_id"]
                insert_data4 = {
                    'keggp_id': main_id,
                    "pathway_id": pathway_id.replace(pathway_id[0:-5], species),
                    "description": df_metab_keggp_table.loc[j, "description"],
                    "first_category": df_metab_keggp_table.loc[j, "first_category"],
                    "second_category":df_metab_keggp_table.loc[j, "second_category"],
                    "metabolite_list": metab_list,
                    "metabolite_id": metab_id,
                    "metab_num": df_metab_keggp_table.loc[j, "count"],
                    "gene_id": empty_list,
                    "gene_name": empty_list,
                    "gene_num": 0
                }
                data_son1 = SON(insert_data4)
                data_list.append(data_son1)
                
        for i in df_trans_keggp_table.index:
            pathway_id = df_trans_keggp_table.loc[i, "pathway"]
            map_id = pathway_id.replace(pathway_id[0:-5], "map")
            gene_id = df_trans_keggp_table.loc[i, "gene_id"].split(';')
            gene_name_list = []
            for z in gene_id:
                if gene_name2id:
                     gene_name_list.append(dict_trans[z])
                pass
            if map_id in set2:
                self.bind_object.logger.info("map_id为：{}".format(map_id))
                insert_data5 = {
                    'keggp_id': main_id,
                    "pathway_id": pathway_id.replace(pathway_id[0:-5], species),
                    "description": df_trans_keggp_table.loc[i, "description"],
                    "first_category": df_trans_keggp_table.loc[i, "first_category"],
                    "second_category": df_trans_keggp_table.loc[i, "second_category"],
                    "metabolite_list": empty_list,
                    "metabolite_id": empty_list,
                    "metab_num": 0,
                    "gene_id": gene_id,
                    "gene_name": gene_name_list,
                    "gene_num": df_trans_keggp_table.loc[i, "count"]
                }
                data_son2 = SON(insert_data5)
                data_list.append(data_son2)
                
        for k in set3:
            for j in df_metab_keggp_table.index:
                metab_id = df_metab_keggp_table.loc[j, "metab_id"].split(';')
                metab_list = []
                for x in metab_id:
                    metab_list.append(dict_metab[x])
                if df_metab_keggp_table.loc[j, "pathway_id"] == k:
                    pathway_id = df_metab_keggp_table.loc[j, "pathway_id"]
                    insert_data6 = {
                            'keggp_id': main_id,
                            "pathway_id": pathway_id.replace(pathway_id[0:-5], species),
                            "description": df_metab_keggp_table.loc[j, "description"],
                            "first_category": df_metab_keggp_table.loc[j, "first_category"],
                            "second_category": df_metab_keggp_table.loc[j, "second_category"],
                            "metabolite_list": metab_list,
                            "metabolite_id": metab_id,
                            "metab_num": df_metab_keggp_table.loc[j, "count"]}
            for i in df_trans_keggp_table.index:
                pathway_id = df_trans_keggp_table.loc[i, "pathway"]
                map_id = pathway_id.replace(pathway_id[0:-5], "map")
                gene_id = df_trans_keggp_table.loc[i, "gene_id"].split(';')
                gene_name_list = []
                for z in gene_id:
                    if gene_name2id:
                        gene_name_list.append(dict_trans[z])
                    pass
                if map_id == k:
                    insert_data6["gene_id"] = gene_id
                    insert_data6["gene_name"] = gene_name_list
                    insert_data6["gene_num"] = df_trans_keggp_table.loc[i, "count"]
            data_son3 = SON(insert_data6)
            data_list.append(data_son3)
        try:
            collection = self.db['relation_keggp_detail']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格relation_keggp_detail信息出错:%s", variables=(e), code="54700107")
        else:
            pass
        self.bind_object.logger.info("导入表格relation_keggp_detail信息成功!")


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

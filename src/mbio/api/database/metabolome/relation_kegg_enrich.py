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
import numpy as np
from mbio.packages.metabolome.common import check_metab


class RelationKeggEnrich(Base):
    def __init__(self, bind_object):
        super(RelationKeggEnrich, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_relation_keggp_enrich(self, name=None, main_id=None, params =None):
        if not main_id:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'RelationKeggEnrich',
                'params': params if params else '',
                'status': 'end',
                'main_id': ''
            }
            try:
                collection = self.db['relation_keggp_enrich']
                relation_kegg_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", relation_kegg_id, relation_kegg_id)
            except Exception, e:
                self.bind_object.set_error('导入relation_keggp_enrich主表异常:%s', variables=(e), code="54700101")
        else:
            self.update_table("main_id", main_id, main_id)
            relation_kegg_id = main_id
        return relation_kegg_id

    @report_check
    def add_relation_keggp_enrich_detail(self, relation_kegg_id, metab_enrich_result=None, trans_enrich_result=None, metab_desc=None, gene_name2id=None, species_name=None, version=None):
        if version in ['v94.2', '94.2']:
            ref_db = self.ref_db["kegg_v94.2_organisms"]
        elif version in ['kegg', '']:
            ref_db = self.ref_db["kegg_organisms"]
        else:
            ref_db = self.ref_db["kegg_v2021.09.18_organisms"]
        df_metab_desc = pd.read_table(metab_desc, '\t')
        df_metab_desc = df_metab_desc.fillna("-")
        dict_metab = df_metab_desc.set_index(["metab_id"])["Metabolite"].to_dict()
        if gene_name2id:
            df_gene_name2id = pd.read_table(gene_name2id, '\t')
            df_trans_desc = df_gene_name2id.fillna("-")
            dict_trans = df_trans_desc.set_index(["gene_id"])["gene_name"].to_dict()                
        metab_enrich = pd.read_table(metab_enrich_result, '\t')
        trans_enrich = pd.read_table(trans_enrich_result, '\t')
        metab_pathway_list = metab_enrich["ID"].tolist()
        trans_pathway_list = trans_enrich["ID"].tolist()
        trans_newpathway_list = []
        for pathway in trans_pathway_list:
            if not pathway.startswith("map"):
                species = pathway[0:-5]  # 获取转录的物种缩写
                rename = pathway.replace(species, "map")  # 将物种替换成map
                trans_newpathway_list.append(rename)
            else:
                trans_newpathway_list.append(pathway)
                species = "map"
        sp_species = species_name.split(";")
        self.bind_object.logger.info("sp_species为{}".format(sp_species))
        if len(sp_species) == 1:
            pass
        else:
            res = ref_db.find_one({"second_category": sp_species[1]})
            species = res["name"]
        self.bind_object.logger.info("species为{}".format(species))
        set1 = set(metab_pathway_list)
        self.bind_object.logger.info("set1为{}".format(set1))
        set2 = set(trans_newpathway_list)
        set_metab = set1-set2
        set_trans = set2-set1
        set_intersection = set1 & set2
        self.bind_object.logger.info("set_metab为{}".format(set_metab))
        self.bind_object.logger.info("set_metab为{}".format(len(set_metab)))
        self.bind_object.logger.info("set_trans为{}".format(set_trans))
        self.bind_object.logger.info("set_trans为{}".format(len(set_trans)))
        self.bind_object.logger.info("set_intersection为{}".format(set_intersection))
        self.bind_object.logger.info("set_intersection为{}".format(len(set_intersection)))
        data_list = []
        empty_list = []
        for i in metab_enrich.index:
            metab_enrich_factor = float(metab_enrich.loc[i, "Ratio_in_study"].split("/")[0])/float(metab_enrich.loc[i, "Ratio_in_pop"].split("/")[0])
            metab_list = metab_enrich.loc[i, "Metab_ids"].split('|')
            metab_name_list = []
            for a in metab_list:
                metab_name_list.append(dict_metab[a])
            if metab_enrich.loc[i, "ID"] in set_metab:
                pathway_id = metab_enrich.loc[i, "ID"]
                insert_data1 = {
                    "kegg_id":relation_kegg_id,
                    "pathway_id": pathway_id.replace(pathway_id[0:-5], species),
                    "description":metab_enrich.loc[i, "Term"],
                    "database": metab_enrich.loc[i, "Database"],
                    "first_category": metab_enrich.loc[i, "typeI"],
                    "second_category": metab_enrich.loc[i, "typeII"],
                    "metab_id": metab_list,
                    "metab_name": metab_name_list,
                    "metab_num": len(metab_list),
                    "metab_pvalue": metab_enrich.loc[i, "P-Value"],
                    "metab_qvalue": metab_enrich.loc[i, "Corrected P-Value"],
                    "metab_log_p": -np.log10(float(metab_enrich.loc[i, "P-Value"])),
                    "metab_log_q": -np.log10(float(metab_enrich.loc[i, "Corrected P-Value"])),
                    "metab_ratio_in_study": metab_enrich.loc[i, "Ratio_in_study"],
                    "metab_ratio_in_pop": metab_enrich.loc[i, "Ratio_in_pop"],
                    "metab_enrich_factor": metab_enrich_factor,
                    "gene_id": empty_list,
                    "gene_name": empty_list,
                    "gene_num": 0,
                    "gene_pvalue": -1,
                    "gene_qvalue": -1,
                    "gene_log_p": -1,
                    "gene_log_q": -1,
                    "gene_ratio_in_study": "-",
                    "gene_ratio_in_pop": "-",
                    "gene_enrich_factor": 0
                }
                data_son1 = SON(insert_data1)
                data_list.append(data_son1)
        for j in trans_enrich.index:
            pathway_id = trans_enrich.loc[j, "ID"]
            map_id = pathway_id.replace(pathway_id[0:-5], "map")
            trans_enrich_factor = float(trans_enrich.loc[j, "Ratio_in_study"].split("/")[0])/float(trans_enrich.loc[j, "Ratio_in_pop"].split("/")[0])
            trans_list = trans_enrich.loc[j, "Genes"].split('|')
            trans_name_list = []
            if gene_name2id:
                for b in trans_list:
                    trans_name_list.append(dict_trans[b])
            else:
                pass
            if map_id in set_trans:
                pathway_id = trans_enrich.loc[j, "ID"]
                insert_data2 = {
                    "kegg_id": relation_kegg_id,
                    "pathway_id": pathway_id.replace(pathway_id[0:-5], species),
                    "description":trans_enrich.loc[j, "Term"],
                    "database": trans_enrich.loc[j, "Database"],
                    "first_category": trans_enrich.loc[j, "typeI"],
                    "second_category": trans_enrich.loc[j, "typeII"],
                    "metab_name": empty_list,
                    "metab_id": empty_list,
                    "metab_num": 0,
                    "metab_pvalue": -1,
                    "metab_qvalue": -1,
                    "metab_log_p": -1,
                    "metab_log_q": -1,
                    "metab_ratio_in_study": '-',
                    "metab_ratio_in_pop": '-',
                    "metab_enrich_factor": 0,
                    "gene_id": trans_list,
                    "gene_name": trans_name_list,
                    "gene_num": len(trans_list),
                    "gene_pvalue": trans_enrich.loc[j, "P-Value"],
                    "gene_qvalue": trans_enrich.loc[j, "Corrected P-Value"],
                    "gene_log_p": -np.log10(float(trans_enrich.loc[j, "P-Value"])),
                    "gene_log_q": -np.log10(float(trans_enrich.loc[j, "Corrected P-Value"])),
                    "gene_ratio_in_study": trans_enrich.loc[j, "Ratio_in_study"],
                    "gene_ratio_in_pop": trans_enrich.loc[j, "Ratio_in_pop"],
                    "gene_enrich_factor": trans_enrich_factor
                }
                data_son2 = SON(insert_data2)
                data_list.append(data_son2)
        for k in set_intersection:
            for i in metab_enrich.index:
                pathway_id = metab_enrich.loc[i, "ID"]
                metab_enrich_factor = float(metab_enrich.loc[i, "Ratio_in_study"].split("/")[0])/float(metab_enrich.loc[i, "Ratio_in_pop"].split("/")[0])
                metab_list = metab_enrich.loc[i, "Metab_ids"].split('|')
                metab_name_list = []
                for a in metab_list:
                    metab_name_list.append(dict_metab[a])
                if metab_enrich.loc[i, "ID"] == k:
                    insert_data3 = {
                                "kegg_id": relation_kegg_id,
                                "pathway_id": pathway_id.replace(pathway_id[0:-5], species),
                                "description": metab_enrich.loc[i, "Term"],
                                "database": metab_enrich.loc[i, "Database"],
                                "first_category": metab_enrich.loc[i, "typeI"],
                                "second_category": metab_enrich.loc[i, "typeII"],
                                "metab_id": metab_list,
                                "metab_name": metab_name_list,
                                "metab_num": len(metab_list),
                                "metab_pvalue": metab_enrich.loc[i, "P-Value"],
                                "metab_qvalue": metab_enrich.loc[i, "Corrected P-Value"],
                                "metab_log_p": -np.log10(float(metab_enrich.loc[i, "P-Value"])),
                                "metab_log_q": -np.log10(float(metab_enrich.loc[i, "Corrected P-Value"])),
                                "metab_ratio_in_study": metab_enrich.loc[i, "Ratio_in_study"],
                                "metab_ratio_in_pop": metab_enrich.loc[i, "Ratio_in_pop"],
                                "metab_enrich_factor": metab_enrich_factor
                                }
            for j in trans_enrich.index:
                pathway_id = trans_enrich.loc[j, "ID"]
                map_id = pathway_id.replace(pathway_id[0:-5], "map")
                trans_enrich_factor = float(trans_enrich.loc[j, "Ratio_in_study"].split("/")[0])/float(trans_enrich.loc[j, "Ratio_in_pop"].split("/")[0])
                trans_list = trans_enrich.loc[j, "Genes"].split('|')
                trans_name_list = []
                if gene_name2id:
                    for b in trans_list:
                        trans_name_list.append(dict_trans[b])
                else:
                    pass
                if map_id == k:
                    insert_data3["gene_id"] = trans_list
                    insert_data3["gene_name"] = trans_name_list
                    insert_data3["gene_num"] = len(trans_list)
                    insert_data3["gene_pvalue"] = trans_enrich.loc[j, "P-Value"]
                    insert_data3["gene_qvalue"] = trans_enrich.loc[j, "Corrected P-Value"]
                    insert_data3["gene_log_p"] = -np.log10(float(trans_enrich.loc[j, "P-Value"]))
                    insert_data3["gene_log_q"] = -np.log10(float(trans_enrich.loc[j, "Corrected P-Value"]))
                    insert_data3["gene_ratio_in_study"] = trans_enrich.loc[j, "Ratio_in_study"]
                    insert_data3["gene_ratio_in_pop"] = trans_enrich.loc[j, "Ratio_in_pop"]
                    insert_data3["gene_enrich_factor"] = trans_enrich_factor
            data_son3 = SON(insert_data3)
            data_list.append(data_son3)
        try:
            collection = self.db['relation_keggp_enrich_detail']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格relation_kegg_enrich_detail信息出错:%s", variables=(e), code="54700107")
        else:
            pass
        self.bind_object.logger.info("导入表格relation_kegg_enrich_detail信息成功!")

    @report_check
    def update_table(self, str, name, main_table_id):
        try:
            self.db['relation_corr_network'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {str: name}})
        except Exception as e:
            self.bind_object.set_error('relation_corr_network%s字段出错:%s', variables=(str,e), code="54700109")


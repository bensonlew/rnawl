# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'
# last modify date: 2022.01.25
# last modified: zhaoyuzhuo

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import numpy as np
import commands
from mbio.packages.metabolome.get_data import dump_trans_data
from biocluster.config import Config
from bson.objectid import ObjectId


class RelationKeggpAgent(Agent):
    """
    代谢转录关联分析--筛选转录kegg表
    """

    def __init__(self, parent):
        super(RelationKeggpAgent, self).__init__(parent)
        options = [
            {"name": "metab_set_table_id", "type": "string"},
            {"name": "trans_keggp_main_id", "type": "string"},
            {"name": "trans_geneset_main_id", "type": "string"},
            {"name": "task_id", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(RelationKeggpAgent, self).end()


class RelationKeggpTool(Tool):
    def __init__(self, config):
        super(RelationKeggpTool, self).__init__(config)
        self.dele_med_gloabl = ['map01100', 'map01110', 'map01120', 'map01130', 'map01200', 'map01210', 'map01212', 'map01230', 'map01220', 'hsa01100',
                                 'hsa01110', 'hsa01120', 'hsa01130', 'hsa01200', 'hsa01210', 'hsa01212', 'hsa01230', 'hsa01220', 'mmu01100', 'mmu01110',
                                'mmu01120', 'mmu01130', 'mmu01200', 'mmu01210', 'mmu01212', 'mmu01230', 'mmu01220', 'rna01100', 'rna01110', 'rna01120',
                                 'rna01130', 'rna01200', 'rna01210', 'rna01212', 'rna01230', 'rna01220']

    def run(self):
        super(RelationKeggpTool, self).run()
        self.get_metab_keggp_table()
        self.get_trans_keggp_table()
        self.set_output()
        self.end()

    def get_metab_keggp_table(self):
        client = Config().get_mongo_client(mtype="metabolome")
        db = client[Config().get_mongo_dbname("metabolome")]
        metab_keggp_file = self.output_dir + "/metab_keggp_table.xls"
        head_list = ['pathway_id', 'description', 'first_category', 'second_category', 'compound_id', 'metab_id',
                     'count']
        f = open(metab_keggp_file, 'w+')
        f.write("\t".join(head_list) + "\n")
        search_result = db['metab_set'].find_one({"task_id": self.option('task_id'), "_id": ObjectId(self.option('metab_set_table_id'))})
        try:
            metab_set_name = search_result["name"]
            search_main_id_result = db['metabset_keggp'].find_one(
                {"task_id": self.option('task_id'), "metabset_list": metab_set_name})
            main_id = search_main_id_result["main_id"]
            result = db['metabset_keggp_level'].find({"kegg_id": main_id})
            for one in result:
                mongo_name_list = ['pathway_id', 'description', 'first_category', 'second_category', 'compound_id', 'metab_id', 'count']
                write_list = map(lambda x: str(one[x]), mongo_name_list)
                f.write("\t".join(write_list) + "\n")
        except:  # 如果没有注释结果，重新进行注释分析
            anno_overview_id = db["anno_overview"].find_one({"task_id": self.option('task_id')})
            metabset_main_id = search_result["main_id"]
            metabset_detail = db["metab_set_detail"].find_one({"set_id": metabset_main_id})
            metabset_list = metabset_detail["set_list"]
            pathway_list = []
            for i in metabset_list:
                anno_overview_res = db["anno_overview_detail"].find_one(
                    {"overview_id": anno_overview_id["main_id"], "metab_id": i})
                pathway_res = anno_overview_res["pathway_id"].split(";")
                pathway_list += pathway_res
            self.logger.info("pathway_list_len为{}".format(len(pathway_list)))
            pathway_list_treat = list(set(pathway_list))
            pathway_list_treat.remove("-")
            self.logger.info("pathway_list_treat_len为{}".format(len(pathway_list_treat)))
            for j in pathway_list_treat:
                self.logger.info("j为{}".format(j))
                try:
                    anno_overview_ko_res = db["anno_overview_ko"].find_one({"overview_id": anno_overview_id["main_id"], "pathway_id": j})
                    self.logger.info("anno_keggp_level_res为{}".format(anno_overview_ko_res))
                    mongo_name_list = ['pathway_id', 'description', 'first_category', 'second_category', 'compound_id', 'metab_id']
                    write_list = map(lambda x: str(anno_overview_ko_res[x]), mongo_name_list)
                    write_list.append(str(len(anno_overview_ko_res["metab_id"].split(";"))))
                    self.logger.info("write_list为{}".format(write_list))
                    f.write("\t".join(write_list) + "\n")
                except:
                    continue
        f.close()

    def get_trans_keggp_table(self):
        self.logger.info("start get_keggp_table")
        metab_client = Config().get_mongo_client(mtype="metabolome")
        relation_db = metab_client[Config().get_mongo_dbname("metabolome")]
        relation_info = relation_db['sg_relation_analysis'].find_one(
            {"task_id": self.option('task_id'), "delete_by": ""})
        relate_task_id = relation_info["relate_task_id"]
        relate_project_type = relation_info["relate_project_type"]
        try:
            db_version = relation_info["relate_db_version"]
        except:
            db_version = 1
        self.logger.info("relate_task_id为{}".format(relate_task_id))
        self.logger.info("relate_project_type为{}".format(relate_project_type))
        self.logger.info("db_version为{}".format(db_version))
        trans_client = Config().get_mongo_client(mtype=relate_project_type)
        trans_db = trans_client[Config().get_mongo_dbname(relate_project_type)]
        # 基因集表
        geneset_file = self.work_dir + "/geneset_table.xls"
        dump_trans_data(proj_type=relate_project_type, task_id=relate_task_id, col_type="geneset",
                        main_id=self.option("trans_geneset_main_id"), db_version=db_version, outfile=geneset_file)
        # kegg总表
        trans_keggp_file = self.work_dir + "/trans_keggp_total_table.xls"
        dump_trans_data(proj_type=relate_project_type, task_id=relate_task_id, col_type="kegg_anno",
                        main_id=self.option("trans_keggp_main_id"), db_version=db_version, outfile=trans_keggp_file)
        # 基因name和id对应表
        genename_table = self.work_dir + "/gene_id2name.xls"
        dump_trans_data(proj_type=relate_project_type, task_id=relate_task_id, col_type="gene_name",
                        db_version=db_version, outfile=genename_table)
        df_geneset = pd.read_table(geneset_file, '\t')
        geneset_list = df_geneset["gene_id"].tolist()
        df_kegg_table = pd.read_table(trans_keggp_file, '\t')
        df_kegg_table = df_kegg_table.dropna()
        index_list = []
        for i in df_kegg_table.index:
            if df_kegg_table.loc[i, "gene_id"] in geneset_list:
                index_list.append(i)
        self.logger.info('index_list={}'.format(len(index_list)))
        df_sort_table = df_kegg_table.ix[index_list, :]
        gene_id_list = []
        pathway_list = []
        for i in df_sort_table.index:
            pathway = df_sort_table.loc[i, "pathways"].split(';')
            for j in pathway:
                pathway_list.append(j)
                gene_id_list.append(df_sort_table.loc[i, "gene_id"])
        self.logger.info('pathway_set_len={}'.format(len(set(pathway_list))))
        self.logger.info('pathway_set={}'.format(set(pathway_list)))
        self.logger.info('gene_id_list={}'.format(len(gene_id_list)))
        dict_p2id = {}
        for x in range(len(pathway_list)):
            if pathway_list[x] in dict_p2id:
                dict_p2id[pathway_list[x]] += [gene_id_list[x]]
            else:
                dict_p2id[pathway_list[x]] = [gene_id_list[x]]
        trans_sort_keggp_table = self.output_dir + "/trans_keggp_table.xls"
        with open(trans_sort_keggp_table, 'w+') as outfile:
            outfile.write("pathway\tgene_id\tcount\tdescription\tfirst_category\tsecond_category\n")
            for pathway_id, y in dict_p2id.items():
                self.logger.info('pathwayid={}'.format(pathway_id))
                if relate_project_type == "whole_transcriptome":
                    kegg_anno = trans_db["annotation_kegg"].find_one({"task_id": relate_task_id})
                    res = trans_db["annotation_kegg_level"].find_one(
                        {"kegg_id": kegg_anno["main_id"], "pathway_id": pathway_id})
                    outfile.write(
                        pathway_id + "\t" + ";".join(y) + "\t" + str(len(y)) + "\t" + res["pathway_definition"] + "\t" +
                        res["first_category"] + "\t" + res["second_category"] + "\n")
                else:
                    kegg_anno = trans_db["sg_annotation_kegg"].find_one({"task_id": relate_task_id})
                    try:  # 在医学有参中一些项目会筛除一些pathwayid，获取不到注释结果的直接跳过
                        res = trans_db["sg_annotation_kegg_level"].find_one(
                            {"kegg_id": kegg_anno["main_id"], "pathway_id": pathway_id})
                        self.logger.info('res={}'.format(res))
                        outfile.write(pathway_id + "\t" + ";".join(y) + "\t" + str(len(y)) + "\t" + res[
                            "pathway_definition"] + "\t" + res["first_category"] + "\t" + res["second_category"] + "\n")
                    except:
                        continue

    def set_output(self):
        self.logger.info('开始设置输出结果文件')
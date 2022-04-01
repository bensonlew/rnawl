# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# last_modify:20180309

import os
import pandas as pd
import numpy as np
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import json
import re
import gridfs
import glob
from collections import defaultdict
from mbio.api.database.prok_rna.api_base import ApiBase
import math
import time


class Geneset(ApiBase):
    def __init__(self, bind_object):
        super(Geneset, self).__init__(bind_object)
        self._project_type = 'prok_rna'

    @report_check
    def add_main_table(self, collection_name, params, name):
        """
        添加主表的导表函数
        :param collection_name: 主表的collection名字
        :param params: 主表的参数
        :param name: 主表的名字
        :return:
        """
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "status": "end",
            "name": name,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "params": json.dumps(params, sort_keys=True, separators=(',', ':'))
        }

        collection = self.db[collection_name]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_geneset_cog_detail(self, geneset_cog_table, geneset_cog_id):
        """
        cog详情表导表函数
        :param geneset_cog_table:cog结果表
        :param geneset_cog_id:主表ID
        :return:
        """
        data_list = []
        geneset_name = []
        cogclass = defaultdict(list)
        name2list = defaultdict(list)
        with open(geneset_cog_table, 'r') as f:
            gene_list = list()
            first_line = f.readline().strip("\n").split("\t")
            for gn in first_line[2:]:
                if "list" in gn or "LIST" in gn:
                    continue
                elif not gn[:-4] in geneset_name:
                    geneset_name.append(gn[:-4])
            self.bind_object.logger.info(geneset_name)
            for line in f:
                line = line.strip().split("\t")
                data = {
                    'geneset_cog_id': ObjectId(geneset_cog_id),
                    'type': line[0],
                    'function_categories': line[1]
                }
                if line[0] not in cogclass:
                    cogclass[line[0]] = list()
                if line[1] not in cogclass[line[0]]:
                    cogclass[line[0]].append(line[1])
                for n, gn in enumerate(geneset_name):
                    data[gn + "_cog"] = int(line[2*n+2])
                    # data[gn + "_nog"] = int(line[6*n+3])
                    # data[gn + "_kog"] = int(line[6*n+4])
                    data[gn + "_cog_list"] = line[2*n+3].split(";")
                    if data[gn + "_cog_list"] == ["none"]:
                        data[gn + "_cog_list"] = []
                    # gene_list += line[2*n+3].split(";")
                    name2list[gn] += line[2*n+3].split(";")
                    # data[gn + "_nog_list"] = line[6*n+6].split(";")
                    # data[gn + "_kog_list"] = line[6*n+7].split(";")
                    # data[gn + "_cog_str"] = line[6*n+5]
                    # data[gn + "_nog_str"] = line[6*n+6]
                    # data[gn + "_kog_str"] = line[6*n+7]
                data_list.append(data)
        # self.insert_genelist_forquery(list(set(gene_list)), geneset_cog_id)
        try:
            collection = self.db['sg_geneset_cog_class_detail']
            main_collection = self.db['sg_geneset_cog_class']
            if len(data_list) != 0:
                collection.insert_many(data_list)
            main_collection.update({"_id": ObjectId(geneset_cog_id)},
                                   {"$set": {"table_columns": geneset_name, "status": "end", "cog_class": cogclass}})
            self.bind_object.logger.info(geneset_name)
        except Exception as e:
            self.bind_object.set_error("导入cog表格：%s出错:%s" % (geneset_cog_table, e))
        else:
            self.bind_object.logger.info("导入cog表格：%s成功!" % (geneset_cog_table))

# -----添加基因对应关系表------
        data_list = []
        all_gene = set(sum(name2list.values(), []))
        for gene in all_gene:
            genesets = set()
            gene_name = description = ''
            cog_id = list()
            cog_description = list()
            cog_type = list()
            type_des = list()
            category = list()
            geneset_buer = dict()
            for name in name2list:
                if gene in name2list[name]:
                    # genesets.add(name)
                    geneset_buer[name] = True
                else:
                    geneset_buer[name] = False
            query_info = self.get_gene_query_info(gene)
            if query_info:
                gene_name = query_info["gene_name"]
                description = query_info["description"]
                cog_id = [x.strip().split('(')[0] for x in query_info['cog'].split(';')]
                cog_description = [x.strip().split('(')[1].strip(')') for x in query_info["cog_description"].split(';')]
                cog_type = [x.strip().split('(')[1].split(':')[0] for x in query_info["cog"].split(';')]
                type_des = [x.strip().split('(')[1].split(':')[1].strip(')') for x in query_info["cog"].split(';')]
                category = query_info['cog_class1'].split(';')
            type_log = set()
            for n, i in enumerate(cog_id):
                if cog_type[n] in type_log:
                    dup = True
                else:
                    dup = False
                    type_log.add(cog_type[n])
                data = {
                    'geneset_cog_id': ObjectId(geneset_cog_id),
                    'gene_id': gene,
                    'gene_name': gene_name,
                    'cog_id': i,
                    'cog_description': cog_description[n],
                    'cog_type': cog_type[n],
                    'type_description': type_des[n],
                    'category': category[n],
                    # 'genesets': ';'.join(genesets),
                    'description': description,
                    'dup': dup
                }
                data.update(geneset_buer)
                data_list.append(data)
        try:
            collection = self.db['sg_geneset_cog_class_gene']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入cog基因对应关系出错:" )
        else:
            self.bind_object.logger.info("导入cog基因对应关系成功!" )

        #-----在geneset_cog_class主表中添加注释表中的cog_class字段------

    @report_check
    def add_go_enrich_detail(self, go_enrich_id, go_enrich_dir):
        """
        GO富集详情导表函数
        :param go_enrich_id: 主表ID
        :param go_enrich_dir: 结果文件（不是文件夹）
        :return:
        """
        if not isinstance(go_enrich_id, ObjectId):
            if isinstance(go_enrich_id, types.StringTypes):
                go_enrich_id = ObjectId(go_enrich_id)
            else:
                raise Exception('go_enrich_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(go_enrich_dir):
            raise Exception('{}所指定的路径不存在。请检查！'.format(go_enrich_dir))
        data_list = []
        go_type = []
        gene_list = list()
        enrichgo_gene = dict()
        with open(go_enrich_dir, 'r') as f:
            f.readline()
            for line in f:
                line = line.strip().split('\t')
                if line[2] == 'p':
                    continue
                gene_list += line[-2].split(";")
                enrichgo_gene[line[0]] = line[-2]
                data = [
                    ('go_enrich_id', go_enrich_id),
                    ('go_id', line[0]),
                    ('go_type', line[1]),
                    ('enrichment', line[2]),
                    ('description', line[3]),
                    ('ratio_in_study', line[4]),
                    ('ratio_in_pop', line[5]),
                    ('p_uncorrected', float(line[6])),
                    ('p_corrected', float(line[-1])),
                    ('enrich_factor', float(line[4].split("/")[0])/float(line[5].split("/")[0])),
                    ('depth', int(line[7])),
                    ('study_count', int(line[4].split("/")[0])),
                    ('pop_count', int(line[5].split("/")[0])),
                    ('seq_list', line[-2]),
                    ('seq_num', len(line[-2].split(";")))
                ]
                go_type.append(line[1])
                data = SON(data)
                data_list.append(data)
            self.insert_genelist_forquery(gene_list, go_enrich_id)
        if data_list:
            # 插入-logpvalue -logpcorrected 值相关字段
            pvalues = [dict(son)['p_uncorrected'] for son in data_list]
            if len([x for x in pvalues if x > 0]) > 0:
                pvalues_min = min([x for x in pvalues if x > 0])/10
            else:
                pvalues_min = 0.0001
            pvalues_min = - math.log10(pvalues_min)
            log10x = [-math.log10(x) if x>0 else pvalues_min for x in pvalues]
            for i in range(0,len(log10x)):
                data_list[i]['neg_log10p_uncorrected'] = log10x[i]

            pvalues = [dict(son)['p_corrected'] for son in data_list]
            if len([x for x in pvalues if x > 0]) > 0:
                pvalues_min = min([x for x in pvalues if x > 0])/10
            else:
                pvalues_min = 0.0001
            pvalues_min = - math.log10(pvalues_min)
            log10x = [-math.log10(x) if x>0 else min(pvalues)/10 for x in pvalues]
            for i in range(0,len(log10x)):
                data_list[i]['neg_log10p_corrected'] = log10x[i]

            try:
                collection = self.db['sg_geneset_go_enrich_detail']
                collection.insert_many(data_list)
                coll = self.db['sg_geneset_go_enrich']
                go_type = list(set(go_type))
                coll.update({'_id': go_enrich_id}, {'$set': {'categories': go_type, 'status': 'end'}})
                # main_collection = self.db['sg_geneset_go_enrich']
                # main_collection.update({"_id": ObjectId(go_enrich_id)}, {"$set": {"status": "end"}})
            except Exception as e:
                print("导入go富集信息：%s出错:%s" % (go_enrich_dir, e))
            else:
                print("导入go富集信息：%s成功!" % (go_enrich_dir))
        else:
            raise Exception('GO富集没有结果')

        data_list = []
        for gene in set(gene_list):
            gene_name = description = go_level2 = ''
            go_description = list()
            go_id = list()
            query_info = self.get_gene_query_info(gene)
            if query_info:
                try:
                    gene_name = query_info["gene_name"]
                    description = query_info["description"]
                    go_level2 = query_info['go_lev2']
                    go_id = [x.split("(")[0].strip() for x in query_info['go'].split(';')]
                    print(go_id)
                    go_description = [x.split("(")[1].strip(')').strip() for x in query_info['go'].split(';')]
                    enrichgo = list()
                    for go in enrichgo_gene:
                        if gene in enrichgo_gene[go].split(';'):
                            enrichgo.append(go)
                except:
                    pass
            data = {
                'go_enrich_id': go_enrich_id,
                'gene_id': gene,
                'gene_name': gene_name,
                'description': description,
                'go_level2': go_level2,
                'go_id': ';'.join(go_id),
                'go_description': ';'.join(go_description),
                'enrichgo': ';'.join(list(set(enrichgo))),
            }
            data_list.append(data)
        try:
            collection = self.db['sg_geneset_go_enrich_gene']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入go富集基因对应关系出错:")
        else:
            self.bind_object.logger.info("导入go富集基因对应关系成功!")

    @report_check
    def update_directed_graph(self, go_enrich_id, go_graph_png, go_graph_pdf):
        collection = self.db['sg_geneset_go_enrich']
        fs = gridfs.GridFS(self.db)
        gra = fs.put(open(go_graph_png, 'rb'))
        gra_pdf = fs.put(open(go_graph_pdf, 'rb'))
        try:
            collection.update({"_id": ObjectId(go_enrich_id)},
                              {"$set": {'go_directed_graph': gra, "graph_pdf": gra_pdf}})
        except Exception as e:
            print("导入%s信息出错：%s" % (go_graph_png, e))
        else:
            print("导入%s信息成功！" % (go_graph_png))

    @report_check
    def add_kegg_enrich_detail(self, enrich_id, kegg_enrich_table):
        """
        KEGG富集详情表导表函数
        :param enrich_id: 主表id
        :param kegg_enrich_table: 结果表
        :return:
        """
        if not isinstance(enrich_id, ObjectId):
            if isinstance(enrich_id, types.StringTypes):
                enrich_id = ObjectId(enrich_id)
            else:
                raise Exception('kegg_enrich_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(kegg_enrich_table):
            raise Exception('kegg_enrich_table所指定的路径:{}不存在，请检查！'.format(kegg_enrich_table))
        data_list = []
        # geneset_length = len(open(geneset_list_path, "r").readlines())
        # all_list_length = len(open(all_list_path, "r").readlines())
        kegg_type1 = []
        gene_list = list()
        with open(kegg_enrich_table, 'rb') as r:
            for line in r:
                if re.match(r'\w', line):
                    line = line.strip('\n').split('\t')
                    gene_list += line[8].split('|')
                    insert_data = {
                        'kegg_enrich_id': enrich_id,
                        'term': line[1],
                        'database': line[2],
                        'id': line[3].split("path:")[1] if "path:" in line[3] else line[3],
                        'discription': line[1],
                        'study_count': int(line[0]),
                        "background_number": line[5].split("/")[1],
                        'ratio_in_study': line[4],
                        'ratio_in_pop': line[5],
                        'enrich_factor': float(line[0])/float(line[5].split("/")[0]),
                        'pvalue': round(float(line[6]), 4),
                        'corrected_pvalue': round(float(line[7]), 4) if not line[7] == "None" else "None",
                        'seq_list': line[8],
                        'seq_num': len(line[8].split('|')),
                        'hyperlink': line[9],
                        'kegg_type': "".join([x[0] for x in line[11].split(' ')])
                    }
                    kegg_type1.append("".join([x[0] for x in line[11].split(' ')]))
                    data_list.append(insert_data)
            # self.insert_genelist_forquery(list(set(gene_list)), enrich_id)
            if data_list:
                # 插入-logpvalue -logpcorrected 值相关字段
                pvalues = [dict(son)['pvalue'] for son in data_list]
                if len([x for x in pvalues if x > 0]) > 0:
                    pvalues_min = min([x for x in pvalues if x > 0])/10
                else:
                    pvalues_min = 0.0001
                pvalues_min = - math.log10(pvalues_min)
                log10x = [-math.log10(x) if x>0 else pvalues_min for x in pvalues]
                for i in range(0,len(log10x)):
                    data_list[i]['neg_log10p_uncorrected'] = log10x[i]

                pvalues = [dict(son)['corrected_pvalue'] for son in data_list]
                if len([x for x in pvalues if x > 0]) > 0:
                    pvalues_min = min([x for x in pvalues if x > 0])/10
                else:
                    pvalues_min = 0.0001
                pvalues_min = - math.log10(pvalues_min)

                log10x = [-math.log10(x) if x>0 else pvalues_min for x in pvalues]
                for i in range(0,len(log10x)):
                    data_list[i]['neg_log10p_corrected'] = log10x[i]

                kegg_type1 = list(set(kegg_type1))
                try:
                    collection = self.db['sg_geneset_kegg_enrich_detail']
                    collection.insert_many(data_list)
                    coll = self.db['sg_geneset_kegg_enrich']

                    coll.update({'_id': enrich_id}, {'$set': {'categories': kegg_type1, "status": "end"}})
                    # main_collection = self.db['sg_geneset_kegg_enrich']
                    # main_collection.update({"_id": ObjectId(enrich_id)}, {"$set": {"status": "end"}})
                except Exception as e:
                    self.bind_object.set_error("导入kegg富集统计表：%s信息出错:%s" % (kegg_enrich_table, e))
                else:
                    self.bind_object.logger.info("导入kegg富集统计表:%s信息成功!" % kegg_enrich_table)
            else:
                coll = self.db['sg_geneset_kegg_enrich']
                coll.update({'_id': enrich_id}, {'$set': {'desc': 'no_result'}})
                # self.bind_object.logger.info("kegg富集统计表没结果：" % kegg_enrich_table)
                raise Exception("kegg富集统计表没结果")

        data_list = []
        for gene in set(gene_list):
            gene_name = description = ko_id = ko_description = ''
            pathway_id = list()
            pathway_description = list()
            query_info = self.get_gene_query_info(gene)
            if query_info:
                gene_name = query_info["gene_name"]
                description = query_info["description"]
                ko_id = query_info['ko_id']
                ko_name = query_info['ko_name']
                # ko_description = query_info['ko_des']
                pathway_id = [x.split('(')[0].strip() for x in query_info['pathways'].split(';')]
                pathway_description = [x.split('(')[1].strip(')') for x in query_info['pathways'].split(';')]
            data = {
                'kegg_enrich_id': enrich_id,
                'gene_id': gene,
                'gene_name': gene_name,
                'description': description,
                'ko_id': ko_id,
                'ko_name': ko_name,
                'pathway_description': ';'.join(pathway_description),
                'pathway_id': ';'.join(pathway_id),
            }
            data_list.append(data)
        try:
            collection = self.db['sg_geneset_kegg_enrich_gene']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入kegg基因对应关系出错:")
        else:
            self.bind_object.logger.info("导入kegg基因对应关系成功!")

    @report_check
    def add_go_regulate_detail(self, go_regulate_dir, go_regulate_id, gene_go=None):
        """
        :param go_regulate_id: 主表ID
        :param go_regulate_dir: GO上下调结果
        :return:
        """
        data_list = []
        if not isinstance(go_regulate_id, ObjectId):
            if isinstance(go_regulate_id, types.StringTypes):
                go_regulate_id = ObjectId(go_regulate_id)
            else:
                raise Exception('go_enrich_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(go_regulate_dir):
            raise Exception('{}所指定的路径不存在，请检查！'.format(go_regulate_dir))
        with open(go_regulate_dir, 'r') as f:
            # gene_list = list()
            first_line = f.readline().strip().split("\t")
            doc_keys = []
            for l in first_line[3:]:
                name = l.split(" ")[0]
                if name not in doc_keys:
                    doc_keys.append(name)
            geneset_name = set(doc_keys)
            name2list = defaultdict(list)
            for line in f:
                line = line.strip().split('\t')
                data = {
                    'go_regulate_id': go_regulate_id,
                    'go_type': line[0],
                    'go': line[1],
                    'go_id': line[2]
                }
                for n, dk in enumerate(doc_keys):
                    line4 = line[4+n*3].split("(")
                    data["{}_num".format(dk)] = int(line[3+n*3])
                    data["{}_percent".format(dk)] = float(line4[0])
                    try:
                        gene_l = line[5+n*3]
                        if gene_l == "none":
                            gene_l = ''
                        # data["{}_str".format(dk)] = gene_l
                        # data["{}_genes".format(dk)] = gene_l.split(";")
                        # gene_list += gene_l.split(";")
                        name2list[dk] += gene_l.split(";")
                    except:
                        data["{}_str".format(dk)] = ""
                        # data["{}_genes".format(dk)] = ""
                    if len(line4) > 1:
                        data["{}_percent_str".format(dk)] = line4[1][:-1]
                    else:
                        data["{}_percent_str".format(dk)] = 0
                data_list.append(data)
        all_gene = sum(name2list.values(), [])
        while '' in all_gene:
            all_gene.remove('')
        # self.insert_genelist_forquery(list(set(gene_list)), go_regulate_id)
        try:
            collection = self.db['sg_geneset_go_class_detail']
            main_collection = self.db['sg_geneset_go_class']
            if len(data_list) != 0:
                collection.insert_many(data_list)
            main_collection.update({"_id": ObjectId(go_regulate_id)}, {"$set": {"table_columns": list(geneset_name), "status": "end"}})
            self.bind_object.logger.info(geneset_name)
            self.bind_object.logger.info(ObjectId(go_regulate_id))
        except Exception as e:
            self.bind_object.logger.info("导入go调控信息：%s出错:%s" % (go_regulate_dir, e))
        else:
            self.bind_object.logger.info("导入go调控信息：%s成功!" % (go_regulate_dir))

    #--------导入基因注释对应关系表---------
        data_list = []
        if gene_go:
            with open(gene_go, 'r') as f:
                f.readline()
                for line in f:
                    cols = line.strip().split("\t")
                    query_info = self.get_gene_query_info(cols[0])
                    if query_info:
                        gene_name = query_info["gene_name"]
                        description = query_info["description"]
                    else:
                        gene_name = ""
                        description = ""

                    geneset_buer = dict()
                    for name in name2list:
                        if cols[0] in name2list[name]:
                            # genesets.add(name)
                            geneset_buer[name] = True
                        else:
                            geneset_buer[name] = False
                    data = {
                        'go_regulate_id': ObjectId(go_regulate_id),
                        'gene_id': cols[0],
                        'gene_name': gene_name,
                        'description': description,
                        'go_level2': cols[1],
                        'gos_name': cols[2]
                    }
                    data.update(geneset_buer)
                    data_list.append(data)

        else:
            for gene in set(all_gene):
                genesets = set()
                gene_name = description = go_level2 = ''
                geneset_buer = dict()
                for name in name2list:
                    if gene in name2list[name]:
                        # genesets.add(name)
                        geneset_buer[name] = True
                    else:
                        geneset_buer[name] = False
                query_info = self.get_gene_query_info(gene)
                if query_info:
                    gene_name = query_info["gene_name"]
                    description = query_info["description"]
                    go_level2 = query_info['go_lev2']
                data = {
                    'go_regulate_id': ObjectId(go_regulate_id),
                    'gene_id': gene,
                    'gene_name': gene_name,
                    'description': description,
                    'go_level2': go_level2,
                    # 'genesets': ';'.join(genesets),
                }
                data.update(geneset_buer)
                data_list.append(data)
        try:
            collection = self.db['sg_geneset_go_class_gene']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入go基因对应关系出错:")
        else:
            self.bind_object.logger.info("导入go基因对应关系成功!")

    @report_check
    def add_go_regulate_detail2(self, go_regulate_dir, go_regulate_id):
        """
        :param go_regulate_id: 主表ID
        :param go_regulate_dir: GO上下调结果
        :return:
        """
        data_list = []
        if not isinstance(go_regulate_id, ObjectId):
            if isinstance(go_regulate_id, types.StringTypes):
                go_regulate_id = ObjectId(go_regulate_id)
            else:
                raise Exception('go_enrich_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(go_regulate_dir):
            raise Exception('{}所指定的路径不存在，请检查！'.format(go_regulate_dir))
        with open(go_regulate_dir, 'r') as f:
            first_line = f.readline().strip().split("\t")
            doc_keys = []
            for l in first_line[3:]:
                name = l.split(" ")[0]
                if name not in doc_keys:
                    doc_keys.append(name)
            geneset_name = set(doc_keys)
            for line in f:
                line = line.strip().split('\t')
                data = {
                    'go_regulate_id': go_regulate_id,
                    'go_type': line[0],
                    'go': line[1],
                    'go_id': line[2]
                }
                for n, dk in enumerate(doc_keys):
                    line4 = line[4+n*3].split("(")
                    data["{}_num".format(dk)] = int(line[3+n*3])
                    data["{}_percent".format(dk)] = float(line4[0])
                    try:
                        data["{}_str".format(dk)] = line[5+n*3]
                        # data["{}_genes".format(dk)] = line[5+n*3].split(";")
                    except:
                        data["{}_str".format(dk)] = ""
                        # data["{}_genes".format(dk)] = ""
                    if len(line4) > 1:
                        data["{}_percent_str".format(dk)] = line4[1][:-1]
                    else:
                        data["{}_percent_str".format(dk)] = 0
                data_list.append(data)
        try:
            time.sleep(30)
            collection = self.db['sg_geneset_go_class2_detail']
            main_collection = self.db['sg_geneset_go_class2']
            collection.insert_many(data_list)
            main_collection.update({"_id": ObjectId(go_regulate_id)}, {"$set": {"table_columns": list(geneset_name)}})
            self.bind_object.logger.info(geneset_name)
            self.bind_object.logger.info(ObjectId(go_regulate_id))
        except Exception as e:
            self.bind_object.logger.info("导入go调控信息：%s出错:%s" % (go_regulate_dir, e))
        else:
            self.bind_object.logger.info("导入go调控信息：%s成功!" % (go_regulate_dir))

    @report_check
    def add_kegg_regulate_pathway(self, pathway_dir, regulate_id):
        """

        :param regulate_id: 主表id
        :param pathway_dir:~/output/pathway 结果图片文件夹
        :return:
        """
        if not isinstance(regulate_id, ObjectId):
            if isinstance(regulate_id, types.StringTypes):
                regulate_id = ObjectId(regulate_id)
            else:
                raise Exception('kegg_regulate_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(pathway_dir):
            raise Exception('pathway_dir所指定的路径:{}不存在，请检查！'.format(pathway_dir))
        data_list = []
        png_files = glob.glob("{}/*.png".format(pathway_dir))
        pdf_files = glob.glob("{}/*.pdf".format(pathway_dir))
        fs = gridfs.GridFS(self.db)
        for f in png_files:
            # png_id = fs.put(open(os.path.join(pathway_dir, f), 'rb'))
            f_name = os.path.basename(f).split(".")[0]
            png_id = fs.put(open(f, 'rb'))
            pdf_id = fs.put(open(os.path.join(pathway_dir, f_name + ".pdf"), 'rb'))
            insert_data = {
                'kegg_id': regulate_id,
                'pathway_png': png_id,
                'pathway_pdf': pdf_id,
                'pathway_id': f_name
            }
            data_list.append(insert_data)
        try:
            collection = self.db['sg_geneset_kegg_class_pathway']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.info("导入kegg调控pathway：%s信息出错:%s" % (pathway_dir, e))
        else:
            self.bind_object.logger.info("导入kegg调控pathway:%s信息成功!" % pathway_dir)

    # 这个函数来自api_base.py，用在创建集合，也可以通过以前的传统方式

    def add_kegg_regulate_pic(self, main_table_id, level_path, png_dir, source=None):
        # 导入图片信息数据
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                kegg_id = ObjectId(main_table_id)
            else:
                raise Exception('main_table_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(level_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(level_path))
        if not os.path.exists(png_dir):
            raise Exception('{}所指定的路径不存在，请检查！'.format(png_dir))
        data_list = []
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                pid = re.sub('path:', '', line[0])
                if os.path.exists(png_dir + '/' + line[0] + '.html.mark'):
                    with open(png_dir + '/' + line[0] + '.html.mark', 'r') as mark_f:
                        for line_mark in mark_f.readlines():
                            # print len(line_mark.strip("\n").split("\t"))
                            if len(line_mark.strip("\n").split("\t")) >= 10:
                                [png, shape, bg_color, fg_color, coords, title, kos, href, gene_list, geneset_list] = line_mark.strip("\n").split("\t")
                                title = title.replace("\\n", "\n")
                            else:
                                continue

                            if bg_color == "":
                                continue

                            insert_data = {
                                'kegg_id': kegg_id,
                                'pathway_id': line[0],
                                'shape': shape,
                                'bg_colors': bg_color,
                                'fg_colors': fg_color,
                                'coords': coords,
                                'href': href,
                                'kos': kos,
                                'title': title
                            }

                            if source == "diff_exp":
                                geneset_reg_list = geneset_list.split("|")
                                geneset_list = ["_".join(x.split("_")[:-1]) for x in geneset_reg_list]
                                reg_list = [x.split("_")[-1] for x in geneset_reg_list]
                                insert_data.update({
                                    'gene_list': gene_list.split("|"),
                                    'geneset_list': geneset_list,
                                    'reg_list': reg_list
                                })
                            else:
                                insert_data.update({
                                    'gene_list': gene_list.split("|"),
                                    'geneset_list': geneset_list.split("|")
                                })

                        
                            if bg_color != "" and len(bg_color.split(",")) > 0:
                                insert_data.update({'bg_type': len(bg_color.split(","))})
                            if fg_color != "" and len(fg_color.split(",")) > 0:
                                insert_data.update({'fg_type': len(fg_color.split(","))})
                            data_list.append(insert_data)
                else:
                    self.bind_object.logger.info("kegg 图片{} 不存在html标记!".format(line[0]))

        if data_list:
            try:
                collection = self.db['sg_geneset_kegg_class_pic']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.set_error("导入kegg注释图片信息：%s、%s出错!" % (level_path, png_dir))
            else:
                self.bind_object.logger.info("导入kegg注释图片信息：%s、%s 成功!" % (level_path, png_dir))
    
    def create_db_table(self, table_name, content_dict_list, tag_dict=None):
        """
        Create main/detail table in database system.
        :param table_name: table name
        :param content_dict_list: list with dict as elements
        :param tag_dict: a dict to be added into each record in content_dict_list.
        :return: None or main table id
        """
        table_id = None
        conn = self.db[table_name]
        if tag_dict:
            for row_dict in content_dict_list:
                row_dict.update(tag_dict)
        record_num = len(content_dict_list)
        try:
            # split data and dump them to db separately
            if record_num > 5000:
                for i in range(0, record_num, 3000):
                    tmp_list = content_dict_list[i: i+3000]
                    conn.insert_many(tmp_list)
            else:
                if record_num >= 2:
                    conn.insert_many(content_dict_list)
                else:
                    table_id = conn.insert_one(content_dict_list[0]).inserted_id
        except Exception as e:
            if record_num >= 2:
                self.bind_object.set_error("Failed to insert records into table {} as: {}".format(table_name, e))
            else:
                self.bind_object.set_error("Failed to insert record into table {} as: {}".format(table_name, e))
        else:
            if record_num >= 2:
                self.bind_object.logger.info("Success to insert records into table {}".format(table_name))
            else:
                self.bind_object.logger.info("Success to insert record into table {}".format(table_name))

            return table_id

    def update_db_record(self, table_name, record_id, **kwargs):
        """
        根据记录的唯一标识record_id找到table_name的相关记录，并用kwargs更新记录
        :param table_name:
        :param record_id:
        :param kwargs: kwargs to add
        :return:
        """
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
           record_id = record_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        conn = self.db[table_name]
        conn.update({"_id": record_id}, {"$set": kwargs}, upsert=True)

    @report_check
    # 最后通过更新插入kegg_class主表的geneset的名字;gene_kegg_level_table_xls是交互workflowkegg_table_2对应的值
    # kegg_stat_xls是kegg_class这个tool产生的，也是通过这个文件来更新sg_geneset_kegg_class这个主表的table_columns字段
    def add_ipath_detail(self, main_table_id, ipath_input, geneset):
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(ipath_input):
            raise Exception('{}所指定的路径不存在，请检查！'.format(ipath_input))

        ipath = pd.read_table(ipath_input, header=0)
        ipath.columns = ['seq_id', 'ko', 'color', 'width']
        ipath['ipath_id'] = main_table_id
        row_dict_list = ipath.to_dict('records')
        main_collection = self.db['sg_geneset_ipath']

        doc_keys = []
        with open(geneset, 'r') as f:
            for l in f.readlines():
                name = l.strip().split('\t')[0]
                if name not in doc_keys:
                    doc_keys.append(name)
        geneset_name = list(set(doc_keys))


        try:
            self.create_db_table('sg_geneset_ipath_detail', row_dict_list)
            self.bind_object.logger.info("主表id：{} 基因集：{}".format(main_table_id, geneset_name))
            main_collection.update({"_id": main_table_id},
                                   {"$set": {"table_columns": geneset_name, "status": "end"}})
        except Exception as e:
            self.bind_object.set_error("导入ipath：%s出错!" % (ipath_input))
        else:
            self.bind_object.logger.info("导入ipath：%s出错!" % (ipath_input))


    def add_circ_detail(self, main_table_id, circ_input, enrich_type):
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                raise Exception('main_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(circ_input):
            raise Exception('{}所指定的路径不存在，请检查！'.format(circ_input))

        circ = pd.read_table(circ_input, header=0)
        circ.columns = ['seq_id', 'id', 'term', 'log2fc', 'gene_name']
        circ = circ.fillna('-')
        circ['circ_id'] = main_table_id
        row_dict_list = circ.to_dict('records')
        main_collection = self.db['sg_geneset_circ']

        try:
            self.create_db_table('sg_geneset_circ_detail', row_dict_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" % (circ_input))
        else:
            self.bind_object.logger.info("导入circ_detail：%s出错!" % (circ_input))

    def add_circ_graph(self, main_table_id, circ_input, enrich_type):
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(circ_input):
            raise Exception('{}所指定的路径不存在，请检查！'.format(circ_input))

        circ = pd.read_table(circ_input, header=0)
        circ = circ.fillna('-')
        columns=list(circ.columns)
        columns[0]='seq_id'
        columns[-2]='log2fc'
        columns[-1]='gene_name'

        circ.columns = columns
        circ['circ_id'] = main_table_id
        row_dict_list = circ.to_dict('records')
        main_collection = self.db['sg_geneset_circ']

        self.bind_object.logger.info("导入circ：%s出错!" % (row_dict_list))
        try:
            self.create_db_table('sg_geneset_circ_graph', row_dict_list)
        except Exception as e:
            self.bind_object.set_error("导入circ_graph：%s出错!" % (circ_input))
        else:
            self.bind_object.logger.info("导入circ：%s出错!" % (circ_input))

    def update_circ_main(self, main_table_id, circ_zscore_input, enrich_type):
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(circ_zscore_input):
            raise Exception('{}所指定的路径不存在，请检查！'.format(circ_zscore_input))

        circ_zscore = pd.read_table(circ_zscore_input, header=None)
        term_ids = list(circ_zscore[0])
        term_des = list(circ_zscore[1])
        term_zscores = list(circ_zscore[2])
        term  = [{term_ids[x]:[term_des[x], term_zscores[x]]} for x in range(0,len(term_ids))]
        main_collection = self.db['sg_geneset_circ']
        try:
            main_collection.update({"_id": main_table_id},
                                   {"$set": {"terms": term , "status": "end"}})
        except Exception as e:

            self.bind_object.set_error("更新circ主表：%s出错!" % (circ_zscore_input))
        else:
            self.bind_object.logger.info("导入ipath：%s出错!" % (main_table_id))
    def add_kegg_regulate_new2(self, main_table_id, geneset_file, kegg_stat_file, gene_kegg_level_table_xls, source=None):
        # kegg 分类导表
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                raise Exception('kegg_main_table_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(gene_kegg_level_table_xls):
            raise Exception('gene_kegg_level_table_xls所指定的路径:{}不存在，请检查！'.format(gene_kegg_level_table_xls))

        # 读入基因集列表
        with open(geneset_file, 'r') as pset:
            pset_list = [line.split("\t")[0] for line in pset.readlines()]
        stat = pd.read_table(kegg_stat_file, header=0)
        stat.columns=stat.columns.str.replace('Unnamed.*','link')
        level = pd.read_table(gene_kegg_level_table_xls, header=0)
        stat_class = pd.merge(stat, level, on='Pathway_id')

        # 按照kegg官网进行一级分类的排序
        list_custom = ['Metabolism',
                       'Genetic Information Processing',
                       'Environmental Information Processing',
                       'Cellular Processes',
                       'Organismal Systems',
                       'Human Diseases',
                       'Drug Development']
        first_class_index = dict(zip(list_custom,range(len(list_custom))))
        stat_class['first_rank']= stat_class['first_category'].map(first_class_index)
        stat_class.sort_values(['first_rank', 'second_category'], ascending = [True, True], inplace = True)

        stat_class.drop(['graph_id', 'hyperlink', 'graph_png_id', 'first_rank'], axis=1, inplace=True)
        stat_class.rename(columns={'Ko_ids':'ko_ids','Pathway_id':'pathway_id'}, inplace=True)
        stat_class.replace(np.nan, '', regex=True, inplace=True)
        stat_class['kegg_id'] = main_table_id

        # gene_list = list()
        name2list = defaultdict(list)
        for gene_set in pset_list:
            stat_class.rename(columns={ gene_set + '_genes' : gene_set + '_geneko'}, inplace=True)
            stat_class[gene_set + '_str']=stat_class[gene_set + '_geneko'].replace(r'\([^\)]*\)','',regex=True)
            # stat_class[gene_set + '_genes']=stat_class[gene_set + '_str'].map(lambda x: x.split(";"))
            stat_class[gene_set + '_genes'] = stat_class[gene_set + '_str'].map(lambda x: x.split(";") if x else [])  # modified by zhangyitong on 20211012
            stat_class[gene_set + '_str']=stat_class[gene_set + '_genes'].map(lambda x: ",".join(x))
            # gene_list += stat_class[gene_set + '_genes'].tolist()
            name2list[gene_set] += sum(stat_class[gene_set + '_genes'].tolist(), [])
        # gene_list = sum(gene_list, [])
        # self.insert_genelist_forquery(list(set(gene_list)), main_table_id)

        kegg_class_detail = stat_class.to_dict('records')
        self.create_db_table('sg_geneset_kegg_class_detail', kegg_class_detail)
        self.update_db_record('sg_geneset_kegg_class', ObjectId(main_table_id), main_id=ObjectId(main_table_id))

        # 导入统计数据
        all_second_cate = list(stat_class['second_category'])
        kegg_class_list = [j for i,j in enumerate(all_second_cate) if all_second_cate.index(j) == i]
        kegg_2to1 = dict(zip(stat_class['second_category'],stat_class['first_category']))

        data_list = list()
        for kegg_class in kegg_class_list:
            data = [
                ('first_category', kegg_2to1[kegg_class]),
                ('second_category', kegg_class),
                ('kegg_id', main_table_id)
            ]
            for gene_set in pset_list:
                class_genes = []
                genes_list = list(stat_class[stat_class['second_category'] == kegg_class][gene_set + '_genes'])
                for genes in genes_list:
                    class_genes.extend(genes)
                while '' in class_genes:
                    class_genes.remove('')
                class_genes = list(set(class_genes))
                data.extend([
                    (gene_set + '_genes', class_genes),
                    (gene_set + '_genes_num', len(class_genes))
                ])
            data = SON(data)
            data_list.append(data)
        categories = ["".join(map(lambda y:y[0], x.split(' '))) for x in list(set(stat_class['first_category']))]
        try:
            collection = self.db['sg_geneset_kegg_class_statistic']
            if len(data_list) > 0:
                collection.insert_many(data_list)
            self.update_db_record('sg_geneset_kegg_class', main_table_id, categories=categories, table_columns=pset_list, status="end")
        except Exception as e:
            raise Exception("导入kegg注释分类信息：%s出错 %s!" % (kegg_stat_file, e))
        else:
            self.bind_object.logger.info("导入kegg注释分类信息：%s 成功!" % kegg_stat_file)

        # --------导入基因注释对应关系表---------
        data_list = []
        all_gene = sum(name2list.values(), [])
        while '' in all_gene:
            all_gene.remove('')
        for gene in set(all_gene):
            genesets = set()
            gene_name = description = ko_id = ko_description = ''
            pathway_id = list()
            geneset_buer = dict()
            for name in name2list:
                if gene in name2list[name]:
                    # genesets.add(name)
                    geneset_buer[name] = True
                else:
                    geneset_buer[name] = False
            query_info = self.get_gene_query_info(gene)
            if query_info:
                gene_name = query_info["gene_name"]
                description = query_info["description"]
                pathways_class1 = query_info["pathways_class1"]
                pathways_class2 = query_info["pathways_class2"]
                ko_id = query_info['ko_id']
                ko_description = query_info['ko_des']
                ko_name = query_info['ko_name']
                pathway_id = [x.split('(')[0].strip() for x in query_info['pathways'].split(';')]
            data = {
                'kegg_id': main_table_id,
                'gene_id': gene,
                'gene_name': gene_name,
                'description': description,
                'ko_id': ko_id,
                # 'genesets': ';'.join(genesets),
                'ko_description': ko_description,
                'first_category': pathways_class1,
                'second_category': pathways_class2,
                'pathway_id': ';'.join(pathway_id),
                'ko_name': ko_name,
            }
            data.update(geneset_buer)
            data_list.append(data)
        try:
            collection = self.db['sg_geneset_kegg_class_gene']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入kegg基因对应关系出错:")
        else:
            self.bind_object.logger.info("导入kegg基因对应关系成功!")

    def add_kegg_regulate_new(self, main_table_id, geneset_id, kegg_stat_xls, gene_kegg_level_table_xls, work_dir):
        # 通过判断传入的geneset_id的个数来确认取数据的位置，确认是一个还是两个基因集，然后现在分情况讨论
        # 以后mongo出现NaN的时候，通过fillna更改的时候，尽量靠近插入mongo库那一步，测试发现二者之间如果还进行读写操作，会导致
        # NaN改不过来的情况
        stat = pd.read_table(kegg_stat_xls, header=0)
        level = pd.read_table(gene_kegg_level_table_xls, header=0)
        stat_level = pd.merge(stat, level, on='Pathway_id')
        stat_level.to_csv(work_dir + "/" + "stat_level", sep = '\t', index=False)

        # 按照kegg官网进行一级分类的排序
        list_custom = ['Metabolism', 'Genetic Information Processing', 'Environmental Information Processing', 'Cellular Processes', 'Organismal Systems', 'Human Diseases', 'Drug Development']
        appended_data = []
        for i in list_custom:
            if i in list(stat_level.first_category):
                data = stat_level.loc[stat_level['first_category']==i]
                appended_data.append(data)

        appended_data = pd.concat(appended_data)
        appended_data.drop(['graph_id', 'hyperlink', 'graph_png_id'], axis=1, inplace=True)
        appended_data.to_csv(work_dir + "/" + "kegg_annotation_analysis", sep = '\t', index=False)

        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                raise Exception('kegg_main_table_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(gene_kegg_level_table_xls):
            raise Exception('gene_kegg_level_table_xls所指定的路径:{}不存在，请检查！'.format(gene_kegg_level_table_xls))

        if len(geneset_id.split(",")) == 1:
            with open(work_dir + "/" + "kegg_annotation_analysis") as f_kaa, open(work_dir + "/" + "kegg_annotation_analysis_new", "w") as fw_kaa:
                header_kaa = f_kaa.readline()
                hkl = header_kaa.strip().split("\t")
                geneko_name1 =  hkl[3][:-1] + "ko"

                fw_kaa.write(hkl[0] + "\t" + hkl[1] + "\t" + hkl[2] + "\t" + geneko_name1 + "\t" + hkl[3] + "\t" + hkl[4] + "\t"
                             + hkl[5] + "\t" +
                             hkl[6] + "\t" + hkl[7] + "\t" + hkl[8] + "\t" + hkl[9] + "\t" + hkl[10] + "\n" )

                for line in f_kaa:
                    genes_1 = []
                    genes_2 = []
                    line_list = line.strip().split("\t")
                    line_list_3 = line_list[3]
                    name1_genes = line_list_3.split(');')
                    genes_1 += [x.split('(')[0] for x in name1_genes]

                    fw_kaa.write(line_list[0] + "\t" + line_list[1] + "\t" + line_list[2] + "\t" + line_list[3] + "\t" +
                                 ",""".join(genes_1) + "\t" + line_list[4] + "\t" + line_list[5] + "\t" + line_list[6] + "\t" +
                                 line_list[7] + "\t" + line_list[8] + "\t" + line_list[9] + "\t" + line_list[10] + "\n" )

            kaa = pd.read_table(work_dir + "/" + "kegg_annotation_analysis_new", sep = '\t', header=0)
            kaa.rename(columns={kaa.columns[0]: "pathway_id", kaa.columns[5]: "link"},inplace=True)
            kaa.fillna("", inplace=True)
            kaa.drop(['seq_list', 'number_of_seqs'], axis=1, inplace=True)
            kaa.to_csv(work_dir + "/" + "kegg_analysis_of_anotate", sep = '\t', index=False)

            with open(work_dir + "/" + "kegg_analysis_of_anotate") as r_kaa, open(work_dir + "/" + "new_data_rkaa", "w") as wkaa:
                head_r_kaa = r_kaa.readline().strip().split("\t")
                geneset_name_r_kaa = head_r_kaa[2].split("numbers")[0].rstrip("_")
                str_name = geneset_name_r_kaa + "_str"
                head_r_kaa.insert(5,str_name)
                wkaa.write("\t".join(head_r_kaa) + "\n")
                for line in r_kaa:
                    line = line.strip().split("\t")
                    new_ele = line[4].split(",")
                    new_ele = str(new_ele)
                    line.insert(5,new_ele)
                    wkaa.write("\t".join(line) + "\n")
            new_data_rkaa = pd.read_table(work_dir + "/" + "new_data_rkaa", header=0, sep="\t")
            new_data_rkaa.rename(columns={new_data_rkaa.columns[1]:"ko_ids",
                                          new_data_rkaa.columns[4]: new_data_rkaa.columns[5], new_data_rkaa.columns[5]: new_data_rkaa.columns[4]},inplace=True)
            new_data_rkaa['kegg_id'] = ObjectId(main_table_id)
            new_data_rkaa.fillna("", inplace=True)
            new_data_rkaa_list = new_data_rkaa.to_dict('records')
            target_col1 = new_data_rkaa.columns[5]
            for each in new_data_rkaa_list:
                each[target_col1] = eval(each[target_col1])
            # kaa['kegg_id'] = ObjectId(main_table_id)
            # kaa_list = kaa.to_dict('records')


            self.create_db_table('sg_geneset_kegg_class_detail', new_data_rkaa_list)
            # self.create_db_table('sg_geneset_kegg_class_detail', kaa_list)
            self.update_db_record('sg_geneset_kegg_class', ObjectId(main_table_id), main_id=ObjectId(main_table_id))

            new_data = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep = '\t', header=0)
            new_data.groupby("second_category")
            group_obj = new_data.groupby("second_category")
            groups = group_obj.groups.keys()
            # 做一个基因集
            genesets = new_data.columns[3].split()
            result = defaultdict(dict)
            for each in groups:
                first = new_data.loc[new_data["second_category"] ==each]['first_category']
                first = first.to_dict().values()[0]
                for geneset in genesets:
                    group_detail = group_obj.get_group(each)
                    genes = list()
                    for g in group_detail[geneset]:

                        if not pd.isnull(g):
                            tmp = g.split(');')
                            genes += [x.split('(')[0] for x in tmp]
                            #genes = [i for i in genes if genes.count(i) == 1]
                            genes = list(set(genes))
                        else:
                            genes = []
                    result[geneset][each] = [len(genes), first]

            try:
                a = pd.DataFrame(result)
                a.reset_index(inplace=True)
                a.rename(columns={a.columns[0]: "second_category"}, inplace=True)
                a.to_csv(work_dir + "/" + "k", sep = '\t', index=False)
                with open(work_dir + "/" + "k") as f1, open(work_dir + "/" + "kegg_statistic", "w") as fw:
                    header = f1.readline()
                    geneset_name1 = header.strip().split("\t")[1]
                    geneset_name1_num = geneset_name1 + "_num"
                    fw.write("first_category" + "\t" + header.strip().split("\t")[0] + "\t" + geneset_name1_num + "\n")
                    for line in f1:
                        line_split = line.strip().split("\t")
                        sec = line_split[0]
                        num1 = line_split[1].strip("[]").split(",")[0]
                        first_cate = line_split[1].strip("[]").split(",")[1].strip().strip("'")
                        fw.write(first_cate + "\t" + sec + "\t" + num1 + "\n")
                df_a = pd.read_table(work_dir + "/" + "kegg_statistic", header=0, sep="\t")

                list_custom = ['Metabolism', 'Genetic Information Processing',
                               'Environmental Information Processing',
                               'Cellular Processes', 'Organismal Systems',
                               'Human Diseases',
                               'Drug Development']
                appended_data_new_1 = []
                for i in list_custom:
                    if i in list(df_a.first_category):
                        data = df_a.loc[df_a['first_category'] == i]
                        appended_data_new_1.append(data)

                appended_data_new_1 = pd.concat(appended_data_new_1)

                appended_data_new_1["kegg_id"] = ObjectId(main_table_id)
                appended_data_new_1['geneset_id'] = ObjectId(geneset_id)
                data_new = appended_data_new_1.to_dict('records')
                appended_data_new_1.to_csv(work_dir + "/" +  "kegg_statistic", sep='\t', index=False)
                # data_new = a.to_dict('records')
                self.create_db_table('sg_geneset_kegg_class_statistic', data_new)

                with open(kegg_stat_xls, 'rb') as r:
                    # 获取numbers和genesets的列
                    first_line = r.readline().strip().split("\t")[2:]
                    # print r.next()
                    genesets_name = []
                    for fl in first_line:
                        if "numbers" in fl:
                            # 获取geneset的name，
                            genesets_name.append(fl[:-8])

                main_collection = self.db['sg_geneset_kegg_class']
                main_collection.update({"_id": ObjectId(main_table_id)}, {"$set": {"table_columns": genesets_name}})
                self.bind_object.logger.info("成功更新kegg主表的基因集名字信息")
                df_b = pd.read_table(work_dir + "/" + "kegg_statistic", header=0, sep="\t")
                df_b.drop_duplicates(['first_category'], inplace=True)
                df_control = pd.DataFrame({'first_category': ['Cellular Processes', 'Human Diseases',
                                                              'Genetic Information Processing',
                                                              'Environmental Information Processing',
                                                              'Organismal Systems', 'Metabolism', 'Drug Development'],
                                           'categories': ['CP', 'HD', 'GIP', 'EIP', 'OS', 'M', 'DD']})
                df_short = pd.merge(df_b, df_control, on="first_category")
                categories = list(df_short['categories'])
                main_collection.update({"_id": ObjectId(main_table_id)}, {"$set": {"categories": categories}})
                self.bind_object.logger.info("成功更新kegg主表的一级分类信息缩写")
            except Exception as e:
                self.bind_object.set_error("导入kegg统计信息出错")
            else:
                self.bind_object.logger.info("导入kegg统计信息成功")

        if len(geneset_id.split(",")) == 2:
            # kaa = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep = '\t', header=0)
            # kaa.rename(columns={kaa.columns[0]: "pathway_id", kaa.columns[6]: "link"},inplace=True)
            # 这个替换如果放在写"kegg_annotation_analysis"这个文件的前面，然后读到这里，填充进去还会是NaN,所以要靠近导表前一步才可以

            with open(work_dir + "/" + "kegg_annotation_analysis") as f_kaa, open(work_dir + "/" + "kegg_annotation_analysis_new", "w") as fw_kaa:
                header_kaa = f_kaa.readline()
                hkl = header_kaa.strip().split("\t")
                geneko_name1 =  hkl[3][:-1] + "ko"
                geneko_name2 = hkl[5][:-1] + "ko"
                fw_kaa.write(hkl[0] + "\t" + hkl[1] + "\t" + hkl[2] + "\t" + geneko_name1 + "\t" + hkl[3] + "\t" + hkl[4] + "\t" +
                geneko_name2 + "\t" + hkl[5] + "\t" + hkl[6] + "\t" + hkl[7] + "\t" + hkl[8] + "\t" + hkl[9] + "\t" + hkl[10] +
                             "\t" + hkl[11] + "\t" + hkl[12] + "\n" )

                for line in f_kaa:
                    genes_1 = []
                    genes_2 = []
                    line_list = line.strip().split("\t")
                    line_list_3 = line_list[3]
                    name1_genes = line_list_3.split(');')
                    genes_1 += [x.split('(')[0] for x in name1_genes]

                    line_list_5 = line_list[5]
                    name2_genes = line_list_5.split(');')
                    genes_2 += [x.split('(')[0] for x in name2_genes]
                    fw_kaa.write(line_list[0] + "\t" + line_list[1] + "\t" + line_list[2] + "\t" + line_list[3] + "\t" +
                                 ",""".join(genes_1) + "\t" +line_list[4] + "\t" + line_list[5] + "\t" + ",".join(genes_2) + "\t" +
                                 line_list[6] + "\t" + line_list[7] + "\t" + line_list[8] + "\t" + line_list[9] + "\t" +
                                 line_list[10] + "\t" + line_list[11] + "\t" + line_list[12] + "\n")

            kaa = pd.read_table(work_dir + "/" + "kegg_annotation_analysis_new", sep = '\t', header=0)
            kaa.rename(columns={kaa.columns[0]: "pathway_id", kaa.columns[8]: "link"},inplace=True)
            kaa.fillna("", inplace=True)
            kaa.drop(['seq_list', 'number_of_seqs'], axis=1, inplace=True)
            kaa.to_csv(work_dir + "/" + "kegg_analysis_of_anotate", sep = '\t', index=False)

            with open(work_dir + "/" + "kegg_analysis_of_anotate") as r_kaa, open(work_dir + "/" + "new_data_rkaa", "w") as wkaa:
                head_r_kaa = r_kaa.readline().strip().split("\t")
                geneset_name_r_kaa_1 = head_r_kaa[2].split("numbers")[0].rstrip("_")
                str_name_1 = geneset_name_r_kaa_1 + "_str"

                geneset_name_r_kaa_2 = head_r_kaa[5].split("numbers")[0].rstrip("_")
                str_name_2 = geneset_name_r_kaa_2 + "_str"
                head_r_kaa.insert(5,str_name_1)
                head_r_kaa.insert(9,str_name_2)

                wkaa.write("\t".join(head_r_kaa) + "\n")
                for line in r_kaa:
                    line = line.strip().split("\t")
                    new_ele_1 = line[4].split(",")
                    new_ele_1 = str(new_ele_1)
                    line.insert(5,new_ele_1)

                    new_ele_2 = line[8].split(",")
                    new_ele_2 = str(new_ele_2)
                    line.insert(9,new_ele_2)
                    wkaa.write("\t".join(line) + "\n")
            new_data_rkaa = pd.read_table(work_dir + "/" + "new_data_rkaa", header=0, sep="\t")
            new_data_rkaa.rename(columns={new_data_rkaa.columns[1]:"ko_ids",
            new_data_rkaa.columns[4]: new_data_rkaa.columns[5],
            new_data_rkaa.columns[5]: new_data_rkaa.columns[4],
            new_data_rkaa.columns[8]: new_data_rkaa.columns[9],
            new_data_rkaa.columns[9]: new_data_rkaa.columns[8]},inplace=True)
            new_data_rkaa['kegg_id'] = ObjectId(main_table_id)
            new_data_rkaa.fillna("", inplace=True)
            target_col1 = new_data_rkaa.columns[5]
            target_col2 = new_data_rkaa.columns[9]
            new_data_rkaa_list = new_data_rkaa.to_dict('records')
            for each in new_data_rkaa_list:
                each[target_col1] = eval(each[target_col1])
                each[target_col2] = eval(each[target_col2])
            self.create_db_table('sg_geneset_kegg_class_detail', new_data_rkaa_list)
            # self.create_db_table('sg_geneset_kegg_class_detail', kaa_list)
            self.update_db_record('sg_geneset_kegg_class', ObjectId(main_table_id), main_id=ObjectId(main_table_id))

            new_data = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep = '\t', header=0)
            self.bind_object.logger.info("开始进行class分类导表")
            new_data.groupby("second_category")
            group_obj = new_data.groupby("second_category")
            groups = group_obj.groups.keys()
            # 做2个基因集
            genesets = new_data.columns[3], new_data.columns[5]
            result = defaultdict(dict)
            for each in groups:
                first = new_data.loc[new_data["second_category"] == each]['first_category']
                first = first.to_dict().values()[0]
                for geneset in genesets:
                    group_detail = group_obj.get_group(each)
                    genes = list()
                    for g in group_detail[geneset]:

                        # isnull支持的数据类型更多，相比isnan
                        if not pd.isnull(g):
                            tmp = g.split(');')
                            genes += [x.split('(')[0] for x in tmp]
                            # 用set会弹出不知名的错误
                            genes = [i for i in genes if genes.count(i) == 1]
                        else:
                            genes = []
                    result[geneset][each] = [len(genes), first]
            # try:
            a = pd.DataFrame(result)
            a.reset_index(inplace=True)
            a.rename(columns={a.columns[0]: "second_category"}, inplace=True)
            a.to_csv(work_dir + "/" + "k", sep = '\t', index=False)
            with open(work_dir + "/" + "k") as f1, open(work_dir + "/" + "kegg_statistic", "w") as fw:
                header = f1.readline()
                geneset_name1 = header.strip().split("\t")[1]
                geneset_name1_num = geneset_name1 + "_num"
                geneset_name2 = header.strip().split("\t")[2]
                geneset_name2_num = geneset_name2 + "_num"
                fw.write("first_category" + "\t" + header.strip().split("\t")[0] + "\t" + geneset_name1_num + "\t" + \
                         geneset_name2_num +"\n")
                for line in f1:
                    line_split = line.strip().split("\t")
                    sec = line_split[0]
                    num1 = line_split[1].strip("[]").split(",")[0]
                    num2 = line_split[2].strip("[]").split(",")[0]
                    first_cate = line_split[1].strip("[]").split(",")[1].strip().strip("'")
                    fw.write(first_cate + "\t" + sec + "\t" + num1 + "\t" + num2 + "\n")
            df_a = pd.read_table(work_dir + "/" + "kegg_statistic", header=0, sep="\t")

            list_custom = ['Metabolism', 'Genetic Information Processing',
                           'Environmental Information Processing',
                           'Cellular Processes', 'Organismal Systems',
                           'Human Diseases',
                           'Drug Development']
            appended_data_new_2 = []
            for i in list_custom:
                if i in list(df_a.first_category):
                    data = df_a.loc[df_a['first_category'] == i]
                    appended_data_new_2.append(data)

            appended_data_new_2 = pd.concat(appended_data_new_2)

            appended_data_new_2["kegg_id"] = ObjectId(main_table_id)
            # appended_data_new_2['geneset_type'] = geneset_type
            # appended_data_new_2['geneset_id'] = ObjectId(geneset_id)
            appended_data_new_2['geneset_id'] = geneset_id


            appended_data_new_2.to_csv(work_dir + "/" +
                                       "kegg_statistic",sep='\t', index=False)
            data_new = appended_data_new_2.to_dict('records')
            # data_new = a.to_dict('records')
            self.create_db_table('sg_geneset_kegg_class_statistic', data_new)
            self.bind_object.logger.info("完成class分类导表")

            with open(kegg_stat_xls, 'rb') as r:
                self.bind_object.logger.info("开始kegg主表的基因集名字信息更新")
                # 获取numbers和genesets的列
                first_line = r.readline().strip().split("\t")[2:]
                # print r.next()
                genesets_name = []
                for fl in first_line:
                    if "numbers" in fl:
                        # 获取geneset的name，
                        genesets_name.append(fl[:-8])

            main_collection = self.db['sg_geneset_kegg_class']
            main_collection.update({"_id": ObjectId(main_table_id)}, {"$set": {"table_columns": genesets_name}})
            self.bind_object.logger.info("成功更新kegg主表的基因集名字信息")
            df_b = pd.read_table(work_dir + "/" + "kegg_statistic", header=0, sep="\t")
            df_b.drop_duplicates(['first_category'], inplace=True)
            df_control = pd.DataFrame({'first_category': ['Cellular Processes', 'Human Diseases',
                                                          'Genetic Information Processing',
                                                          'Environmental Information Processing',
                                                          'Organismal Systems', 'Metabolism', 'Drug Development'],
                                       'categories': ['CP', 'HD', 'GIP', 'EIP', 'OS', 'M', 'DD']})
            df_short = pd.merge(df_b, df_control, on="first_category")
            categories = list(df_short['categories'])
            main_collection.update({"_id": ObjectId(main_table_id)}, {"$set": {"categories": categories}})
            self.bind_object.logger.info("成功更新kegg主表的一级分类信息缩写")
            # except Exception as e:
            #     self.bind_object.logger.info("导入kegg统计信息出错")
            # else:
            self.bind_object.logger.info("导入kegg统计信息成功")

    @report_check
    # 最后通过更新插入kegg_class主表的geneset的名字
    def add_kegg_regulate_detail(self, regulate_id, kegg_regulate_table):
        """
        :param regulate_id: 主表ID
        :param kegg_regulate_table: kegg_stat.xls统计结果文件
        :return:
        """
        # 会导表3张
        main_collection = self.db['sg_geneset_kegg_class']
        kegg_main = self.db['sg_annotation_kegg']
        kegg_level_coll = self.db['sg_annotation_kegg_level']
        if not isinstance(regulate_id, ObjectId):
            if isinstance(regulate_id, types.StringTypes):
                regulate_id = ObjectId(regulate_id)
            else:
                raise Exception('kegg_regulate_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(kegg_regulate_table):
            raise Exception('kegg_regulate_table所指定的路径:{}不存在，请检查！'.format(kegg_regulate_table))
        task_id = main_collection.find_one({"_id": regulate_id})['task_id']
        kegg_main_id = kegg_main.find_one({"task_id": task_id})['_id']
        kegg_result = kegg_level_coll.find({"kegg_id": kegg_main_id})
        path_def = {}
        for kr in kegg_result:
            # print kr['pathway_definition']
            path_def[kr['pathway_id'].split(":")[-1]] = kr['pathway_definition']
        # print path_def
        # print task_id
        data_list = []
        with open(kegg_regulate_table, 'rb') as r:
            # 获取numbers和genesets的列
            first_line = r.readline().strip().split("\t")[2:]
            # print r.next()
            genesets_name = []
            for fl in first_line:
                if "numbers" in fl:
                    # 获取geneset的name，
                    genesets_name.append(fl[:-8])
            for line in r:
                line = line.strip('\n').split('\t')
                # regulate_id是表 "sg_geneset_kegg_class"的_id
                insert_data = {
                    'kegg_id': regulate_id,
                    'pathway_id': line[0],
                    'ko_ids': line[1],
                    # 'pathway_definition': path_def[line[0]],
                    'link': line[-1]
                }
                try:
                    insert_data.update({'pathway_definition': path_def[line[0]]})
                except:
                    insert_data.update({'pathway_definition': ''})
                # print path_def[line[0]]
                for n, gn in enumerate(genesets_name):
                    # 因为ko号之间和基因之间都是;分割，这样会导致'Ko12)'.split("(")没有(的字符串分割得到本身，这一点很容易错；Python的基础split方法也支持多个分隔符
                    # gene_list = line[3+2*n].split(";")
                    gene_list = line[3+2*n].split(");")
                    gene_list = [ x.split("(")[0] for x in gene_list ]
                    insert_data["{}_geneko".format(gn)] = line[3+2*n]
                    insert_data["{}_numbers".format(gn)] = line[2+2*n]
                    insert_data["{}_genes".format(gn)] = gene_list
                    insert_data["{}_str".format(gn)] = ";".join(set(gene_list))
                data_list.append(insert_data)
            try:
                collection = self.db['sg_geneset_kegg_class_detail']
                # main_collection = self.db['sg_geneset_kegg_class']
                collection.insert_many(data_list)
                main_collection.update({"_id": ObjectId(regulate_id)}, {"$set": {"table_columns": genesets_name}})
            except Exception as e:
                self.bind_object.set_error("导入kegg调控统计表：%s信息出错:%s" % (kegg_regulate_table, e))
            else:
                self.bind_object.logger.info("导入kegg调控统计表:%s信息成功!" % kegg_regulate_table)

    def add_geneset(self, diff_work_dir, group_id=None, task_id=None, project_sn=None):
        if not isinstance(group_id, ObjectId):
            if isinstance(group_id, types.StringTypes):
                group_id = ObjectId(group_id)
            else:
                raise Exception('group_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(diff_work_dir):
            raise Exception('diff_work_dir所指的路径:{}不存在'.format(diff_work_dir))
        diff_exp_files = os.listdir(diff_work_dir)
        for f in diff_exp_files:
            if f.startswith('diffcmp'):
                self.bind_object.logger.info("开始导入{}的基因集".format(f))
                compare_all, gene_set_all = self.get_diff_list(diffcmp=diff_work_dir + "/" + f, up_down='all')
                compare_up, gene_set_up = self.get_diff_list(diffcmp=diff_work_dir + "/" + f, up_down='up')
                compare_down, gene_set_down = self.get_diff_list(diffcmp=diff_work_dir + "/" + f, up_down='down')
                if len(gene_set_all) >=1:
                    data_all = {
                        'group_id': group_id,
                        'task_id': task_id,
                        'name': compare_all,
                        'desc': '{}基因集'.format(compare_all),
                        'project_sn': project_sn,
                        'geneset_length': len(gene_set_all),
                        'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                        'is_use': 1,
                    }
                    try:
                        collection = self.db["sg_geneset"]
                        geneset_up_id_all = collection.insert_one(data_all).inserted_id
                        self.update_db_record('sg_geneset', ObjectId(geneset_up_id_all), main_id=ObjectId(geneset_up_id_all))
                        if geneset_up_id_all:
                            self.add_geneset_detail(geneset_up_id_all, gene_set_all)
                    except Exception as e:
                        self.bind_object.set_error("导入基因集集：%s信息出错:%s" % (f, e))
                    else:
                        self.bind_object.logger.info("导入基因集：%s信息成功!" % f)
                if len(gene_set_up) >=1:
                    data_up = {
                        'group_id': group_id,
                        'task_id': task_id,
                        'name': compare_up,
                        'desc': '{}基因集'.format(compare_up),
                        'project_sn': project_sn,
                        'geneset_length': len(gene_set_up),
                        'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                        'is_use': 1,
                    }
                    try:
                        collection = self.db["sg_geneset"]
                        geneset_up_id_up = collection.insert_one(data_up).inserted_id
                        self.update_db_record('sg_geneset', ObjectId(geneset_up_id_up), main_id=ObjectId(geneset_up_id_up))
                        if geneset_up_id_up:
                            self.add_geneset_detail(geneset_up_id_up, gene_set_up)
                    except Exception as e:
                        self.bind_object.set_error("导入基因集集：%s信息出错:%s" % (f, e))
                    else:
                        self.bind_object.logger.info("导入基因集：%s信息成功!" % f)
                if len(gene_set_down) >=1:
                    data_down = {
                        'group_id': group_id,
                        'task_id': task_id,
                        'name': compare_down,
                        'desc': '{}基因集'.format(compare_down),
                        'project_sn': project_sn,
                        'geneset_length': len(gene_set_down),
                        'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                        'is_use': 1,
                    }
                    try:
                        collection = self.db["sg_geneset"]
                        geneset_up_id_down = collection.insert_one(data_down).inserted_id
                        self.update_db_record('sg_geneset', ObjectId(geneset_up_id_down), main_id=ObjectId(geneset_up_id_down))
                        if geneset_up_id_down:
                            self.add_geneset_detail(geneset_up_id_down, gene_set_down)
                    except Exception as e:
                        self.bind_object.set_error("导入基因集集：%s信息出错:%s" % (f, e))
                    else:
                        self.bind_object.logger.info("导入基因集：%s信息成功!" % f)
            else:
                continue

    def add_geneset_detail(self, geneset_id, sequence):
        if not isinstance(geneset_id, ObjectId):
            if isinstance(geneset_id, types.StringTypes):
                geneset_id = ObjectId(geneset_id)
            else:
                raise Exception('geneset_id必须为ObjectId对象或其对应的字符串！')
        data = [
            ("geneset_id", ObjectId(geneset_id)),
            ("seq_list", sequence)
        ]
        data = SON(data)
        try:
            collection = self.db["sg_geneset_detail"]
            collection.insert_one(data)
        except Exception as e:
            self.bind_object.set_error("导入基因集detail表出错:%s" % e)
        else:
            self.bind_object.logger.info("导入基因集detail表成功!")

    def get_diff_list(self, diffcmp, up_down=None):
        if not os.path.exists(diffcmp):
            raise Exception("{}文件不存在，无法对up和down差异基因进行分类！".format(diffcmp))
        with open(diffcmp, 'r') as f1:
            header = f1.readline()
            sequence = []
            for lines in f1:
                line = lines.strip().split("\t")
                seq_id = line[0]
                significant = line[-3]
                regulate = line[-4]
                cmp = line[-1].split("|")
                if significant == 'yes':
                    if up_down == 'all' and (regulate == 'up' or regulate == 'down'):
                        sequence.append(seq_id)
                    if up_down == 'up' and regulate == 'up':
                        sequence.append(seq_id)
                    if up_down == 'down' and regulate == 'down':
                        sequence.append(seq_id)
                else:
                    pass
        compare = '{}_vs_{}_{}'.format(cmp[0], cmp[1], up_down)
        if sequence:
            return compare, sequence
        else:
            return compare, sequence

    def insert_genelist_forquery(self, gene_list, main_id):
        seq_dict = {
            'genelist_id': main_id,
            'seq_list': gene_list
        }
        try:
            collection = self.db["sg_geneset_querylist"]
            collection.insert_one(seq_dict)
        except Exception as e:
            self.bind_object.set_error("挑选出来的基因列表导入sg_geneset_querylist表出错:%s" % e)
        else:
            self.bind_object.logger.info("挑选出来的基因列表导入sg_geneset_querylist表成功!")

#     -------为了原核基因集中基因与注释关系对应表编写的通用函数，可以获得特定基因在query表里的注释信息-----
    def get_gene_query_info(self, gene_id):
        """为了给差异分析结果添加gene_name"""
        task_id = self.bind_object.sheet.option('task_id')
        type_ = self.bind_object.sheet.option('type')
        # print(task_id,type_)
        query_table = self.db['sg_annotation_query']
        query_main = query_table.find_one({"task_id": task_id, "type": type_})
        if not query_main:
            self.bind_object.set_error("找不到对应的query表",task_id,type_)
        if "main_id" not in query_main:
            query_main_id = query_main['_id']
        else:
            query_main_id = query_main['main_id']
        query_detail = self.db['sg_annotation_query_detail']
        query_info = ''
        query_info = query_detail.find_one({"query_id": query_main_id, "gene_id" : gene_id})
        return query_info


    def add_kegg_enrich_pic(self, main_table_id, level_path, png_dir, source):
        # 导入图片信息数据
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                kegg_id = ObjectId(main_table_id)
            else:
                raise Exception('main_table_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(level_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(level_path))
        if not os.path.exists(png_dir):
            raise Exception('{}所指定的路径不存在，请检查！'.format(png_dir))
        data_list = []
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                pid = re.sub('path:', '', line[0])
                if os.path.exists(png_dir + '/' + line[0] + '.html.mark'):
                    with open(png_dir + '/' + line[0] + '.html.mark', 'r') as mark_f:
                        for line_mark in mark_f.readlines():
                            # print len(line_mark.strip("\n").split("\t"))
                            if len(line_mark.strip("\n").split("\t")) >= 10:
                                [png, shape, bg_color, fg_color, coords, title, kos, href, gene_list, geneset_list] = line_mark.strip("\n").split("\t")
                                title = title.replace("\\n", "\n")
                            else:
                                continue

                            if bg_color == "":
                                continue

                            insert_data = {
                                'kegg_enrich_id': kegg_id,
                                'pathway_id': line[0],
                                'shape': shape,
                                'bg_colors': bg_color,
                                'fg_colors': fg_color,
                                'coords': coords,
                                'href': href,
                                'kos': kos,
                                'title': title
                            }

                            if source == "diff_exp":
                                geneset_reg_list = geneset_list.split("|")
                                geneset_list = ["_".join(x.split("_")[:-1]) for x in geneset_reg_list]
                                reg_list = [x.split("_")[-1] for x in geneset_reg_list]
                                insert_data.update({
                                    'gene_list': gene_list.split("|"),
                                    'geneset_list': geneset_list,
                                    'reg_list': reg_list
                                })
                            else:
                                insert_data.update({
                                    'gene_list': gene_list.split("|"),
                                    'geneset_list': geneset_list.split("|")
                                })

                            if bg_color != "" and len(bg_color.split(",")) > 0:
                                insert_data.update({'bg_type': len(bg_color.split(","))})
                            if fg_color != "" and len(fg_color.split(",")) > 0:
                                insert_data.update({'fg_type': len(fg_color.split(","))})
                            data_list.append(insert_data)

                else:
                    self.bind_object.logger.info("kegg 图片{} 不存在html标记!".format(line[0]))

        if data_list:
            try:
                collection = self.db['sg_geneset_kegg_enrich_pic']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.set_error("导入kegg注释图片信息：%s、%s出错!" % (level_path, png_dir))
            else:
                self.bind_object.logger.info("导入kegg注释图片信息：%s、%s 成功!" % (level_path, png_dir))

if __name__ == "__main__":
    pass

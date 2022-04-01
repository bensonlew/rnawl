# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# last_modify:20161205

import os
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import pandas as pd
import json
import re
import gridfs
import glob
from collections import defaultdict
from mbio.api.database.denovo_rna_v2.api_base import ApiBase
import math


class RefRnaV2Geneset(ApiBase):
    def __init__(self, bind_object):
        super(RefRnaV2Geneset, self).__init__(bind_object)
        self._project_type = 'smallrna'
        #self._db_name = Config().MONGODB + '_ref_rna'

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
            "status": "start",
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
        with open(geneset_cog_table, 'r') as f:
            first_line = f.readline().strip().split("\t")
            print(first_line)
            # print f.next().split("\t")

            for gn in first_line[2:]:
                if "list" in gn:
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
                for n, gn in enumerate(geneset_name):
                    data[gn + "_cog"] = int(line[2*n+2])
                    # data[gn + "_nog"] = int(line[6*n+3])
                    # data[gn + "_kog"] = int(line[6*n+4])
                    data[gn + "_cog_list"] = line[2*n+3].split(";")
                    # data[gn + "_nog_list"] = line[6*n+6].split(";")
                    # data[gn + "_kog_list"] = line[6*n+7].split(";")
                    # data[gn + "_cog_str"] = line[6*n+5]
                    # data[gn + "_nog_str"] = line[6*n+6]
                    # data[gn + "_kog_str"] = line[6*n+7]
                data_list.append(data)
        try:
            collection = self.db['sg_geneset_cog_class_detail']
            main_collection = self.db['sg_geneset_cog_class']
            if len(data_list) != 0:
                collection.insert_many(data_list)
            main_collection.update({"_id": ObjectId(geneset_cog_id)}, {"$set": {"table_columns": geneset_name, "status": "end"}})
            self.bind_object.logger.info(geneset_name)
        except Exception as e:
            self.bind_object.set_error("导入cog表格：%s出错:%s" % (geneset_cog_table, e))
        else:
            self.bind_object.logger.info("导入cog表格：%s成功!" % (geneset_cog_table))

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
        with open(go_enrich_dir, 'r') as f:
            f.readline()
            for line in f:
                line = line.strip().split('\t')
                if line[2] == 'p':
                    continue
                data = [
                    ('go_enrich_id', go_enrich_id),
                    ('go_id', line[0]),
                    ('go_type', line[1]),
                    ('enrichment', line[2]),
                    ('discription', line[3]),
                    ('ratio_in_study', line[4]),
                    ('ratio_in_pop', line[5]),
                    ('p_uncorrected', float(line[6])),
                    ('p_corrected', float(line[9])),
                    ('enrich_factor', float(line[4].split("/")[0]) / float(line[5].split("/")[0])),
                    ('depth', int(line[7])),
                    ('study_count', int(line[4].split("/")[0])),
                    ('pop_count', int(line[5].split("/")[0])),
                    ('seq_list', line[-1]),
                    ('seq_str', line[-1].split(";"))
                ]
                go_type.append(line[1])
                data = SON(data)
                data_list.append(data)
        if data_list:
            # 插入-logpvalue -logpcorrected 值相关字段
            pvalues = [dict(son)['p_uncorrected'] for son in data_list]
            if len([x for x in pvalues if x > 0]) > 0:
                pvalues_min = min([x for x in pvalues if x > 0]) / 10
            else:
                pvalues_min = 0.0001
            pvalues_min = - math.log10(pvalues_min)
            log10x = [-math.log10(x) if x > 0 else pvalues_min for x in pvalues]
            for i in range(0, len(log10x)):
                data_list[i]['neg_log10p_uncorrected'] = log10x[i]

            pvalues = [dict(son)['p_corrected'] for son in data_list]
            if len([x for x in pvalues if x > 0]) > 0:
                pvalues_min = min([x for x in pvalues if x > 0]) / 10
            else:
                pvalues_min = 0.0001
            pvalues_min = - math.log10(pvalues_min)
            log10x = [-math.log10(x) if x > 0 else min(pvalues) / 10 for x in pvalues]
            for i in range(0, len(log10x)):
                data_list[i]['neg_log10p_corrected'] = log10x[i]

            try:
                collection = self.db['sg_geneset_go_enrich_detail']
                collection.insert_many(data_list)
                coll = self.db['sg_geneset_go_enrich']
                go_type = list(set(go_type))
                coll.update({'_id': go_enrich_id}, {'$set': {'categories': go_type}})
                # main_collection = self.db['sg_geneset_go_enrich']
                # main_collection.update({"_id": ObjectId(go_enrich_id)}, {"$set": {"status": "end"}})
            except Exception as e:
                print("导入go富集信息：%s出错:%s" % (go_enrich_dir, e))
            else:
                print("导入go富集信息：%s成功!" % (go_enrich_dir))
        else:
            raise Exception('GO富集没有结果')

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
    def add_kegg_enrich_detail(self, enrich_id, kegg_enrich_table, geneset_list_path, all_list_path):
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
        with open(kegg_enrich_table, 'rb') as r:
            for line in r:
                if re.match(r'\w', line):
                    line = line.strip('\n').split('\t')
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
                        'hyperlink': line[9],
                        'kegg_type': "".join([x[0] for x in line[11].split(' ')])
                    }
                    kegg_type1.append("".join([x[0] for x in line[11].split(' ')]))
                    data_list.append(insert_data)
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

                    coll.update({'_id': enrich_id}, {'$set': {'categories': kegg_type1}})
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

    @report_check
    def add_go_regulate_detail(self, go_regulate_dir, go_regulate_id):
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
                        data["{}_genes".format(dk)] = line[5+n*3].split(";")
                    except:
                        data["{}_str".format(dk)] = ""
                        data["{}_genes".format(dk)] = ""
                    if len(line4) > 1:
                        data["{}_percent_str".format(dk)] = line4[1][:-1]
                    else:
                        data["{}_percent_str".format(dk)] = 0
                data_list.append(data)
        try:
            collection = self.db['sg_geneset_go_class_detail']
            main_collection = self.db['sg_geneset_go_class']
            if len(data_list) != 0:
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
        files = os.listdir(pathway_dir)
        fs = gridfs.GridFS(self.db)
        for f in files:
            png_id = fs.put(open(os.path.join(pathway_dir, f), 'rb'))
            insert_data = {
                'kegg_id': regulate_id,
                'pathway_png': png_id,
                'pathway_id': f.split('.pdf')[0]
            }
            data_list.append(insert_data)
        try:
            collection = self.db['sg_geneset_kegg_class_pathway']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.info("导入kegg调控pathway：%s信息出错:%s" % (pathway_dir, e))
        else:
            self.bind_object.logger.info("导入kegg调控pathway:%s信息成功!" % pathway_dir)

    @report_check
    def add_kegg_regulate_detail(self, regulate_id, kegg_regulate_table):
        """

        :param regulate_id: 主表ID
        :param kegg_regulate_table: kegg_stat.xls统计结果文件
        :return:
        """
        if not isinstance(regulate_id, ObjectId):
            if isinstance(regulate_id, types.StringTypes):
                regulate_id = ObjectId(regulate_id)
            else:
                raise Exception('kegg_regulate_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(kegg_regulate_table):
            raise Exception('kegg_regulate_table所指定的路径:{}不存在，请检查！'.format(kegg_regulate_table))
        data_list = []
        with open(kegg_regulate_table, 'rb') as r:
            first_line = r.readline().strip().split("\t")[2:]
            print(r.next())
            genesets_name = []
            for fl in first_line:
                if "numbers" in fl:
                    genesets_name.append(fl[:-8])
            print(genesets_name)
            print(first_line)
            for line in r:
                line = line.strip('\n').split('\t')
                insert_data = {
                    'kegg_id': regulate_id,
                    'pathway_id': line[0],
                    'ko_ids': line[1]
                }
                # print line
                for n, gn in enumerate(genesets_name):
                    gene_list = re.findall(r"(.*?)\(.*?\);", line[3+2*n])
                    insert_data["{}_numbers".format(gn)] = line[2+2*n]
                    insert_data["{}_genes".format(gn)] = gene_list
                    insert_data["{}_str".format(gn)] = ";".join(gene_list)
                data_list.append(insert_data)
            try:
                collection = self.db['sg_geneset_kegg_class_detail']
                main_collection = self.db['sg_geneset_kegg_class']
                collection.insert_many(data_list)
                main_collection.update({"_id": ObjectId(regulate_id)}, {"$set": {"table_columns": genesets_name}})
            except Exception as e:
                self.bind_object.logger.info("导入kegg调控统计表：%s信息出错:%s" % (kegg_regulate_table, e))
            else:
                self.bind_object.logger.info("导入kegg调控统计表:%s信息成功!" % kegg_regulate_table)

    def add_kegg_regulate_new(self, main_table_id, geneset_id, kegg_stat_xls, gene_kegg_level_table_xls, work_dir,
                              geneset_type):
        # 通过判断传入的geneset_id的个数来确认取数据的位置，确认是一个还是两个基因集，然后现在分情况讨论
        # 以后mongo出现NaN的时候，通过fillna更改的时候，尽量靠近插入mongo库那一步，测试发现二者之间如果还进行读写操作，会导致
        # NaN改不过来的情况
        stat = pd.read_table(kegg_stat_xls, header=0)
        level = pd.read_table(gene_kegg_level_table_xls, header=0)
        stat_level = pd.merge(stat, level, on='Pathway_id')
        stat_level.to_csv(work_dir + "/" + "stat_level", sep='\t', index=False)

        # 按照kegg官网进行一级分类的排序
        list_custom = ['Metabolism', 'Genetic Information Processing', 'Environmental Information Processing',
                       'Cellular Processes', 'Organismal Systems', 'Human Diseases',
                       'Drug Development']
        appended_data = []
        for i in list_custom:
            if i in list(stat_level.first_category):
                data = stat_level.loc[stat_level['first_category'] == i]
                appended_data.append(data)

        appended_data = pd.concat(appended_data)
        appended_data.drop(['graph_id', 'hyperlink', 'graph_png_id'], axis=1, inplace=True)
        appended_data.to_csv(work_dir + "/" + "kegg_annotation_analysis", sep='\t', index=False)

        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                raise Exception('kegg_main_table_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(gene_kegg_level_table_xls):
            raise Exception('gene_kegg_level_table_xls所指定的路径:{}不存在，请检查！'.format(gene_kegg_level_table_xls))

        if len(geneset_id.split(",")) == 1:
            with open(work_dir + "/" + "kegg_annotation_analysis") as f_kaa, open(
                    work_dir + "/" + "kegg_annotation_analysis_new", "w") as fw_kaa:
                header_kaa = f_kaa.readline()
                hkl = header_kaa.strip().split("\t")
                geneko_name1 = hkl[3][:-1] + "ko"

                fw_kaa.write(
                    hkl[0] + "\t" + hkl[1] + "\t" + hkl[2] + "\t" + geneko_name1 + "\t" + hkl[3] + "\t" + hkl[4] + "\t"
                    + hkl[5] + "\t" +
                    hkl[6] + "\t" + hkl[7] + "\t" + hkl[8] + "\t" + hkl[9] + "\t" + hkl[10] + "\n")

                for line in f_kaa:
                    genes_1 = []
                    genes_2 = []
                    line_list = line.strip().split("\t")
                    line_list_3 = line_list[3]
                    name1_genes = line_list_3.split(');')
                    genes_1 += [x.split('(')[0] for x in name1_genes]

                    fw_kaa.write(line_list[0] + "\t" + line_list[1] + "\t" + line_list[2] + "\t" + line_list[3] + "\t" +
                                 ",""".join(genes_1) + "\t" + line_list[4] + "\t" + line_list[5] + "\t" + line_list[
                                     6] + "\t" +
                                 line_list[7] + "\t" + line_list[8] + "\t" + line_list[9] + "\t" + line_list[10] + "\n")

            kaa = pd.read_table(work_dir + "/" + "kegg_annotation_analysis_new", sep='\t', header=0)
            kaa.rename(columns={kaa.columns[0]: "pathway_id", kaa.columns[5]: "link"}, inplace=True)
            kaa.fillna("", inplace=True)
            kaa.drop(['seq_list', 'number_of_seqs'], axis=1, inplace=True)
            kaa.to_csv(work_dir + "/" + "kegg_analysis_of_anotate", sep='\t', index=False)

            with open(work_dir + "/" + "kegg_analysis_of_anotate") as r_kaa, open(work_dir + "/" + "new_data_rkaa",
                                                                                  "w") as wkaa:
                head_r_kaa = r_kaa.readline().strip().split("\t")
                geneset_name_r_kaa = head_r_kaa[2].split("numbers")[0].rstrip("_")
                str_name = geneset_name_r_kaa + "_str"
                head_r_kaa.insert(5, str_name)
                wkaa.write("\t".join(head_r_kaa) + "\n")
                for line in r_kaa:
                    line = line.strip().split("\t")
                    new_ele = line[4].split(",")
                    new_ele = str(new_ele)
                    line.insert(5, new_ele)
                    wkaa.write("\t".join(line) + "\n")
            new_data_rkaa = pd.read_table(work_dir + "/" + "new_data_rkaa", header=0, sep="\t")
            new_data_rkaa.rename(columns={new_data_rkaa.columns[1]: "ko_ids",
                                          new_data_rkaa.columns[4]: new_data_rkaa.columns[5],
                                          new_data_rkaa.columns[5]: new_data_rkaa.columns[4]}, inplace=True)
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

            new_data = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep='\t', header=0)
            new_data.groupby("second_category")
            group_obj = new_data.groupby("second_category")
            groups = group_obj.groups.keys()
            # 做一个基因集
            genesets = new_data.columns[3].split()
            result = defaultdict(dict)
            for each in groups:
                first = new_data.loc[new_data["second_category"] == each]['first_category']
                first = first.to_dict().values()[0]
                for geneset in genesets:
                    group_detail = group_obj.get_group(each)
                    genes = list()
                    for g in group_detail[geneset]:

                        if not pd.isnull(g):
                            tmp = g.split(');')
                            genes += [x.split('(')[0] for x in tmp]
                            genes = [i for i in genes if genes.count(i) == 1]
                        else:
                            genes = []
                    result[geneset][each] = [len(genes), first]

            try:
                a = pd.DataFrame(result)
                a.reset_index(inplace=True)
                a.rename(columns={a.columns[0]: "second_category"}, inplace=True)
                a.to_csv(work_dir + "/" + "k", sep='\t', index=False)
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
                appended_data_new_1['geneset_type'] = geneset_type
                appended_data_new_1['geneset_id'] = ObjectId(geneset_id)
                data_new = appended_data_new_1.to_dict('records')
                appended_data_new_1.to_csv(work_dir + "/" + "kegg_statistic", sep='\t', index=False)
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
                self.bind_object.logger.info("导入kegg统计信息出错")
            else:
                self.bind_object.logger.info("导入kegg统计信息成功")

        if len(geneset_id.split(",")) == 2:
            # kaa = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep = '\t', header=0)
            # kaa.rename(columns={kaa.columns[0]: "pathway_id", kaa.columns[6]: "link"},inplace=True)
            # 这个替换如果放在写"kegg_annotation_analysis"这个文件的前面，然后读到这里，填充进去还会是NaN,所以要靠近导表前一步才可以

            with open(work_dir + "/" + "kegg_annotation_analysis") as f_kaa, open(
                    work_dir + "/" + "kegg_annotation_analysis_new", "w") as fw_kaa:
                header_kaa = f_kaa.readline()
                hkl = header_kaa.strip().split("\t")
                geneko_name1 = hkl[3][:-1] + "ko"
                geneko_name2 = hkl[5][:-1] + "ko"
                fw_kaa.write(hkl[0] + "\t" + hkl[1] + "\t" + hkl[2] + "\t" + geneko_name1 + "\t" + hkl[3] + "\t" + hkl[
                    4] + "\t" +
                             geneko_name2 + "\t" + hkl[5] + "\t" + hkl[6] + "\t" + hkl[7] + "\t" + hkl[8] + "\t" + hkl[
                                 9] + "\t" + hkl[10] +
                             "\t" + hkl[11] + "\t" + hkl[12] + "\n")

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
                                 ",""".join(genes_1) + "\t" + line_list[4] + "\t" + line_list[5] + "\t" + ",".join(
                        genes_2) + "\t" +
                                 line_list[6] + "\t" + line_list[7] + "\t" + line_list[8] + "\t" + line_list[9] + "\t" +
                                 line_list[10] + "\t" + line_list[11] + "\t" + line_list[12] + "\n")

            kaa = pd.read_table(work_dir + "/" + "kegg_annotation_analysis_new", sep='\t', header=0)
            kaa.rename(columns={kaa.columns[0]: "pathway_id", kaa.columns[8]: "link"}, inplace=True)
            kaa.fillna("", inplace=True)
            kaa.drop(['seq_list', 'number_of_seqs'], axis=1, inplace=True)
            kaa.to_csv(work_dir + "/" + "kegg_analysis_of_anotate", sep='\t', index=False)

            with open(work_dir + "/" + "kegg_analysis_of_anotate") as r_kaa, open(work_dir + "/" + "new_data_rkaa",
                                                                                  "w") as wkaa:
                head_r_kaa = r_kaa.readline().strip().split("\t")
                geneset_name_r_kaa_1 = head_r_kaa[2].split("numbers")[0].rstrip("_")
                str_name_1 = geneset_name_r_kaa_1 + "_str"

                geneset_name_r_kaa_2 = head_r_kaa[5].split("numbers")[0].rstrip("_")
                str_name_2 = geneset_name_r_kaa_2 + "_str"
                head_r_kaa.insert(5, str_name_1)
                head_r_kaa.insert(9, str_name_2)

                wkaa.write("\t".join(head_r_kaa) + "\n")
                for line in r_kaa:
                    line = line.strip().split("\t")
                    new_ele_1 = line[4].split(",")
                    new_ele_1 = str(new_ele_1)
                    line.insert(5, new_ele_1)

                    new_ele_2 = line[8].split(",")
                    new_ele_2 = str(new_ele_2)
                    line.insert(9, new_ele_2)
                    wkaa.write("\t".join(line) + "\n")
            new_data_rkaa = pd.read_table(work_dir + "/" + "new_data_rkaa", header=0, sep="\t")
            new_data_rkaa.rename(columns={new_data_rkaa.columns[1]: "ko_ids",
                                          new_data_rkaa.columns[4]: new_data_rkaa.columns[5],
                                          new_data_rkaa.columns[5]: new_data_rkaa.columns[4],
                                          new_data_rkaa.columns[8]: new_data_rkaa.columns[9],
                                          new_data_rkaa.columns[9]: new_data_rkaa.columns[8]}, inplace=True)
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

            new_data = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep='\t', header=0)
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
            a.to_csv(work_dir + "/" + "k", sep='\t', index=False)
            with open(work_dir + "/" + "k") as f1, open(work_dir + "/" + "kegg_statistic", "w") as fw:
                header = f1.readline()
                geneset_name1 = header.strip().split("\t")[1]
                geneset_name1_num = geneset_name1 + "_num"
                geneset_name2 = header.strip().split("\t")[2]
                geneset_name2_num = geneset_name2 + "_num"
                fw.write("first_category" + "\t" + header.strip().split("\t")[0] + "\t" + geneset_name1_num + "\t" + \
                         geneset_name2_num + "\n")
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
            appended_data_new_2['geneset_type'] = geneset_type
            # appended_data_new_2['geneset_id'] = ObjectId(geneset_id)
            appended_data_new_2['geneset_id'] = geneset_id

            appended_data_new_2.to_csv(work_dir + "/" +
                                       "kegg_statistic", sep='\t', index=False)
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

#------------------蛋白相关复制过来的函数---------------------------

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
            self.bind_object.logger.info("主表id：{} 蛋白集：{}".format(main_table_id, geneset_name))
            main_collection.update({"_id": main_table_id},
                                   {"$set": {"table_columns": geneset_name, "status": "end"}})
        except Exception as e:
            raise Exception("导入ipath：%s出错!" % (ipath_input))
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
        circ['circ_id'] = main_table_id
        row_dict_list = circ.to_dict('records')
        main_collection = self.db['sg_geneset_circ']

        try:
            self.create_db_table('sg_geneset_circ_detail', row_dict_list)
        except Exception as e:
            raise Exception("导入main: %s出错!" % (circ_input))
        else:
            self.bind_object.logger.info("导入circ_detail：%s出错!" % (circ_input))

    def add_circ_graph(self, main_table_id, circ_input, enrich_type):
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(circ_input):
            raise Exception('{}所指定的路径不存在，请检查！'.format(ipath_input))

        circ = pd.read_table(circ_input, header=0)
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
            raise Exception("导入circ_graph：%s出错!" % (circ_input))
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

            raise Exception("更新circ主表：%s出错!" % (circ_zscore_input))
        else:
            self.bind_object.logger.info("导入ipath：%s出错!" % (main_table_id))

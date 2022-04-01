# -*- coding: utf-8 -*-
# __author__ = 'qindanhua, qinjincheng'

from mbio.api.database.small_rna.api_base import ApiBase
from biocluster.api.database.base import report_check
import os
from bson.objectid import ObjectId
import pandas as pd
import types
from collections import defaultdict
import re
import unittest

class GenesetClass(ApiBase):
    def __init__(self, bind_object):
        super(GenesetClass, self).__init__(bind_object)

    @report_check
    def add_geneset_cog_class(self, cog_class_table, main_id):
        # check cog_class_table
        if not os.path.exists(cog_class_table):
            self.bind_object.set_error('{} does not exist, abord'.format(go_class_table))
        # transform main_id
        main_id = ObjectId(main_id)
        data_list = list()
        geneset_name = list()
        list_key = list()
        with open(cog_class_table) as f:
            first_line = f.readline().strip().split('\t')
            for cell in first_line[2:]:
                if 'list' in cell:
                    continue
                elif not cell[:-4] in geneset_name:
                    geneset_name.append(cell[:-4])
            for line in f:
                line = line.strip("\n").split('\t')
                data = {
                    'geneset_cog_id': main_id,
                    'type': line[0],
                    'function_categories': line[1]
                }
                for n, gn in enumerate(geneset_name):
                    data['{}_cog'.format(gn)] = int(line[2*n+2])
                    data['{}_cog_list'.format(gn)] = line[2*n+3].split(';')
                    if data['{}_cog_list'.format(gn)] == [""]:
                        data['{}_cog_list'.format(gn)] = []
                    list_key.append('{}_cog_list'.format(gn))
                    data['{}_cog_str'.format(gn)] = line[2*n+3]
                data_list.append(data)

        # add a key-value pair contain all name in geneset
        for idx, each in enumerate(data_list):
            list_all = list()
            for k in set(list_key):
                list_all.extend(each[k])
            list_all = list(set(list_all))
            str_all = ';'.join(list_all)
            data_list[idx]['list_all'] = list_all
            data_list[idx]['str_all'] = str_all

        # try to insert documents to sg_geneset_cog_class_detail
        try:
            self.create_db_table('sg_geneset_cog_class_detail', data_list)
            self.bind_object.logger.info('succeed in creating documents in sg_geneset_cog_class_detail')
        except Exception as e:
            self.bind_object.set_error('fail to create documents with exception: {}'.format(e))
        # try to update document in sg_geneset_cog_class
        try:
            self.update_db_record('sg_geneset_cog_class', main_id, table_columns=list(geneset_name), status='end')
            self.bind_object.logger.info('succeed in updating document in sg_geneset_cog_class')
        except Exception as e:
            self.bind_object.set_error('fail to update document with exception: {}'.format(e))

    @report_check
    def add_geneset_go_class(self, go_class_table, main_id):
        # check go_class_table
        if not os.path.exists(go_class_table):
            self.bind_object.set_error('{} does not exist, abord'.format(go_class_table))
        # transform main_id
        main_id = ObjectId(main_id)
        # prepare data_list
        data_list = list()
        list_key = list()
        with open(go_class_table) as f:
            first_line = f.readline().strip().split('\t')
            doc_keys = list()
            for cell in first_line[3:]:
                name = cell.split(" ")[0]
                if name not in doc_keys:
                    doc_keys.append(name)
            geneset_name = set(doc_keys)
            for line in f:
                line = line.strip().split('\t')
                data = {
                    'go_regulate_id': main_id,
                    'go_type': line[0],
                    'go': line[1],
                    'go_id': line[2]
                }
                for n, dk in enumerate(doc_keys):
                    pc_cell = line[4+n*3].split('(')
                    data['{}_num'.format(dk)] = int(line[3+n*3])
                    data['{}_percent'.format(dk)] = float(pc_cell[0])
                    try:
                        data['{}_str'.format(dk)] = line[5+n*3]
                        data['{}_genes'.format(dk)] = line[5+n*3].split(';')
                        if data['{}_genes'.format(dk)] == [""]:
                            data['{}_genes'.format(dk)] = []
                    except:
                        data['{}_str'.format(dk)] = str()
                        data['{}_genes'.format(dk)] = []
                    list_key.append('{}_genes'.format(dk))
                    if len(pc_cell) > 1:
                        data['{}_percent_str'.format(dk)] = pc_cell[1][:-1]
                    else:
                        data['{}_percent_str'.format(dk)] = 0
                data_list.append(data)

        # add a key-value pair contain all name in geneset
        for idx, each in enumerate(data_list):
            list_all = list()
            for k in set(list_key):
                list_all.extend(each[k])
            list_all = list(set(list_all))
            str_all = ';'.join(list_all)
            data_list[idx]['list_all'] = list_all
            data_list[idx]['str_all'] = str_all

        # try to insert documents to sg_geneset_go_class_detail
        try:
            self.create_db_table('sg_geneset_go_class_detail', data_list)
            self.bind_object.logger.info('succeed in creating documents in sg_geneset_go_class_detail')
        except Exception as e:
            self.bind_object.set_error('fail to create documents with exception: {}'.format(e))
        # try to update document in sg_geneset_go_class
        try:
            self.update_db_record('sg_geneset_go_class', main_id, table_columns=list(geneset_name), status='end')
            self.bind_object.logger.info('succeed in updating document in sg_geneset_go_class')
        except Exception as e:
            self.bind_object.set_error('fail to update document with exception: {}'.format(e))

    @report_check
    def add_kegg_regulate_new(self, main_table_id, geneset_id, kegg_stat_xls, gene_kegg_level_table_xls, work_dir, geneset_type):
        '''
        通过判断传入的geneset_id的个数来确认取数据的位置，确认是一个还是两个基因集，然后现在分情况讨论，
        以后mongo出现NaN的时候，通过fillna更改的时候，尽量靠近插入mongo库那一步，
        测试发现二者之间如果还进行读写操作，会导致NaN改不过来的情况
        '''
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
                    new_ele = list(set(line[4].split(",")))
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
            list_key = list()

            # important change: *_str = "GENE1,GENE2,..." --> *_str = "GENE1;GENE2;..."
            # important change: anno_type = "G" --> geneset_type = "G"
            for index, record in enumerate(new_data_rkaa_list):
                for key in record.keys():
                    if key[-4:] == '_str':
                        new_data_rkaa_list[index][key] = record[key].replace(',', ';')
                    if key == 'anno_type':
                        new_data_rkaa_list[index]['geneset_type'] = 'G'
                        del new_data_rkaa_list[index]['anno_type']
                    if key[-6:] == '_genes':
                        list_key.append(key)

            # add a key-value pair contain all name in geneset
            for idx, each in enumerate(new_data_rkaa_list):
                list_all = list()
                for k in set(list_key):
                    list_all.extend(each[k])
                list_all = list(set(list_all))
                str_all = ';'.join(list_all)
                new_data_rkaa_list[idx]['list_all'] = list_all
                new_data_rkaa_list[idx]['str_all'] = str_all

            self.create_db_table('sg_geneset_kegg_class_detail', new_data_rkaa_list)
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
                            #genes = [i for i in genes if genes.count(i) == 1]
                            genes = list(set(genes))
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
            list_key = list()

            # important change: *_str = "GENE1,GENE2,..." --> *_str = "GENE1;GENE2;..."
            # important change: anno_type = "G" --> geneset_type = "G"
            for index, record in enumerate(new_data_rkaa_list):
                for key in record.keys():
                    if key[-4:] == '_str':
                        new_data_rkaa_list[index][key] = record[key].replace(',', ';')
                    if key == 'anno_type':
                        new_data_rkaa_list[index]['geneset_type'] = 'G'
                        del new_data_rkaa_list[index]['anno_type']
                    if key[-6:] == '_genes':
                        list_key.append(key)

            # add a key-value pair contain all name in geneset
            for idx, each in enumerate(new_data_rkaa_list):
                list_all = list()
                for k in set(list_key):
                    list_all.extend(each[k])
                list_all = list(set(list_all))
                str_all = ';'.join(list_all)
                new_data_rkaa_list[idx]['list_all'] = list_all
                new_data_rkaa_list[idx]['str_all'] = str_all

            self.create_db_table('sg_geneset_kegg_class_detail', new_data_rkaa_list)
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

    @report_check
    def add_kegg_regulate_pic(self, main_table_id, level_path, png_dir):
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
        data_list = list()
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                pid = re.sub('path:', '', line[0])
                if os.path.exists(png_dir + '/' + line[0] + '.html.mark'):
                    with open(png_dir + '/' + line[0] + '.html.mark', 'r') as mark_f:
                        for line_mark in mark_f.readlines():
                            # print len(line_mark.strip("\n").split("\t"))
                            if len(line_mark.strip("\n").split("\t")) == 10:
                                [png, shape, bg_color, fg_color, coords, title, kos, href, gene_list, geneset_list] = line_mark.strip("\n").split("\t")
                                title = title.replace("\\n", "\n")
                            else:
                                continue
                            ''' 无参的为 #00CD00 不修改
                            if bg_color == "#00CD00":
                                bg_color = "#FFFF00"
                            '''
                            insert_data = {
                                'kegg_id': kegg_id,
                                'pathway_id': line[0],
                                'shape': shape,
                                'bg_colors': bg_color,
                                'fg_colors': fg_color,
                                'coords': coords,
                                'href': href,
                                'kos': kos,
                                'title': title,
                                'gene_list': gene_list.split("|"),
                                'geneset_list': geneset_list.split("|")
                            }

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
            except Exception, e:
                raise Exception("导入kegg注释图片信息：%s、%s出错!" % (level_path, png_dir))
            else:
                self.bind_object.logger.info("导入kegg注释图片信息：%s、%s 成功!" % (level_path, png_dir))

class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''
    # TODO: test functions
    pass

if __name__ == '__main__':
    unittest.main()

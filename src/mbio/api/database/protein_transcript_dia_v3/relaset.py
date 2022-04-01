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
from mbio.api.database.dia.api_base import ApiBase
import math
import time
from collections import OrderedDict


class Relaset(ApiBase):
    def __init__(self, bind_object):
        super(Relaset, self).__init__(bind_object)
        self._project_type = 'dia'

    def add_relaset_cluster(self, cluster_output_dir, main_id=None, project_sn='dia', task_id='dia',
                            params=None):
        # prepare main_table data
        results = os.listdir(cluster_output_dir)
        gene_cluster, sample_cluster = False, False
        genes, samples = list(), list()
        gene_tree, sample_tree = "", ""
        if "seq.cluster_tree.txt" in results:
            target_file = os.path.join(cluster_output_dir, "seq.cluster_tree.txt")
            with open(target_file) as f:
                gene_cluster = True
                gene_tree = f.readline().strip()
                genes = f.readline().strip().split(";")
            exp = os.path.join(cluster_output_dir, "expression_matrix.xls")
            exp_matrix = pd.read_table(exp, header=0, index_col=None, sep='\t')
            seq_id = exp_matrix['accession_id'].tolist()
            for i in seq_id:
                if ':' in i:
                    if i.replace(':', '-') in genes:
                        genes[genes.index(i.replace(':', '-'))] = i
        #
        if "sample.cluster_tree.txt" in results:
            target_file = os.path.join(cluster_output_dir, "sample.cluster_tree.txt")
            with open(target_file) as f:
                sample_cluster = True
                sample_tree = f.readline().strip()
                samples = f.readline().strip().split(";")
        #
        if "seq.kmeans_cluster.txt" in results:
            gene_cluster = True
            target_file = os.path.join(cluster_output_dir, "seq.kmeans_cluster.txt")
            with open(target_file) as f:
                genes = list()
                for line in f:
                    if not line.strip():
                        continue
                    genes += line.strip().split('\t')[1].split(";")
        #
        detail_info = list()
        trend_dict = dict()
        if ("seq.cluster_tree.txt" in results) or ("seq.kmeans_cluster.txt" in results) :
            sub_clusters = [x for x in results if x.startswith('seq.subcluster')]
            number_order = [(x, int(x.split('_')[1])) for x in sub_clusters]
            tmp = sorted(number_order, key=lambda x: x[1])
            sub_clusters = [x[0] for x in tmp]
            for sub in sub_clusters:
                target_file = os.path.join(cluster_output_dir, sub)
                tmp_df = pd.read_table(target_file, header=0)
                sub_cluster_id = int(sub.split('_')[1])
                tmp_df["sub_cluster"] = sub_cluster_id
                detail_info += json.loads(tmp_df.to_json(orient="records"))
                mean_dict = tmp_df.iloc[:, 1:-1].mean().to_dict()
                trend_dict[str(sub_cluster_id)] = mean_dict
        #
        target_file = os.path.join(cluster_output_dir, "expression_matrix.xls")
        exp_pd = pd.read_table(target_file, header=0)
        if not detail_info:
            detail_info = exp_pd.to_dict('records')
        detail_info = [dict(d, **{'accession_id_protein':d['accession_id'].split('|')[0],'accession_id_gene':d['accession_id'].split('|')[1]}) for d in detail_info]
        if not genes:
            genes = list(exp_pd['accession_id'])
        if not samples:
            samples = list(exp_pd.columns)[1:]
        # add main table info'
        if main_id is None:
            name = "RelaSet_Cluster" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if params is None:
                params_dict = dict()
            elif type(params) == dict:
                params_dict = params
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            else:
                params_dict = json.loads(params)

            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='proteinset cluster main table',
                status="start",
                params=params
            )
            main_id = self.create_db_table('sg_relaset_cluster', [main_info])
        else:
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
        # update main table
        self.update_db_record('sg_relaset_cluster', main_id,
                              trend_dict=trend_dict,
                              samples=samples,
                              gene_cluster=gene_cluster,
                              sample_cluster=sample_cluster, )
        # add detail info
        tree_info = dict(
            relas=genes,
            rela_tree=gene_tree,
            sample_tree=sample_tree,
            cluster_id=main_id,
        )
        self.create_db_table('sg_relaset_cluster_tree', [tree_info])
        self.create_db_table('sg_relaset_cluster_detail', detail_info, tag_dict=dict(cluster_id=main_id))
        self.update_db_record('sg_relaset_cluster', main_id, status="end", main_id=main_id, )
        return main_id

    @report_check
    def add_relaset_corr(self, corr_output_dir, main_id=None, workflow_out=None):
        # prepare main_table data
        results = os.listdir(corr_output_dir)
        if 'protein_rna_density_matrix.xls' and 'corr_result.xls' in results:
            density_file = os.path.join(corr_output_dir, 'protein_rna_density_matrix.xls')
            corr_file = os.path.join(corr_output_dir, 'corr_result.xls')
        else:
            self.bind_object.set_error('找不到密度文件')
        result_pd = pd.read_table(density_file)
        result_pd = result_pd.fillna('_')
        result_pd['corr_id'] = ObjectId(main_id)
        result_dict_list = result_pd.to_dict("records")
        result_dict_list = [dict(d, **{'rela_id_protein':d['rela_id'].split('|')[0],'rela_id_gene':d['rela_id'].split('|')[1]}) for d in result_dict_list]
        muilty_list = list()
        with open(corr_file) as c:
            headers = c.readline().strip().split('\t')
            for line in c:
                corr_dict = OrderedDict()
                corr_dict['corr_id'] = ObjectId(main_id)
                line = line.strip().split('\t')
                corr_dict['sample'] = line[0]
                for n, value in enumerate(line[1:]):
                    corr_dict[headers[n]] = float(value)
                muilty_list.append(corr_dict)

        muilty_pdf = os.path.join(workflow_out, 'multily_corr.pdf')
        muilty_svg = os.path.join(workflow_out, 'multily_corr.svg')
        muilty_png = os.path.join(workflow_out, 'multily_corr.png')
        with open(os.path.join(corr_output_dir, 'corr_info')) as cor_r:
            corr_value = cor_r.readline().strip().split('\t')[1]
            p_value = cor_r.readline().strip().split('\t')[1]
        self.create_db_table("sg_relaset_corr_detail", result_dict_list)
        self.create_db_table("sg_relaset_corr_muilty", muilty_list)
        self.update_db_record('sg_relaset_corr', main_id, status="end", main_id=ObjectId(main_id), muilty_pdf= muilty_pdf, corr_value=corr_value, p_value=p_value, samples=headers, muilty_svg=muilty_svg, muilty_png=muilty_png)

    @report_check
    def add_pr_go_class(self, go_merge_file, main_id=None):
        result_pd = pd.read_table(go_merge_file)
        result_pd.columns = ['term_type', 'term', 'go_id', 'level', 'protein_num', 'protein_percent',
'protein_list', 'gene_num', 'gene_percent', 'gene_list']
        result_pd.loc[(result_pd['protein_list'].str.contains(";", na=False)), 'protein_str'] = result_pd['protein_list'].str.split(';')
        isna = result_pd['protein_str'].isna()
        result_pd.loc[isna, 'protein_str'] = pd.Series([[]] * isna.sum()).values
        result_pd['protein_str']=result_pd['protein_str'].fillna(pd.Series([]))
        result_pd = result_pd.fillna('_')
        result_pd['pr_go_id'] = ObjectId(main_id)
        result_dict_list = result_pd.to_dict("records")
        self.create_db_table("sg_pr_go_class_detail", result_dict_list)
        self.update_db_record('sg_pr_go_class', main_id, status="end", main_id=ObjectId(main_id))

    def add_kegg_regulate_pr(self, main_table_id, proteinset_file, kegg_stat_file, gene_kegg_level_table_xls, workflow_out):
        # kegg 分类导表
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                raise Exception('kegg_main_table_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(gene_kegg_level_table_xls):
            raise Exception('gene_kegg_level_table_xls所指定的路径:{}不存在，请检查！'.format(gene_kegg_level_table_xls))

        # 读入基因集列表
        with open(proteinset_file, 'r') as pset:
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

        def len_(x):
            while '' in x:
                x.remove('')
            if not x:
                return 0
            return len(x)
        for gene_set in pset_list:
            stat_class.rename(columns={gene_set + '_genes': gene_set + '_geneko'}, inplace=True)
            stat_class[gene_set + '_str'] = stat_class[gene_set + '_geneko'].replace(r'\([^\)]*\)', '', regex=True)
            stat_class[gene_set + '_genes'] = stat_class[gene_set + '_str'].map(lambda x: list(set(x.split(";"))))
            stat_class[gene_set + '_numbers'] = stat_class[gene_set + '_genes'].map(lambda x: len_(x))
            stat_class[gene_set + '_str'] = stat_class[gene_set + '_genes'].map(lambda x: ",".join(x))

        kegg_class_detail = stat_class.to_dict('records')
        self.create_db_table('sg_pr_kegg_class_detail', kegg_class_detail)
        self.update_db_record('sg_pr_kegg_class', ObjectId(main_table_id), main_id=ObjectId(main_table_id))

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
                class_genes = list(set(class_genes))
                while '' in class_genes:
                    class_genes.remove('')
                data.extend([
                    (gene_set + '_genes', class_genes),
                    (gene_set + '_genes_num', len(class_genes))
                ])
            data = SON(data)
            data_list.append(data)
        categories = ["".join(map(lambda y:y[0], x.split(' '))) for x in list(set(stat_class['first_category']))]
        try:
            collection = self.db['sg_pr_kegg_class_statistic']
            collection.insert_many(data_list)
            self.update_db_record('sg_pr_kegg_class', main_table_id, categories=categories, table_columns=pset_list, result_dir=workflow_out)
        except Exception as e:
            raise Exception("导入kegg注释分类信息：%s出错!" % (kegg_stat_file))
        else:
            self.bind_object.logger.info("导入kegg注释分类信息：%s 成功!" % kegg_stat_file)

    def add_kegg_regulate_pic(self, main_table_id, level_path, png_dir, bg_dict, gene_ko, species_protein, ko_txt):
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
        col_dict = dict(
            pink='#FFC0CB',
            green='#00CD00',
            orange='#FFA500'
        )
        ko2bg = dict()
        with open(ko_txt) as kr:
            _ = kr.readline()
            for line in kr:
                if line.strip():
                    line = line.strip().split('\t')
                    ko2bg[line[0]] = line[1]
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                pid = re.sub('path:', '', line[0])
                if os.path.exists(png_dir + '/' + line[0] + '.html.mark'):
                    with open(png_dir + '/' + line[0] + '.html.mark', 'r') as mark_f:
                        for line_mark in mark_f.readlines():
                            # print len(line_mark.strip("\n").split("\t"))
                            if len(line_mark.strip("\n").split("\t")) == 8:
                                [png, shape, bg_color, fg_color, coords, title, kos, href] = line_mark.strip("\n").split("\t")
                                title = title.replace("\\n", "\n")
                                png = png.split('.png')[0]
                                png = png.replace('map', species_protein)
                                kos_ = list()
                                for ko in kos.split(','):
                                    if ko in gene_ko:
                                        # kos_ += gene_ko[ko]
                                        kos_ += ko
                                kos_ = kos.split(',') + kos_
                                # self.bind_object.logger.info(kos_)
                                # if bg_color and png in bg_dict:
                                #     for ko in kos_:
                                #         if ko in bg_dict[png]:
                                #             if bg_dict[png][ko] == 'orange':
                                #                 bg_color = '#FEE090'
                                #                 break
                                #             else:
                                #                 bg_color = col_dict[bg_dict[png][ko]]
                                bgs = set()
                                if bg_color:
                                    for ko in kos_:

                                        if ko in ko2bg:
                                            # if ko2bg[ko] == '#FFA500':
                                            #     bg_color = '#FFA500'
                                            #     break
                                            # else:
                                            #     bg_color = ko2bg[ko]
                                            #     break
                                            #bg_color = ko2bg[ko]
                                            bgs.add(ko2bg[ko])
                                if not bg_color:
                                    kos_ = list()
                                    for ko in title.split(','):
                                        if u' (' in ko:
                                            ko = ko.split(' (')[0].strip()
                                            kos_.append(ko)
                                    for ko in kos_:
                                        # if ko in bg_dict[png]:
                                        #     if bg_dict[png][ko] == 'orange':
                                        #         bg_color = '#FEE090'
                                        #         break
                                        #     else:
                                        #         bg_color = col_dict[bg_dict[png][ko]]
                                        if ko in ko2bg:
                                            # if ko2bg[ko] == '#FFA500':
                                            #     bg_color = '#FFA500'
                                            #     break
                                            # else:
                                            #     bg_color = ko2bg[ko]
                                            #     break
                                            # bg_color = ko2bg[ko]
                                            bgs.add(ko2bg[ko])
                                if bgs:
                                    if len(bgs)>1:
                                        bg_color = '#FFA500'
                                    else:
                                        bg_color = list(bgs)[0]
                            else:
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

                            if bg_color != "" and len(bg_color.split(",")) > 0:
                                insert_data.update({'bg_type': len(bg_color.split(","))})
                            if fg_color != "" and len(fg_color.split(",")) > 0:
                                insert_data.update({'fg_type': len(fg_color.split(","))})
                            data_list.append(insert_data)
                else:
                    self.bind_object.logger.info("kegg 图片{} 不存在html标记!".format(line[0]))

        if data_list:
            try:
                collection = self.db['sg_pr_kegg_class_pic']
                collection.insert_many(data_list)
            except Exception as e:
                raise Exception("导入kegg注释图片信息：%s、%s出错!" % (level_path, png_dir))
            else:
                self.bind_object.logger.info("导入kegg注释图片信息：%s、%s 成功!" % (level_path, png_dir))



    @report_check
    def add_go_enrich_pr(self, go_enrich_id, go_enrich_protein, go_enrich_rna):
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
        if not os.path.exists(go_enrich_protein):
            raise Exception('{}所指定的路径不存在。请检查！'.format(go_enrich_protein))
        data_list = []
        go_type = []
        go_ids = list()
        with open(go_enrich_protein, 'r') as f:
            f.readline()
            for line in f:
                line = line.strip().split('\t')
                if line[2] == 'p':
                    continue
                data = [
                    ('go_enrich_id', go_enrich_id),
                    ('enrich_type', 'protein'),
                    ('go_id', line[0]),
                    ('go_type', line[1]),
                    ('enrichment', line[2]),
                    ('discription', line[3]),
                    ('ratio_in_study', line[4]),
                    ('ratio_in_pop', line[5]),
                    ('p_uncorrected', float(line[6])),
                    ('p_corrected', float(line[-1])),
                    ('enrich_factor', float(line[4].split("/")[0]) / float(line[5].split("/")[0])),
                    ('depth', int(line[7])),
                    ('study_count', int(line[4].split("/")[0])),
                    ('pop_count', int(line[5].split("/")[0])),
                    ('seq_list', line[-2]),
                    ('seq_str', line[-2].split(";"))
                ]
                if line[0] not in go_ids:
                    go_ids.append(line[0])
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
            log10x = [-math.log10(x) if x > 0 else pvalues_min / 10 for x in pvalues]
            for i in range(0, len(log10x)):
                data_list[i]['neg_log10p_corrected'] = log10x[i]

        if not os.path.exists(go_enrich_rna):
            raise Exception('{}所指定的路径不存在。请检查！'.format(go_enrich_rna))
        data_list_rna = []
        go_type_rna = []
        with open(go_enrich_rna, 'r') as f:
            f.readline()
            for line in f:
                line = line.strip().split('\t')
                if line[2] == 'p':
                    continue
                data = [
                    ('go_enrich_id', go_enrich_id),
                    ('enrich_type', 'gene'),
                    ('go_id', line[0]),
                    ('go_type', line[1]),
                    ('enrichment', line[2]),
                    ('discription', line[3]),
                    ('ratio_in_study', line[4]),
                    ('ratio_in_pop', line[5]),
                    ('p_uncorrected', float(line[6])),
                    ('p_corrected', float(line[-1])),
                    ('enrich_factor', float(line[4].split("/")[0]) / float(line[5].split("/")[0])),
                    ('depth', int(line[7])),
                    ('study_count', int(line[4].split("/")[0])),
                    ('pop_count', int(line[5].split("/")[0])),
                    ('seq_list', line[-2]),
                    ('seq_str', line[-2].split(";"))
                ]
                if line[0] not in go_ids:
                    go_ids.append(line[0])
                go_type_rna.append(line[1])
                data = SON(data)
                data_list_rna.append(data)
        if data_list_rna:
            # 插入-logpvalue -logpcorrected 值相关字段
            pvalues = [dict(son)['p_uncorrected'] for son in data_list_rna]
            if len([x for x in pvalues if x > 0]) > 0:
                pvalues_min = min([x for x in pvalues if x > 0]) / 10
            else:
                pvalues_min = 0.0001
            pvalues_min = - math.log10(pvalues_min)
            log10x = [-math.log10(x) if x > 0 else pvalues_min for x in pvalues]
            for i in range(0, len(log10x)):
                data_list_rna[i]['neg_log10p_uncorrected'] = log10x[i]

            pvalues = [dict(son)['p_corrected'] for son in data_list_rna]
            if len([x for x in pvalues if x > 0]) > 0:
                pvalues_min = min([x for x in pvalues if x > 0]) / 10
            else:
                pvalues_min = 0.0001
            pvalues_min = - math.log10(pvalues_min)
            log10x = [-math.log10(x) if x > 0 else pvalues_min / 10 for x in pvalues]
            for i in range(0, len(log10x)):
                data_list_rna[i]['neg_log10p_corrected'] = log10x[i]

            try:
                collection = self.db['sg_pr_go_enrich_detail']
                collection.insert_many(data_list + data_list_rna)
                coll = self.db['sg_pr_go_enrich']
                go_type = list(set(go_type + go_type_rna))
                coll.update({'_id': go_enrich_id}, {'$set': {'categories': go_type, 'main_id': go_enrich_id}})
                # main_collection = self.db['sg_proteinset_go_enrich']
                # main_collection.update({"_id": ObjectId(go_enrich_id)}, {"$set": {"status": "end"}})
            except Exception as e:
                print("导入go富集信息：%s出错:%s" % (go_enrich_protein, e))
            else:
                print("导入go富集信息：%s成功!" % (go_enrich_protein))
        else:
            raise Exception('GO富集没有结果')

        merge_list = list()
        for go in go_ids:
            tmp_dict={'go_enrich_id': go_enrich_id, 'go_id':go}
            for data in data_list:
                data = dict(data)
                if data['go_id'] == go:
                    tmp_dict['go_type'] = data['go_type']
                    tmp_dict['enrichment'] = data['enrichment']
                    tmp_dict['discription'] = data['discription']
                    for k in ['ratio_in_study','ratio_in_pop','p_uncorrected','p_corrected','enrich_factor','depth']:
                        tmp_dict['%s_protein'%k] = data[k]
                        tmp_dict['%s_gene'%k] = -1
                    for k in ['seq_list','seq_str']:
                        tmp_dict['%s_protein' % k] = data[k]
                        tmp_dict['%s_gene' % k] = ''
                    for k in ['study_count','pop_count']:
                        tmp_dict['%s_protein' % k] = data[k]
                        tmp_dict['%s_gene' % k] = 0

            for data in data_list_rna:
                data = dict(data)
                if data['go_id'] == go:
                    tmp_dict['go_type'] = data['go_type']
                    tmp_dict['enrichment'] = data['enrichment']
                    tmp_dict['discription'] = data['discription']
                    for k in ['ratio_in_study','ratio_in_pop','p_uncorrected','p_corrected','enrich_factor','depth']:
                        if '%s_protein'%k not in tmp_dict:
                            tmp_dict['%s_protein' % k] = -1
                        tmp_dict['%s_gene'%k] = data[k]
                    for k in ['seq_list','seq_str']:
                        if '%s_protein'%k not in tmp_dict:
                            tmp_dict['%s_protein' % k] = ''
                        tmp_dict['%s_gene'%k] = data[k]
                    for k in ['study_count','pop_count']:
                        if '%s_protein'%k not in tmp_dict:
                            tmp_dict['%s_protein' % k] = 0
                        tmp_dict['%s_gene'%k] = data[k]
            #add float type for search function in webpage
            def transStr2Float(a):
                if isinstance(a, str) and '/' in a:
                    b = float(a.split('/')[0])/float(a.split('/')[1])
                else:
                    b = a
                return b
            for ratio_str in ['ratio_in_study_protein','ratio_in_pop_protein','ratio_in_study_gene','ratio_in_pop_gene']:
                tmp_dict[ratio_str+'_float'] = transStr2Float(tmp_dict[ratio_str])
            #
            merge_list.append(tmp_dict)

        try:
            collection = self.db['sg_pr_go_enrich_merge']
            collection.insert_many(merge_list)
        except Exception as e:
            print("导入go富集信息：%s出错:%s" % (go_enrich_protein, e))
        else:
            print("导入go富集信息：%s成功!" % (go_enrich_protein))


    @report_check
    def add_kegg_enrich_pr(self, enrich_id, kegg_enrich_protein, kegg_enrich_rna):
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
        if not os.path.exists(kegg_enrich_protein):
            raise Exception('kegg_enrich_table所指定的路径:{}不存在，请检查！'.format(kegg_enrich_protein))
        data_list = []
        kegg_type1 = []
        kegg_list = list()
        try:
            species_prefix = open(kegg_enrich_protein,'r').readlines()[1].split('\t')[3][:-5]
        except:
            species_prefix = "map"
            pass
        with open(kegg_enrich_protein, 'rb') as r:
            for line in r:
                if re.match(r'\w', line):
                    line = line.strip('\n').split('\t')
                    insert_data = {
                        'kegg_enrich_id': enrich_id,
                        'enrich_type': 'protein',
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
                        'seq_str': line[8].split('|'),
                        'hyperlink': line[9],
                        'kegg_type': "".join([x[0] for x in line[11].split(' ')])
                    }
                    id = line[3].split("path:")[1] if "path:" in line[3] else line[3]
                    if not id in kegg_list:
                        kegg_list.append(id)
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

                # kegg_type1 = list(set(kegg_type1))
        if not os.path.exists(kegg_enrich_rna):
            raise Exception('kegg_enrich_table所指定的路径:{}不存在，请检查！'.format(kegg_enrich_rna))
        data_list1 = []
        kegg_type2 = []
        with open(kegg_enrich_rna, 'rb') as r:
            for line in r:
                if re.match(r'\w', line):
                    line = line.strip('\n').split('\t')
                    insert_data = {
                        'kegg_enrich_id': enrich_id,
                        'enrich_type': 'gene',
                        'term': line[1],
                        'database': line[2],
                        'id': line[3].split("path:")[1] if "path:" in line[3] else line[3],
                        'discription': line[1],
                        'study_count': int(line[0]),
                        "background_number": line[5].split("/")[1],
                        'ratio_in_study': line[4],
                        'ratio_in_pop': line[5],
                        'enrich_factor': float(line[0]) / float(line[5].split("/")[0]),
                        'pvalue': round(float(line[6]), 4),
                        'corrected_pvalue': round(float(line[7]), 4) if not line[7] == "None" else "None",
                        'seq_list': line[8],
                        'seq_str': line[8].split('|'),
                        'hyperlink': line[9],
                        'kegg_type': "".join([x[0] for x in line[11].split(' ')])
                    }
                    id = line[3].split("path:")[1] if "path:" in line[3] else line[3]
                    if not id in kegg_list:
                        kegg_list.append(id)
                    kegg_type2.append("".join([x[0] for x in line[11].split(' ')]))
                    data_list1.append(insert_data)
            if data_list1:
                # 插入-logpvalue -logpcorrected 值相关字段
                pvalues = [dict(son)['pvalue'] for son in data_list1]
                if len([x for x in pvalues if x > 0]) > 0:
                    pvalues_min = min([x for x in pvalues if x > 0]) / 10
                else:
                    pvalues_min = 0.0001
                pvalues_min = - math.log10(pvalues_min)
                log10x = [-math.log10(x) if x > 0 else pvalues_min for x in pvalues]
                for i in range(0, len(log10x)):
                    data_list1[i]['neg_log10p_uncorrected'] = log10x[i]

                pvalues = [dict(son)['corrected_pvalue'] for son in data_list1]
                if len([x for x in pvalues if x > 0]) > 0:
                    pvalues_min = min([x for x in pvalues if x > 0]) / 10
                else:
                    pvalues_min = 0.0001
                pvalues_min = - math.log10(pvalues_min)

                log10x = [-math.log10(x) if x > 0 else pvalues_min for x in pvalues]
                for i in range(0, len(log10x)):
                    data_list1[i]['neg_log10p_corrected'] = log10x[i]

                kegg_type1 = list(set(kegg_type1 + kegg_type2))
                try:
                    collection = self.db['sg_pr_kegg_enrich_detail']
                    collection.insert_many(data_list + data_list1)
                    coll = self.db['sg_pr_kegg_enrich']

                    coll.update({'_id': enrich_id}, {'$set': {'categories': kegg_type1}})
                    # main_collection = self.db['sg_proteinset_kegg_enrich']
                    # main_collection.update({"_id": ObjectId(enrich_id)}, {"$set": {"status": "end"}})
                except Exception as e:
                    self.bind_object.set_error("导入kegg富集统计表：%s信息出错:%s" % (kegg_enrich_protein, e))
                else:
                    self.bind_object.logger.info("导入kegg富集统计表:%s信息成功!" % kegg_enrich_protein)
            else:
                coll = self.db['sg_pr_kegg_enrich']
                coll.update({'_id': enrich_id}, {'$set': {'desc': 'no_result'}})
                # self.bind_object.logger.info("kegg富集统计表没结果：" % kegg_enrich_table)
                raise Exception("kegg富集统计表没结果")

        merge_list = list()
        # 为了整合蛋白单物种转录全物种的情况
        kegg_list = [filter(str.isdigit, kegg) for kegg in kegg_list]
        kegg_list = list(set(kegg_list))
        for kegg in kegg_list:
            tmp_dict = {'kegg_enrich_id': enrich_id, 'id': species_prefix + kegg}
            for data in data_list:
                data = dict(data)
                if filter(str.isdigit, data['id']) == kegg:
                    tmp_dict['term'] = data['term']
                    tmp_dict['database'] = data['database']
                    tmp_dict['discription'] = data['discription']
                    tmp_dict['kegg_type'] = data['kegg_type']
                    for k in ['ratio_in_study', 'ratio_in_pop', 'pvalue', 'corrected_pvalue', 'enrich_factor']:
                        tmp_dict['%s_protein' % k] = data[k]
                        tmp_dict['%s_gene' % k] = -1
                    for k in ['seq_list', 'hyperlink']:
                        tmp_dict['%s_protein' % k] = data[k]
                        tmp_dict['%s_gene' % k] = ''
                    for k in ['study_count', 'background_number']:
                        tmp_dict['%s_protein' % k] = data[k]
                        tmp_dict['%s_gene' % k] = 0

            for data in data_list1:
                data = dict(data)
                if filter(str.isdigit, data['id']) == kegg:
                    tmp_dict['term'] = data['term']
                    tmp_dict['database'] = data['database']
                    tmp_dict['discription'] = data['discription']
                    tmp_dict['kegg_type'] = data['kegg_type']
                    for k in ['ratio_in_study', 'ratio_in_pop', 'pvalue', 'corrected_pvalue', 'enrich_factor']:
                        if '%s_protein' % k not in tmp_dict:
                            tmp_dict['%s_protein' % k] = -1
                        tmp_dict['%s_gene' % k] = data[k]
                    for k in ['seq_list', 'hyperlink']:
                        if '%s_protein' % k not in tmp_dict:
                            tmp_dict['%s_protein' % k] = ''
                        tmp_dict['%s_gene' % k] = data[k]
                    for k in ['study_count', 'background_number']:
                        if '%s_protein' % k not in tmp_dict:
                            tmp_dict['%s_protein' % k] = 0
                        tmp_dict['%s_gene' % k] = data[k]
            #add float type for search function in webpage
            def transStr2Float(a):
                if isinstance(a, str) and '/' in a:
                    b = float(a.split('/')[0])/float(a.split('/')[1])
                else:
                    b = a
                return b
            for ratio_str in ['ratio_in_study_protein','ratio_in_pop_protein','ratio_in_study_gene','ratio_in_pop_gene']:
                tmp_dict[ratio_str+'_float'] = transStr2Float(tmp_dict[ratio_str])
            #
            merge_list.append(tmp_dict)

        try:
            collection = self.db['sg_pr_kegg_enrich_merge']
            collection.insert_many(merge_list)
        except Exception as e:
            print("导入kegg富集信息：%s出错:%s" % (kegg_enrich_protein, e))
        else:
            print("导入kegg富集信息：%s成功!" % (kegg_enrich_protein))


    def add_enrich_cluster(self, cluster_output_dir, main_id=None, project_sn='dia', task_id='dia',
                            params=None, type_ = 'go'):
        # prepare main_table data
        results = os.listdir(cluster_output_dir)
        gene_cluster, sample_cluster = False, False
        genes, samples = list(), list()
        gene_tree, sample_tree = "", ""
        if "seq.cluster_tree.txt" in results:
            target_file = os.path.join(cluster_output_dir, "seq.cluster_tree.txt")
            with open(target_file) as f:
                gene_cluster = True
                gene_tree = f.readline().strip()
                genes = f.readline().strip().split(";")
        #
        if "sample.cluster_tree.txt" in results:
            target_file = os.path.join(cluster_output_dir, "sample.cluster_tree.txt")
            with open(target_file) as f:
                sample_cluster = True
                sample_tree = f.readline().strip()
                samples = f.readline().strip().split(";")
        #
        if "seq.kmeans_cluster.txt" in results:
            gene_cluster = True
            target_file = os.path.join(cluster_output_dir, "seq.kmeans_cluster.txt")
            with open(target_file) as f:
                genes = list()
                for line in f:
                    if not line.strip():
                        continue
                    genes += line.strip().split('\t')[1].split(";")
        #
        detail_info = list()
        trend_dict = dict()
        if ("seq.cluster_tree.txt" in results) or ("seq.kmeans_cluster.txt" in results) :
            sub_clusters = [x for x in results if x.startswith('seq.subcluster')]
            number_order = [(x, int(x.split('_')[1])) for x in sub_clusters]
            tmp = sorted(number_order, key=lambda x: x[1])
            sub_clusters = [x[0] for x in tmp]
            for sub in sub_clusters:
                target_file = os.path.join(cluster_output_dir, sub)
                tmp_df = pd.read_table(target_file, header=0)
                sub_cluster_id = int(sub.split('_')[1])
                tmp_df["sub_cluster"] = sub_cluster_id
                # detail_info += json.loads(tmp_df.to_json(orient="records"))
                detail_info += tmp_df.to_dict('records')
                mean_dict = tmp_df.iloc[:, 1:-1].mean().to_dict()
                trend_dict[str(sub_cluster_id)] = mean_dict
        #
        target_file = os.path.join(cluster_output_dir, "cluster_matrix.xls")
        exp_pd = pd.read_table(target_file, header=0)
        detail_info_tmp = exp_pd.to_dict('records')
        for i in detail_info_tmp:
            for j in detail_info:
                if i[type_] == j[type_]:
                    j.update(**i)
                    break
        if not detail_info:
            detail_info = exp_pd.to_dict('records')
        if not genes:
            genes = list(exp_pd['accession_id'])
        if not samples:
            try:
                samples = list(tmp_df.columns)[1:-1]
            except:
                samples = list(exp_pd.columns)[1:]
        for i in ['study_num', 'database', 'term']:
            if i in samples:
                samples.remove(i)
        # add main table info'
        if main_id is None:
            name = "Enrich_Cluster" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if params is None:
                params_dict = dict()
            elif type(params) == dict:
                params_dict = params
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            else:
                params_dict = json.loads(params)

            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='proteinset cluster main table',
                status="start",
                params=params
            )
            main_id = self.create_db_table('sg_pr_%s_enrich_cluster'%type_, [main_info])
        else:
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
        # update main table
        self.update_db_record('sg_pr_%s_enrich_cluster'%type_, main_id,
                              trend_dict=trend_dict,
                              samples=samples,
                              gene_cluster=gene_cluster,
                              sample_cluster=sample_cluster, )
        # add detail info
        if type_ == 'go':
            genes = ['GO:' + go.lstrip('GO-') for go in genes]
            gene_tree = gene_tree.replace('GO', 'GO:').replace('-', '')
        tree_info = dict(
            ids=genes,
            id_tree=gene_tree,
            sample_tree=sample_tree,
            cluster_id=main_id,
        )
        self.create_db_table('sg_pr_%s_enrich_cluster_tree'%type_, [tree_info])
        self.create_db_table('sg_pr_%s_enrich_cluster_detail'%type_, detail_info, tag_dict=dict(cluster_id=main_id))
        self.update_db_record('sg_pr_%s_enrich_cluster'%type_, main_id, status="end", main_id=main_id, )
        return main_id

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
    def add_proteinset_cog_detail(self, proteinset_cog_table, proteinset_cog_id):
        """
        cog详情表导表函数
        :param proteinset_cog_table:cog结果表
        :param proteinset_cog_id:主表ID
        :return:
        """
        data_list = []
        proteinset_name = []
        with open(proteinset_cog_table, 'r') as f:
            first_line = f.readline().strip("\n").split("\t")
            for gn in first_line[2:]:
                if "list" in gn or "LIST" in gn:
                    continue
                elif not gn[:-4] in proteinset_name:
                    proteinset_name.append(gn[:-4])
            self.bind_object.set_error(proteinset_name)
            for line in f:
                line = line.strip().split("\t")
                data = {
                    'proteinset_cog_id': ObjectId(proteinset_cog_id),
                    'type': line[0],
                    'function_categories': line[1]
                }
                for n, gn in enumerate(proteinset_name):
                    self.bind_object.logger.info(n)
                    self.bind_object.logger.info(gn)
                    self.bind_object.logger.info(proteinset_name)
                    # data[gn + "_cog"] = int(line[6*n+2])
                    # data[gn + "_nog"] = int(line[6*n+3])
                    data[gn + "_cog"] = int(line[4*n+2])
                    data[gn + "_nog"] = int(line[4*n+3])
                    # data[gn + "_kog"] = int(line[6*n+4])
                    if data[gn + "_cog"] == 0:
                        data[gn + "_cog_list"] = ""
                        data[gn + "_cog_str"] = ""
                    else:
                        # data[gn + "_cog_list"] = line[6*n+5].split(";")
                        # data[gn + "_cog_str"] = line[6*n+5]
                        data[gn + "_cog_list"] = line[4*n+4].split(";")
                        data[gn + "_cog_str"] = line[4*n+4]
                    if data[gn + "_nog"] == 0:
                        data[gn + "_nog_str"] = ""
                        data[gn + "_nog_list"] = ""
                    else:
                        # data[gn + "_nog_str"] = line[6*n+6]
                        # data[gn + "_nog_list"] = line[6*n+6].split(";")
                        data[gn + "_nog_str"] = line[4*n+5]
                        data[gn + "_nog_list"] = line[4*n+5].split(";")
                    # if data[gn + "_kog"] == 0:
                    #     data[gn + "_kog_list"] = ""
                    #     data[gn + "_kog_str"] = ""
                    # else:
                    #     data[gn + "_kog_list"] = line[6*n+7].split(";")
                    #     data[gn + "_kog_str"] = line[6*n+7]
                data_list.append(data)
        try:
            collection = self.db['sg_proteinset_cog_class_detail']
            main_collection = self.db['sg_proteinset_cog_class']
            collection.insert_many(data_list)
            main_collection.update({"_id": ObjectId(proteinset_cog_id)},
                                   {"$set": {"table_columns": proteinset_name, "status": "end"}})
            self.bind_object.logger.info(proteinset_name)
        except Exception as e:
            self.bind_object.set_error("导入cog表格：%s出错:%s" % (proteinset_cog_table, e))
        else:
            self.bind_object.logger.info("导入cog表格：%s成功!" % (proteinset_cog_table))

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
                    ('p_corrected', float(line[-1])),
                    ('enrich_factor', float(line[4].split("/")[0])/float(line[5].split("/")[0])),
                    ('depth', int(line[7])),
                    ('study_count', int(line[4].split("/")[0])),
                    ('pop_count', int(line[5].split("/")[0])),
                    ('seq_list', line[-2]),
                    ('seq_str', line[-2].split(";"))
                ]
                go_type.append(line[1])
                data = SON(data)
                data_list.append(data)
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
                collection = self.db['sg_proteinset_go_enrich_detail']
                collection.insert_many(data_list)
                coll = self.db['sg_proteinset_go_enrich']
                go_type = list(set(go_type))
                coll.update({'_id': go_enrich_id}, {'$set': {'categories': go_type}})
                # main_collection = self.db['sg_proteinset_go_enrich']
                # main_collection.update({"_id": ObjectId(go_enrich_id)}, {"$set": {"status": "end"}})
            except Exception as e:
                print("导入go富集信息：%s出错:%s" % (go_enrich_dir, e))
            else:
                print("导入go富集信息：%s成功!" % (go_enrich_dir))
        else:
            raise Exception('GO富集没有结果')

    @report_check
    def update_directed_graph(self, go_enrich_id, go_graph_png, go_graph_pdf):
        collection = self.db['sg_proteinset_go_enrich']
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
        # proteinset_length = len(open(proteinset_list_path, "r").readlines())
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
                    collection = self.db['sg_proteinset_kegg_enrich_detail']
                    collection.insert_many(data_list)
                    coll = self.db['sg_proteinset_kegg_enrich']

                    coll.update({'_id': enrich_id}, {'$set': {'categories': kegg_type1}})
                    # main_collection = self.db['sg_proteinset_kegg_enrich']
                    # main_collection.update({"_id": ObjectId(enrich_id)}, {"$set": {"status": "end"}})
                except Exception as e:
                    self.bind_object.set_error("导入kegg富集统计表：%s信息出错:%s" % (kegg_enrich_table, e))
                else:
                    self.bind_object.logger.info("导入kegg富集统计表:%s信息成功!" % kegg_enrich_table)
            else:
                coll = self.db['sg_proteinset_kegg_enrich']
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
            proteinset_name = set(doc_keys)
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
                        data["{}_proteins".format(dk)] = line[5+n*3].split(";")
                    except:
                        data["{}_str".format(dk)] = ""
                        data["{}_proteins".format(dk)] = ""
                    if len(line4) > 1:
                        data["{}_percent_str".format(dk)] = line4[1][:-1]
                    else:
                        data["{}_percent_str".format(dk)] = 0
                data_list.append(data)
        try:
            collection = self.db['sg_proteinset_go_class_detail']
            main_collection = self.db['sg_proteinset_go_class']
            collection.insert_many(data_list)
            main_collection.update({"_id": ObjectId(go_regulate_id)}, {"$set": {"table_columns": list(proteinset_name)}})
            self.bind_object.logger.info(proteinset_name)
            self.bind_object.logger.info(ObjectId(go_regulate_id))
        except Exception as e:
            self.bind_object.logger.info("导入go调控信息：%s出错:%s" % (go_regulate_dir, e))
        else:
            self.bind_object.logger.info("导入go调控信息：%s成功!" % (go_regulate_dir))

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
            proteinset_name = set(doc_keys)
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
                        data["{}_proteins".format(dk)] = line[5+n*3].split(";")
                    except:
                        data["{}_str".format(dk)] = ""
                        data["{}_proteins".format(dk)] = ""
                    if len(line4) > 1:
                        data["{}_percent_str".format(dk)] = line4[1][:-1]
                    else:
                        data["{}_percent_str".format(dk)] = 0
                data_list.append(data)
        try:
            time.sleep(30)
            collection = self.db['sg_proteinset_go_class2_detail']
            main_collection = self.db['sg_proteinset_go_class2']
            collection.insert_many(data_list)
            main_collection.update({"_id": ObjectId(go_regulate_id)}, {"$set": {"table_columns": list(proteinset_name)}})
            self.bind_object.logger.info(proteinset_name)
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
            collection = self.db['sg_proteinset_kegg_class_pathway']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.info("导入kegg调控pathway：%s信息出错:%s" % (pathway_dir, e))
        else:
            self.bind_object.logger.info("导入kegg调控pathway:%s信息成功!" % pathway_dir)

    # 这个函数来自api_base.py，用在创建集合，也可以通过以前的传统方式
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
                print("Failed to insert records into table {} as: {}".format(table_name, e))
            else:
                print("Failed to insert record into table {} as: {}".format(table_name, e))
        else:
            if record_num >= 2:
                print("Success to insert records into table {}".format(table_name))
            else:
                print("Success to insert record into table {}".format(table_name))

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
    # 最后通过更新插入kegg_class主表的proteinset的名字;gene_kegg_level_table_xls是交互workflowkegg_table_2对应的值
    # kegg_stat_xls是kegg_class这个tool产生的，也是通过这个文件来更新sg_proteinset_kegg_class这个主表的table_columns字段
    def add_ipath_detail(self, main_table_id, ipath_input, proteinset):
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(ipath_input):
            raise Exception('{}所指定的路径不存在，请检查！'.format(ipath_input))

        ipath = pd.read_table(ipath_input, header=0)
        ipath.columns = ['accession_id', 'ko', 'color', 'width']
        ipath['ipath_id'] = main_table_id
        row_dict_list = ipath.to_dict('records')
        main_collection = self.db['sg_proteinset_ipath']

        doc_keys = []
        with open(proteinset, 'r') as f:
            for l in f.readlines():
                name = l.strip().split('\t')[0]
                if name not in doc_keys:
                    doc_keys.append(name)
        proteinset_name = list(set(doc_keys))


        try:
            self.create_db_table('sg_proteinset_ipath_detail', row_dict_list)
            self.bind_object.logger.info("主表id：{} 蛋白集：{}".format(main_table_id, proteinset_name))
            main_collection.update({"_id": main_table_id},
                                   {"$set": {"table_columns": proteinset_name, "status": "end"}})
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
            raise Exception('{}所指定的路径不存在，请检查！'.format(ipath_input))

        circ = pd.read_table(circ_input, header=0)
        circ.columns = ['accession_id', 'term', 'log2fc']
        circ['circ_id'] = main_table_id
        row_dict_list = circ.to_dict('records')
        main_collection = self.db['sg_proteinset_circ']

        try:
            self.create_db_table('sg_proteinset_circ_detail', row_dict_list)
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
        columns[0]='acession_id'
        columns[-1]='log2fc'

        circ.columns = columns
        circ['circ_id'] = main_table_id
        row_dict_list = circ.to_dict('records')
        main_collection = self.db['sg_proteinset_circ']

        self.bind_object.logger.info("导入circ：%s出错!" % (row_dict_list))
        try:
            self.create_db_table('sg_proteinset_circ_graph', row_dict_list)
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
        main_collection = self.db['sg_proteinset_circ']
        try:
            main_collection.update({"_id": main_table_id},
                                   {"$set": {"terms": term , "status": "end"}})
        except Exception as e:

            raise Exception("更新circ主表：%s出错!" % (circ_zscore_input))
        else:
            self.bind_object.logger.info("导入ipath：%s出错!" % (main_table_id))

    def add_kegg_regulate_new2(self, main_table_id, proteinset_file, kegg_stat_file, gene_kegg_level_table_xls):
        # kegg 分类导表
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                raise Exception('kegg_main_table_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(gene_kegg_level_table_xls):
            raise Exception('gene_kegg_level_table_xls所指定的路径:{}不存在，请检查！'.format(gene_kegg_level_table_xls))

        # 读入基因集列表
        with open(proteinset_file, 'r') as pset:
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

        for gene_set in pset_list:
            stat_class.rename(columns={ gene_set + '_genes' : gene_set + '_geneko'}, inplace=True)
            stat_class[gene_set + '_str']=stat_class[gene_set + '_geneko'].replace(r'\([^\)]*\)','',regex=True)
            stat_class[gene_set + '_genes']=stat_class[gene_set + '_str'].map(lambda x: x.split(";"))
            stat_class[gene_set + '_str']=stat_class[gene_set + '_genes'].map(lambda x: ",".join(x))

        kegg_class_detail = stat_class.to_dict('records')
        self.create_db_table('sg_proteinset_kegg_class_detail', kegg_class_detail)
        self.update_db_record('sg_proteinset_kegg_class', ObjectId(main_table_id), main_id=ObjectId(main_table_id))

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
                class_genes = list(set(class_genes))
                data.extend([
                    (gene_set + '_genes', class_genes),
                    (gene_set + '_genes_num', len(class_genes))
                ])
            data = SON(data)
            data_list.append(data)
        categories = ["".join(map(lambda y:y[0], x.split(' '))) for x in list(set(stat_class['first_category']))]
        try:
            collection = self.db['sg_proteinset_kegg_class_statistic']
            collection.insert_many(data_list)
            self.update_db_record('sg_proteinset_kegg_class', main_table_id, categories=categories, table_columns=pset_list)
        except Exception as e:
            raise Exception("导入kegg注释分类信息：%s出错!" % (kegg_stat_file))
        else:
            self.bind_object.logger.info("导入kegg注释分类信息：%s 成功!" % kegg_stat_file)

    def add_kegg_regulate_new(self, main_table_id, proteinset_id, kegg_stat_xls, gene_kegg_level_table_xls, work_dir):
        # 通过判断传入的proteinset_id的个数来确认取数据的位置，确认是一个还是两个基因集，然后现在分情况讨论
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

        if len(proteinset_id.split(",")) == 1:
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
                proteinset_name_r_kaa = head_r_kaa[2].split("numbers")[0].rstrip("_")
                str_name = proteinset_name_r_kaa + "_str"
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


            self.create_db_table('sg_proteinset_kegg_class_detail', new_data_rkaa_list)
            # self.create_db_table('sg_proteinset_kegg_class_detail', kaa_list)
            self.update_db_record('sg_proteinset_kegg_class', ObjectId(main_table_id), main_id=ObjectId(main_table_id))

            new_data = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep = '\t', header=0)
            new_data.groupby("second_category")
            group_obj = new_data.groupby("second_category")
            groups = group_obj.groups.keys()
            # 做一个基因集
            proteinsets = new_data.columns[3].split()
            result = defaultdict(dict)
            for each in groups:
                first = new_data.loc[new_data["second_category"] ==each]['first_category']
                first = first.to_dict().values()[0]
                for proteinset in proteinsets:
                    group_detail = group_obj.get_group(each)
                    genes = list()
                    for g in group_detail[proteinset]:

                        if not pd.isnull(g):
                            tmp = g.split(');')
                            genes += [x.split('(')[0] for x in tmp]
                            #genes = [i for i in genes if genes.count(i) == 1]
                            genes = list(set(genes))
                        else:
                            genes = []
                    result[proteinset][each] = [len(genes), first]

            try:
                a = pd.DataFrame(result)
                a.reset_index(inplace=True)
                a.rename(columns={a.columns[0]: "second_category"}, inplace=True)
                a.to_csv(work_dir + "/" + "k", sep = '\t', index=False)
                with open(work_dir + "/" + "k") as f1, open(work_dir + "/" + "kegg_statistic", "w") as fw:
                    header = f1.readline()
                    proteinset_name1 = header.strip().split("\t")[1]
                    proteinset_name1_num = proteinset_name1 + "_num"
                    fw.write("first_category" + "\t" + header.strip().split("\t")[0] + "\t" + proteinset_name1_num + "\n")
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
                appended_data_new_1['proteinset_id'] = ObjectId(proteinset_id)
                data_new = appended_data_new_1.to_dict('records')
                appended_data_new_1.to_csv(work_dir + "/" +  "kegg_statistic", sep='\t', index=False)
                # data_new = a.to_dict('records')
                self.create_db_table('sg_proteinset_kegg_class_statistic', data_new)

                with open(kegg_stat_xls, 'rb') as r:
                    # 获取numbers和proteinsets的列
                    first_line = r.readline().strip().split("\t")[2:]
                    # print r.next()
                    proteinsets_name = []
                    for fl in first_line:
                        if "numbers" in fl:
                            # 获取proteinset的name，
                            proteinsets_name.append(fl[:-8])

                main_collection = self.db['sg_proteinset_kegg_class']
                main_collection.update({"_id": ObjectId(main_table_id)}, {"$set": {"table_columns": proteinsets_name}})
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

        if len(proteinset_id.split(",")) == 2:
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
                proteinset_name_r_kaa_1 = head_r_kaa[2].split("numbers")[0].rstrip("_")
                str_name_1 = proteinset_name_r_kaa_1 + "_str"

                proteinset_name_r_kaa_2 = head_r_kaa[5].split("numbers")[0].rstrip("_")
                str_name_2 = proteinset_name_r_kaa_2 + "_str"
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
            self.create_db_table('sg_proteinset_kegg_class_detail', new_data_rkaa_list)
            # self.create_db_table('sg_proteinset_kegg_class_detail', kaa_list)
            self.update_db_record('sg_proteinset_kegg_class', ObjectId(main_table_id), main_id=ObjectId(main_table_id))

            new_data = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep = '\t', header=0)
            self.bind_object.logger.info("开始进行class分类导表")
            new_data.groupby("second_category")
            group_obj = new_data.groupby("second_category")
            groups = group_obj.groups.keys()
            # 做2个基因集
            proteinsets = new_data.columns[3], new_data.columns[5]
            result = defaultdict(dict)
            for each in groups:
                first = new_data.loc[new_data["second_category"] == each]['first_category']
                first = first.to_dict().values()[0]
                for proteinset in proteinsets:
                    group_detail = group_obj.get_group(each)
                    genes = list()
                    for g in group_detail[proteinset]:

                        # isnull支持的数据类型更多，相比isnan
                        if not pd.isnull(g):
                            tmp = g.split(');')
                            genes += [x.split('(')[0] for x in tmp]
                            # 用set会弹出不知名的错误
                            genes = [i for i in genes if genes.count(i) == 1]
                        else:
                            genes = []
                    result[proteinset][each] = [len(genes), first]
            # try:
            a = pd.DataFrame(result)
            a.reset_index(inplace=True)
            a.rename(columns={a.columns[0]: "second_category"}, inplace=True)
            a.to_csv(work_dir + "/" + "k", sep = '\t', index=False)
            with open(work_dir + "/" + "k") as f1, open(work_dir + "/" + "kegg_statistic", "w") as fw:
                header = f1.readline()
                proteinset_name1 = header.strip().split("\t")[1]
                proteinset_name1_num = proteinset_name1 + "_num"
                proteinset_name2 = header.strip().split("\t")[2]
                proteinset_name2_num = proteinset_name2 + "_num"
                fw.write("first_category" + "\t" + header.strip().split("\t")[0] + "\t" + proteinset_name1_num + "\t" + \
                         proteinset_name2_num +"\n")
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
            # appended_data_new_2['proteinset_type'] = proteinset_type
            # appended_data_new_2['proteinset_id'] = ObjectId(proteinset_id)
            appended_data_new_2['proteinset_id'] = proteinset_id


            appended_data_new_2.to_csv(work_dir + "/" +
                                       "kegg_statistic",sep='\t', index=False)
            data_new = appended_data_new_2.to_dict('records')
            # data_new = a.to_dict('records')
            self.create_db_table('sg_proteinset_kegg_class_statistic', data_new)
            self.bind_object.logger.info("完成class分类导表")

            with open(kegg_stat_xls, 'rb') as r:
                self.bind_object.logger.info("开始kegg主表的基因集名字信息更新")
                # 获取numbers和proteinsets的列
                first_line = r.readline().strip().split("\t")[2:]
                # print r.next()
                proteinsets_name = []
                for fl in first_line:
                    if "numbers" in fl:
                        # 获取proteinset的name，
                        proteinsets_name.append(fl[:-8])

            main_collection = self.db['sg_proteinset_kegg_class']
            main_collection.update({"_id": ObjectId(main_table_id)}, {"$set": {"table_columns": proteinsets_name}})
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
    # 最后通过更新插入kegg_class主表的proteinset的名字
    def add_kegg_regulate_detail(self, regulate_id, kegg_regulate_table):
        """
        :param regulate_id: 主表ID
        :param kegg_regulate_table: kegg_stat.xls统计结果文件
        :return:
        """
        # 会导表3张
        main_collection = self.db['sg_proteinset_kegg_class']
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
            # 获取numbers和proteinsets的列
            first_line = r.readline().strip().split("\t")[2:]
            # print r.next()
            proteinsets_name = []
            for fl in first_line:
                if "numbers" in fl:
                    # 获取proteinset的name，
                    proteinsets_name.append(fl[:-8])
            for line in r:
                line = line.strip('\n').split('\t')
                # regulate_id是表 "sg_proteinset_kegg_class"的_id
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
                for n, gn in enumerate(proteinsets_name):
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
                collection = self.db['sg_proteinset_kegg_class_detail']
                # main_collection = self.db['sg_proteinset_kegg_class']
                collection.insert_many(data_list)
                main_collection.update({"_id": ObjectId(regulate_id)}, {"$set": {"table_columns": proteinsets_name}})
            except Exception as e:
                self.bind_object.logger.info("导入kegg调控统计表：%s信息出错:%s" % (kegg_regulate_table, e))
            else:
                self.bind_object.logger.info("导入kegg调控统计表:%s信息成功!" % kegg_regulate_table)

    def add_proteinset(self, diff_work_dir, group_id=None, task_id=None, project_sn=None):
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
                self.bind_object.logger.info("开始导入{}的蛋白集".format(f))
                compare_all, protein_set_all = self.get_diff_list(diffcmp=diff_work_dir + "/" + f, up_down='all')
                compare_up, protein_set_up = self.get_diff_list(diffcmp=diff_work_dir + "/" + f, up_down='up')
                compare_down, protein_set_down = self.get_diff_list(diffcmp=diff_work_dir + "/" + f, up_down='down')
                if len(protein_set_all) >=1:
                    data_all = {
                        'group_id': group_id,
                        'task_id': task_id,
                        'name': compare_all,
                        'desc': '{}蛋白集'.format(compare_all),
                        'project_sn': project_sn,
                        'proteinset_length': len(protein_set_all),
                        'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                        'is_use': 1,
                    }
                    try:
                        collection = self.db["sg_proteinset"]
                        proteinset_up_id_all = collection.insert_one(data_all).inserted_id
                        self.update_db_record('sg_proteinset', ObjectId(proteinset_up_id_all), main_id=ObjectId(proteinset_up_id_all))
                        if proteinset_up_id_all:
                            self.add_proteinset_detail(proteinset_up_id_all, protein_set_all)
                    except Exception as e:
                        self.bind_object.set_error("导入蛋白集集：%s信息出错:%s" % (f, e))
                    else:
                        self.bind_object.logger.info("导入蛋白集：%s信息成功!" % f)
                if len(protein_set_up) >=1:
                    data_up = {
                        'group_id': group_id,
                        'task_id': task_id,
                        'name': compare_up,
                        'desc': '{}蛋白集'.format(compare_up),
                        'project_sn': project_sn,
                        'proteinset_length': len(protein_set_up),
                        'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                        'is_use': 1,
                    }
                    try:
                        collection = self.db["sg_proteinset"]
                        proteinset_up_id_up = collection.insert_one(data_up).inserted_id
                        self.update_db_record('sg_proteinset', ObjectId(proteinset_up_id_up), main_id=ObjectId(proteinset_up_id_up))
                        if proteinset_up_id_up:
                            self.add_proteinset_detail(proteinset_up_id_up, protein_set_up)
                    except Exception as e:
                        self.bind_object.set_error("导入蛋白集集：%s信息出错:%s" % (f, e))
                    else:
                        self.bind_object.logger.info("导入蛋白集：%s信息成功!" % f)
                if len(protein_set_down) >=1:
                    data_down = {
                        'group_id': group_id,
                        'task_id': task_id,
                        'name': compare_down,
                        'desc': '{}蛋白集'.format(compare_down),
                        'project_sn': project_sn,
                        'proteinset_length': len(protein_set_down),
                        'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                        'is_use': 1,
                    }
                    try:
                        collection = self.db["sg_proteinset"]
                        proteinset_up_id_down = collection.insert_one(data_down).inserted_id
                        self.update_db_record('sg_proteinset', ObjectId(proteinset_up_id_down), main_id=ObjectId(proteinset_up_id_down))
                        if proteinset_up_id_down:
                            self.add_proteinset_detail(proteinset_up_id_down, protein_set_down)
                    except Exception as e:
                        self.bind_object.set_error("导入蛋白集集：%s信息出错:%s" % (f, e))
                    else:
                        self.bind_object.logger.info("导入蛋白集：%s信息成功!" % f)
            else:
                continue

    def add_proteinset_detail(self, proteinset_id, sequence):
        if not isinstance(proteinset_id, ObjectId):
            if isinstance(proteinset_id, types.StringTypes):
                proteinset_id = ObjectId(proteinset_id)
            else:
                raise Exception('geneset_id必须为ObjectId对象或其对应的字符串！')
        data = [
            ("proteinset_id", ObjectId(proteinset_id)),
            ("seq_list", sequence)
        ]
        data = SON(data)
        try:
            collection = self.db["sg_proteinset_detail"]
            collection.insert_one(data)
        except Exception as e:
            self.bind_object.set_error("导入蛋白集detail表出错:%s" % e)
        else:
            self.bind_object.logger.info("导入蛋白集detail表成功!")

    def get_diff_list(self, diffcmp, up_down=None):
        if not os.path.exists(diffcmp):
            raise Exception("{}文件不存在，无法对up和down差异蛋白进行分类！".format(diffcmp))
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




if __name__ == "__main__":
    pass

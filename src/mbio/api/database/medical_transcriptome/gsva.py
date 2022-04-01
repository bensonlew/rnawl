# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
# last_modify:20200812
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import datetime
import json
import os
import pandas as pd
import glob
import re
from mbio.api.database.medical_transcriptome.api_base import ApiBase

class Gsva(ApiBase):
    def __init__(self, bind_object):
        super(Gsva, self).__init__(bind_object)

    def add_gsva(self, main_id, es_matrix, gmx, gsva, id2name, geneset_description, params=None, project_sn='medical_transcriptome', task_id='medical_transcriptome'):
        """
        dump exp data into database
        :param exp_matrix: expression matrix path or an express matrix in pandas DataFrame format.
        :param exp_level: str, transcript or gene
        :param exp_type: str, usually is tpm or fpkm or count. default tpm
        :param group_dict: ordered dict of group info
        :param quant_method: the method to be used to quant expression.
        :param project_sn: project id
        :param task_id: task id
        :param main_id: 主表id，如果提供，则本函数不创建主表
        :param group_id: 包含分组方案信息的主表id
        :param add_distribution: 是否添加表达分布信息到数据库
        :param params: parameters dict for expression quant.
        :return: main table id
        :version:version info
        """
        # add main table info
        if main_id is None:
            name = "Gsva" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Gsva',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_geneset_gsva', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        name_geneid_dict = dict()
        with open(id2name, 'r') as i:
            for line in i.readlines():
                try:
                    gene_id, name = line.strip().split('\t')
                except:
                    continue
                name_geneid_dict[name] = gene_id

        df = pd.read_table(es_matrix, header=0, sep='\t')
        gsva = pd.read_table(gsva, header=0, sep='\t')
        geneid = gsva['geneid']
        gene_name_dict = dict()
        with open(gmx, 'r') as g:
            for line in g.readlines():
                gene_name= line.strip().split('\t')[0]
                gene_list = line.strip().split('\t')[2:]
                gene_name_dict[gene_name]=gene_list

        geneset_in_gsva_dict = dict()
        for i in df['geneset']:
            geneset_in_gsva_dict[i] = list(set(gene_name_dict[i]).intersection(set(geneid)))
        df['gsva_id'] =main_id
        detail = df.to_dict('r')
        for each in detail:
            # each['gene'] = geneset_in_gsva_dict[each['geneset']]
            each['seq_list'] = list()
            for i in geneset_in_gsva_dict[each['geneset']]:
                if i in name_geneid_dict:
                    each['seq_list'].append(name_geneid_dict[i])
                elif i.capitalize() in name_geneid_dict:
                    each['seq_list'].append(name_geneid_dict[i.capitalize()])
                else:
                    each['seq_list'].append(i)
            # each['seq_list'] = [name_geneid_dict[i] for i in geneset_in_gsva_dict[each['geneset']]]
            each['seq_num'] = len(each['seq_list'])
            each['description'] = geneset_description[each['geneset']]
        self.create_db_table('sg_geneset_gsva_detail', detail)

        return main_id

    def add_gsva_cluster(self, cluster_output_dir, main_id=None, project_sn='medical_transcriptome',
                            task_id='medical_transcriptome',
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

        if 'gene_name' in samples:
            samples.remove('gene_name')

        if isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        # update main table
        self.update_db_record('sg_geneset_gsva', main_id, gene_cluster=gene_cluster,
                              sample_cluster=sample_cluster)

        tree_info = dict(
            genesets=genes,
            geneset_tree=gene_tree,
            sample_tree=sample_tree,
            gsva_id=main_id,
            samples=samples,
        )
        self.create_db_table('sg_geneset_gsva_cluster_tree', [tree_info])
        self.update_db_record('sg_geneset_gsva', main_id, samples=samples)
    def add_gsva_diff(self, diff_output, main_id, group):
        group_dict = dict()
        cmp_detail_dict = dict()
        with open(group, 'r') as g:
            for line in g.readlines():
                if line.startswith('#'):
                    continue
                else:
                    sample, group = line.strip().split('\t')
                    if group not in group_dict:
                        group_dict[group]=[sample]
                    else:
                        group_dict[group].append(sample)
        diff_files = glob.glob(os.path.join(diff_output, '*_vs_*.*.xls'))
        df_list = list()
        cmp_list = list()
        for each in diff_files:
            diff_pd = pd.read_table(each, header=0, sep='\t')
            fname = os.path.basename(each)
            ctrl, test = re.match('(.*)_vs_(.*).limma.xls', fname).groups()
            cmp_combine = ctrl + '|' + test
            cmp_list.append(cmp_combine)
            diff_pd['compare'] = '{}|{}'.format(ctrl, test)
            need_cols = ['geneset', 'fc', 'log2fc', 'pvalue', 'padjust', 'compare']
            columns = diff_pd.columns
            samples = list()
            for x in columns:
                _m = re.match(r'(.*)_count$', x)
                if _m:
                    samples.append(_m.groups()[0])
            cmp_detail_dict[cmp_combine] = samples
            diff_pd_new = diff_pd.loc[:,need_cols]
            diff_pd_new['group1'] = diff_pd.loc[:, ["{}_count".format(i) for i in group_dict[ctrl]]].mean(axis=1)
            diff_pd_new['group2'] = diff_pd.loc[:, ["{}_count".format(i) for i in group_dict[test]]].mean(axis=1)
            df_list.append(diff_pd_new)
        df_total = pd.concat(df_list)
        df_total['gsva_id'] = main_id
        detail = df_total.to_dict('r')
        self.create_db_table('sg_geneset_gsva_diff', detail)
        self.update_db_record('sg_geneset_gsva', main_id,
                              cmp_list=cmp_list,
                              cmp_detail=cmp_detail_dict,
                              main_id=main_id,
                              )
        self.update_db_record('sg_geneset_gsva', main_id, status='end', main_id=main_id)



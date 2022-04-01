# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# last_modify:20160919
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import bson.binary
from cStringIO import StringIO
import re
import json
import pandas as pd
import numpy as np
from collections import Counter


class DenovoExpress(Base):
    def __init__(self, bind_object):
        super(DenovoExpress, self).__init__(bind_object)
        self._db_name = Config().MONGODB + '_rna'

    @report_check
    def add_express(self, rsem_dir=None, samples=None, params=None, name=None, express_diff_id=None, bam_path=None, major=True, gene_distri=None, tran_distri=None):
        # 参数express_diff_id只是为了插入差异基因矩阵时才用到，初始化不用
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        if express_diff_id:
            params['express_diff_id'] = str(express_diff_id)
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'ExpressStat_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'desc': '表达量计算主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'specimen': samples,
            'status': 'end',
            'bam_path': bam_path,
        }
        collection = self.db['sg_denovo_express']
        express_id = collection.insert_one(insert_data).inserted_id
        if major:
            rsem_files = os.listdir(rsem_dir)
            for f in rsem_files:
                if re.search(r'^genes\.TMM', f):
                    fpkm_path = rsem_dir + f
                    count_path = rsem_dir + 'genes.counts.matrix'
                    self.add_express_detail(express_id, count_path, fpkm_path, 'gene')
                    self.add_express_ditribution(express_id, gene_distri, 'gene')
                elif re.search(r'^transcripts\.TMM', f):
                    fpkm_path = rsem_dir + f
                    count_path = rsem_dir + 'transcripts.counts.matrix'
                    self.add_express_detail(express_id, count_path, fpkm_path, 'transcript')
                    self.add_express_ditribution(express_id, tran_distri, 'transcript')
                elif re.search(r'\.genes\.results$', f):
                    sample = f.split('.genes.results')[0]
                    file_ = rsem_dir + f
                    self.add_express_specimen_detail(express_id, file_, 'gene', sample)
                elif re.search(r'\.isoforms\.results$', f):
                    sample = f.split('.isoforms.results')[0]
                    file_ = rsem_dir + f
                    self.add_express_specimen_detail(express_id, file_, 'transcript', sample)
        return express_id

    @report_check
    def add_express_detail(self, express_id, count_path, fpkm_path, path_type):
        if not isinstance(express_id, ObjectId):
            if isinstance(express_id, types.StringTypes):
                express_id = ObjectId(express_id)
            else:
                raise Exception('express_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(count_path):
            raise Exception('count_path:{}所指定的路径不存在，请检查！'.format(count_path))
        if not os.path.exists(fpkm_path):
            raise Exception('fpkm_path:{}所指定的路径不存在，请检查！'.format(fpkm_path))
        data_list = list()
        count_dict = {}
        # sample_count = {}
        with open(count_path, 'rb') as c, open(fpkm_path, 'rb') as f:
            samples = c.readline().strip().split('\t')
            # for sam in samples:
            #     sample_count[sam] = 0
            for line in c:
                line = line.strip().split('\t')
                count_dict[line[0]] = line[1:]
                # count = line[1:]
                # for i in range(len(count)):
                #     if float(count[i]) > 0:
                #         sample_count[samples[i]] += 1
            f.readline()
            for l in f:
                l = l.strip().split('\t')
                gene_id = l[0]
                fpkm = l[1:]
                data = [
                    ('gene_id', gene_id),
                    ('type', path_type),
                    ('express_id', express_id),
                ]
                for i in range(len(samples)):
                    data += [
                        ('{}_count'.format(samples[i]), float(count_dict[gene_id][i])), ('{}_fpkm'.format(samples[i]), float(fpkm[i])),
                        # ('{}_sum'.format(samples[i]), sample_count[samples[i]]),
                    ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_denovo_express_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入表达量矩阵信息出错:%s" % e)
        else:
            self.bind_object.logger.info("导入表达量矩阵信息成功!")

    @report_check
    def add_express_specimen_detail(self, express_id, rsem_result, rsem_type, sample=None):
        if not isinstance(express_id, ObjectId):
            if isinstance(express_id, types.StringTypes):
                express_id = ObjectId(express_id)
            else:
                raise Exception('express_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(rsem_result):
            raise Exception('rsem_result所指定的路径：{}不存在，请检查！'.format(rsem_result))
        sample_name = os.path.basename(rsem_result).split('.')[0]
        data_list = []
        with open(rsem_result, 'rb') as f:
            f.readline()
            for line in f:
                line = line.strip().split('\t')
                data = [
                    ('express_id', express_id),
                    ('specimen_name', sample if sample else sample_name),
                    ('type', rsem_type),
                    ('length', float(line[2])),
                    ('effective_length', float(line[3])),
                    ('expected_count', float(line[4])),
                    ('TPM', round(float(line[5]), 4)),
                    ('FPKM', round(float(line[6]), 4)),
                ]
                if rsem_type == 'gene':
                    data += [
                        ('gene_id', line[0]),
                        ('transcript_id', line[1]),
                    ]
                else:
                    data += [
                        ('gene_id', line[1]),
                        ('transcript_id', line[0]),
                        ('IsoPct', float(line[7])),
                    ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_denovo_express_specimen_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入单样本表达量矩阵：%s信息出错:%s" % (rsem_result, e))
        else:
            self.bind_object.logger.info("导入单样本表达量矩阵: %s信息成功!" % rsem_result)

    @report_check
    def add_express_diff(self, params, samples, compare_column, diff_exp_dir=None, express_id=None, name=None, group_id=None, group_detail=None, control_id=None, major=True):
        # group_id, group_detail, control_id只供denovobase初始化时更新param使用
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        params.update({
            'express_id': str(express_id),
            'group_id': str(group_id),
            'group_detail': group_detail,
            'control_id': str(control_id)
        })  # 为更新workflow的params，因为截停
        if group_id == 'all':
            params['group_detail'] = {'all': group_detail}
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'ExpressDiffStat_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'desc': '表达量差异检测主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'specimen': samples,
            'status': 'end',
            'compare_column': compare_column,
            'group_detail': group_detail,
            'express_id': express_id
        }
        if group_id == 'all':
            insert_data['group_detail'] = {'all': group_detail}
        collection = self.db['sg_denovo_express_diff']
        express_diff_id = collection.insert_one(insert_data).inserted_id
        if major:
            diff_exp_files = os.listdir(diff_exp_dir)
            for f in diff_exp_files:
                if re.search(r'_edgr_stat.xls$', f):
                    con_exp = f.split('_edgr_stat.xls')[0].split('_vs_')
                    self.add_express_diff_detail(express_diff_id, con_exp, diff_exp_dir + f)
        return express_diff_id

    @report_check
    def add_express_diff_detail(self, express_diff_id, group, diff_stat_path):
        """
        group:为两两比较的样本或分组名，列表
        """
        if not isinstance(express_diff_id, ObjectId):
            if isinstance(express_diff_id, types.StringTypes):
                express_diff_id = ObjectId(express_diff_id)
            else:
                raise Exception('express_diff_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(diff_stat_path):
            raise Exception('diff_stat_path所指定的路径:{}不存在，请检查！'.format(diff_stat_path))
        data_list = []
        with open(diff_stat_path, 'rb') as f:
            head = f.readline().strip().split('\t')
            for line in f:
                line = line.strip().split('\t')
                data = [
                    ('name', group[0]),
                    ('compare_name', group[1]),
                    ('express_diff_id', express_diff_id),
                ]
                for i in range(len(head)):
                    if re.match(r'\d', line[i]):
                        data.append((head[i], float(line[i])))
                    else:
                        data.append((head[i], line[i]))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_denovo_express_diff_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入基因表达差异统计表：%s信息出错:%s" % (diff_stat_path, e))
        else:
            self.bind_object.logger.info("导入基因表达差异统计表：%s信息成功!" % diff_stat_path)

    @report_check
    def add_express_ditribution(self, express_id, ditribution_path, query_type=None):
        df = pd.read_table(ditribution_path)
        samples = df.columns[1:]
        data_list = []
        logfc = df.columns[0]
        for i in samples:
            insert_data = [
                ('express_id', express_id),
                ('type', query_type),
                ('specimen', i)
            ]
            data = []
            for v in df[i].index:
                data.append({'logfpkm': df[logfc][v], 'density': round(float(df[i][v]), 4)})
            insert_data.append(('data', data))
            insert_data = SON(insert_data)
            data_list.append(insert_data)
        try:
            collection = self.db["sg_denovo_express_gragh"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入表达量矩阵作图数据：%s信息出错:%s" % (ditribution_path, e))
        else:
            self.bind_object.logger.info("导入表达量矩阵作图数据: %s信息成功!" % ditribution_path)

    @report_check
    def add_express_gragh(self, express_id, fpkm_path, query_type=None):
        df = pd.read_table(fpkm_path)
        samples = df.columns[1:]
        data_list = []
        for i in samples:
            insert_data = [
                ('express_id', express_id),
                ('type', query_type),
                ('specimen', i)
            ]
            tmp = []
            data = df[i][df[i].apply(lambda x: x > 0)]
            data = np.log10(data).apply(lambda x: round(x, 2))
            count = Counter(data)
            _sum = float(sum(count.values()))
            for c in count:
                tmp.append({'logfpkm': c, 'density': round((count[c] / _sum), 4)})
            insert_data.append(('data', tmp))
            insert_data = SON(insert_data)
            data_list.append(insert_data)
        try:
            collection = self.db["sg_denovo_express_gragh"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入表达量矩阵作图数据：%s信息出错:%s" % (fpkm_path, e))
        else:
            self.bind_object.logger.info("导入表达量矩阵作图数据: %s信息成功!" % fpkm_path)

    # @report_check
    # def add_cluster(self, params, express_id, sample_tree=None, gene_tree=None, name=None, samples=None, genes=None):
    #     if not isinstance(express_id, ObjectId):
    #         if isinstance(express_id, types.StringTypes):
    #             express_id = ObjectId(express_id)
    #     task_id = self.bind_object.sheet.id
    #     project_sn = self.bind_object.sheet.project_sn
    #     if gene_tree:
    #         with open(gene_tree, 'rb') as g:
    #             gene_tree = g.readlines()[0].strip('\n')
    #     if sample_tree:
    #         with open(sample_tree, 'rb') as s:
    #             sample_tree = s.readlines()[0].strip('\n')
    #     params['diff_fpkm'] = str(express_id)
    #     insert_data = {
    #         'project_sn': project_sn,
    #         'task_id': task_id,
    #         'name': name if name else 'cluster_table_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
    #         'desc': '差异基因聚类分析主表',
    #         'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
    #         'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
    #         'specimen': samples,
    #         'status': 'end',
    #         'express_id': express_id,
    #         'sample_tree': sample_tree,
    #         'gene_tree': gene_tree,
    #         'genes': genes,
    #     }
    #     collection = self.db['sg_denovo_cluster']
    #     cluster_id = collection.insert_one(insert_data).inserted_id
    #     return cluster_id
    #
    # @report_check
    # def add_cluster_detail(self, cluster_id, sub, sub_path):
    #     if not isinstance(cluster_id, ObjectId):
    #         if isinstance(cluster_id, types.StringTypes):
    #             cluster_id = ObjectId(cluster_id)
    #         else:
    #             raise Exception('cluster_id必须为ObjectId对象或其对应的字符串！')
    #     if not os.path.exists(sub_path):
    #         raise Exception('sub_path所指定的路径:{}不存在，请检查！'.format(sub_path))
    #     data_list = []
    #     with open(sub_path, 'rb') as f:
    #         head = f.readline().strip().split('\t')
    #         for line in f:
    #             line = line.strip().split('\t')
    #             data = [
    #                 ('sub_cluster', sub),
    #                 ('cluster_id', cluster_id),
    #                 ('gene_id', line[0])
    #             ]
    #             for i in range(len(head)):
    #                 data.append((head[i], line[i + 1]))
    #             data = SON(data)
    #             data_list.append(data)
    #     try:
    #         collection = self.db["sg_denovo_cluster_detail"]
    #         collection.insert_many(data_list)
    #     except Exception, e:
    #         self.bind_object.set_error("导入子聚类统计表：%s信息出错:%s" % (sub_path, e))
    #     else:
    #         self.bind_object.logger.info("导入子聚类统计表:%s信息成功!" % sub_path)
    #
    # @report_check
    # def add_network(self, params, express_id, softpower, module, name=None):
    #     if not isinstance(express_id, ObjectId):
    #         if isinstance(express_id, types.StringTypes):
    #             express_id = ObjectId(express_id)
    #         else:
    #             raise Exception('express_matrix_id必须为ObjectId对象或其对应的字符串！')
    #     if not os.path.exists(softpower):
    #         raise Exception('softpower所指定的路径:{}不存在，请检查！'.format(softpower))
    #     if not os.path.exists(module):
    #         raise Exception('module所指定的路径:{}不存在，请检查！'.format(module))
    #     task_id = self.bind_object.sheet.id
    #     project_sn = self.bind_object.sheet.project_sn
    #     collection = self.db['sg_denovo_network']
    #     with open(softpower, 'rb') as s, open(module, 'rb') as m:
    #         softpower_id = StringIO(s.read())
    #         softpower_id = bson.binary.Binary(softpower_id.getvalue())
    #         module_id = StringIO(m.read())
    #         module_id = bson.binary.Binary(module_id.getvalue())
    #     params['diff_fpkm'] = str(express_id)
    #     insert_data = {
    #         'project_sn': project_sn,
    #         'task_id': task_id,
    #         'name': name if name else 'network_table_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
    #         'desc': '差异基因网络分析主表',
    #         'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
    #         'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
    #         'status': 'end',
    #         'express_id': express_id,
    #         'softpower': softpower_id,
    #         'module': module_id,
    #     }
    #     network_id = collection.insert_one(insert_data).inserted_id
    #     return network_id
    #
    # @report_check
    # def add_network_detail(self, network_id, node_path, edge_path):
    #     if not isinstance(network_id, ObjectId):
    #         if isinstance(network_id, types.StringTypes):
    #             network_id = ObjectId(network_id)
    #     if not os.path.exists(node_path):
    #         raise Exception('node_path所指定的路径:{}不存在，请检查！'.format(node_path))
    #     if not os.path.exists(edge_path):
    #         raise Exception('edge_path所指定的路径:{}不存在，请检查！'.format(edge_path))
    #     data_list = []
    #     gene_color = {}
    #     with open(node_path, 'rb') as n, open(edge_path, 'rb') as f:
    #         n.readline()
    #         for line in n:
    #             line = line.strip().split('\t')
    #             gene_color[line[0]] = line[2]
    #         f.readline()
    #         for line in f:
    #             line = line.strip().split('\t')
    #             data = [
    #                 ('network_id', network_id),
    #                 ('gene_id1', {'name': line[0], 'color': gene_color[line[0]]}),
    #                 ('gene_id2', {'name': line[1], 'color': gene_color[line[1]]}),
    #                 ('weight', line[2]),
    #             ]
    #             data = SON(data)
    #             data_list.append(data)
    #     try:
    #         collection = self.db["sg_denovo_network_detail"]
    #         collection.insert_many(data_list)
    #     except Exception, e:
    #         self.bind_object.set_error("导入网络表达统计表：%s，%s信息出错:%s" % (node_path, edge_path, e))
    #     else:
    #         self.bind_object.logger.info("导入网络表达统计表:%s， %s信息成功!" % (node_path, edge_path))
    #
    # @report_check
    # def add_network_module(self, network_id, module_path, module_color):
    #     if not isinstance(network_id, ObjectId):
    #         if isinstance(network_id, types.StringTypes):
    #             network_id = ObjectId(network_id)
    #     if not os.path.exists(module_path):
    #         raise Exception('module_path所指定的路径:{}不存在，请检查！'.format(module_path))
    #     data_list = []
    #     with open(module_path, 'rb') as f:
    #         f.readline()
    #         for line in f:
    #             line = line.strip().split('\t')
    #             data = [
    #                 ('network_id', network_id),
    #                 ('gene_id1', line[0]),
    #                 ('gene_id2', line[1]),
    #                 ('weight', line[2]),
    #                 ('module_color', module_color),
    #             ]
    #             data = SON(data)
    #             data_list.append(data)
    #     try:
    #         collection = self.db["sg_denovo_network_module"]
    #         collection.insert_many(data_list)
    #     except Exception, e:
    #         self.bind_object.set_error("导入网络表达统计表：%s信息出错:%s" % (module_path, e))
    #     else:
    #         self.bind_object.logger.info("导入网络表达统计表:%s信息成功!" % module_path)

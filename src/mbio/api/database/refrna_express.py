# !/usr/bin/python
# -*- coding: utf-8 -*-
# __author__: konghualei 20170424

from pymongo import MongoClient
from bson.objectid import ObjectId
import types
from types import StringTypes
import re
import json, time
import pandas as pd
import numpy as np
import datetime, os
from bson.son import SON
from collections import Counter
import glob
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import math
from math import log10
import pandas as pd
#from gevent.monkey import patch_all
import gevent
from mainapp.libs.param_pack import group_detail_sort
#patch_all()


class RefrnaExpress(Base):
    def __init__(self, bind_object):
        super(RefrnaExpress, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + '_ref_rna'
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]

    def get_gene_name(self, class_code, query_type=None, workflow=False):
        """
        :params: 是否工作流根据class_code信息导入基因/转录本名称
        对转录本都加上相应的gene id信息
        """
        if query_type not in ["gene", "transcript"]:
            raise Exception("query type should be gene or transcript")
        with open(class_code) as f1:
            f1.readline()
            data = dict()
            for lines in f1:
                line = lines.strip('\n').split("\t")
                if workflow:
                    if not line[3]:
                        line[3] = '-'
                    if not line[1]:
                        raise Exception('{} has no gene_id in {}'.format(line[0], class_code))
                    t_id, gene_id, class_code_type, gene_name = line
                    if query_type == 'transcript':
                        data[t_id] = dict(gene_name=gene_name, gene_id=gene_id, class_code=class_code_type)
                    else:
                        data[gene_id] = dict(gene_name=gene_name, class_code=class_code_type)
                else:
                    if len(line) == 4:
                        if not line[3]:
                            line[3] = '-'
                        if not line[1]:
                            raise Exception('{} has no gene_id in {}'.format(line[0], class_code))
                        t_id, gene_name, class_code_type, gene_id = line
                        data[t_id] = dict(gene_name=gene_name, gene_id=gene_id)
                    else:
                        if not line[1]:
                            line[1] = '-'
                        gene_id, gene_name, class_code_type = line
                        data[gene_id] = dict(gene_name=gene_name)
        return data

    # @report_check
    def add_express_gragh(self, express_id, distribution_path_log2, distribution_path_log10, distribution_path,
                          sample_group,
                          query_type=None, value_type=None):
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]

        def stat(fpkm_data, density, log=None):
            tmp = []
            if len(fpkm_data) != len(density):
                raise Exception("density必须和fpkm长度相等！")
            else:
                for i in range(len(fpkm_data)):
                    if not log:
                        tmp.append({"fpkm": round(fpkm_data[i], 6), "density": round(density[i], 6)})
                return tmp

        dflog2 = pd.read_table(distribution_path_log2)
        dflog10 = pd.read_table(distribution_path_log10)
        df = pd.read_table(distribution_path)
        samples = df.columns[1:]
        data_list = []
        for i in samples:
            insert_data = [
                ('express_id', express_id),
                ('type', query_type),
                ('specimen', i),
                ('sample_group', sample_group),
                ('value_type', value_type)
            ]
            tmp = stat(fpkm_data=df["fpkm"], density=df[i])
            tmplog2 = stat(fpkm_data=dflog2["log2fpkm"], density=dflog2[i])
            tmplog10 = stat(fpkm_data=dflog10["log10fpkm"], density=dflog10[i])
            insert_data.append(('data', tmp))
            insert_data.append(('data_log2', tmplog2))
            insert_data.append(('data_log10', tmplog10))
            insert_data = SON(insert_data)
            data_list.append(insert_data)
        try:
            collection = self.db["sg_express_gragh"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入表达量矩阵作图数据：%s信息出错:%s" % (distribution_path_log2, e))
        else:
            self.bind_object.logger.info("导入表达量矩阵作图数据: %s信息成功!" % distribution_path_log2)

    # @report_check
    def add_express(self, rsem_dir=None, group_fpkm_path=None, transcript_fasta_path=None, is_duplicate=None,
                    class_code=None, samples=None, params=None, name=None, express_diff_id=None, bam_path=None,
                    major=True, distri_path=None):

        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        # if not express_diff_id:
        # params={"value_type":value_type, "query_type":query_type,"method":method, "group_id":group_id,"group_detail":group_detail}
        if params:
            params["submit_location"] = "express_rsem"
        value_type = params["type"]
        if 'express_method' not in params.keys():
            raise Exception("请在params中设置express_method参数!")
        express_method = params['express_method']
        print "value_type"
        print value_type
        re_name = "ExpStat_RSEM_{}_".format(value_type.lower()) + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else re_name,
            'desc': '表达量计算主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'specimen': samples,
            'status': 'end',
            'bam_path': bam_path,
            'transcript_fasta_path': transcript_fasta_path,
            'is_duplicate': is_duplicate
        }
        if params:
            insert_data["genes"] = True
            if params["express_method"] == "rsem":
                insert_data["trans"] = True
            if params["express_method"].lower() == "featurecounts":
                insert_data["trans"] = False
            if params['group_detail']:
                if is_duplicate:
                    insert_data["group"] = params["group_detail"].keys()
        if express_diff_id:
            insert_data["express_diff_id"] = express_diff_id
        collection = self.db['sg_express']

        # if params:
        #    params.update({"submit_location": "express_rsem_{}".format(value_type)})
        express_id = collection.insert_one(insert_data).inserted_id
        print "插入主表id是{}".format(express_id)

        sample_group = "sample"
        method = params["express_method"]
        # express_id=ObjectId("58f03a28a4e1af44d4139c79")
        if major:
            gevent_list = []
            rsem_files = os.listdir(rsem_dir)
            # sample_group = "sample"
            for f in rsem_files:
                if re.search(r'^genes\.TMM', f):
                    fpkm_path = rsem_dir + "/" + f
                    count_path = rsem_dir + '/genes.counts.matrix'
                    # query_type=None,value_type=None, method=None, sample_group=None
                    print fpkm_path
                    print count_path
                    gevent_list.append(gevent.spawn(self.add_express_detail, express_id, fpkm_path, class_code, 'gene', value_type, method,
                                            sample_group))
                    gevent_list.append(gevent.spawn(self.add_express_detail, express_id, count_path, class_code, 'gene', 'count', method,
                                            sample_group))
                    try:
                        gevent_list.append(gevent.spawn(self.add_express_gragh, express_id,
                                               distribution_path_log2=distri_path + "/log2gene_distribution.xls",
                                               distribution_path_log10=distri_path + "/log10gene_distribution.xls",
                                               distribution_path=distri_path + "/gene_distribution.xls",
                                               sample_group="sample",
                                               query_type="gene", value_type=value_type))
                    except Exception:
                        print 'error!'
                    print '导入基因{}graph成功'.format(value_type)

                    gevent_list.append(gevent.spawn(self.add_express_gragh, express_id,
                                           distribution_path_log2=distri_path + "/sample_count_distribution/log2gene_distribution.xls",
                                           distribution_path_log10=distri_path + "/sample_count_distribution/log10gene_distribution.xls",
                                           distribution_path=distri_path + "/sample_count_distribution/gene_distribution.xls",
                                           sample_group='sample',
                                           query_type="gene", value_type='count'))
                    print '导入基因{}graph成功'.format('count')

                    gevent_list.append(gevent.spawn(self.add_express_box, express_id, fpkm_path=fpkm_path, sample_group="sample", query_type="gene",
                                         value_type=value_type))
                    gevent_list.append(gevent.spawn(self.add_express_box, express_id, fpkm_path=count_path, sample_group="sample", query_type="gene",
                                         value_type='count'))

                    if is_duplicate:
                        if value_type == 'fpkm':
                            gene_group_fpkm_path = distri_path + "/Group.genes_genes.TMM.fpkm.matrix"
                            gene_group_count_path = distri_path + "/Group.genes_count_genes.TMM.fpkm.matrix"
                        if value_type == 'tpm':
                            gene_group_fpkm_path = distri_path + "/Group.genes_genes.TMM.EXPR.matrix"
                            gene_group_count_path = distri_path + "/Group.genes_count_genes.TMM.EXPR.matrix"
                        print '开始进行group导表'
                        if os.path.exists(gene_group_fpkm_path):
                            gevent_list.append(gevent.spawn(self.add_express_group_detail, express_id, gene_group_fpkm_path, "gene", value_type, "rsem",
                                                          "group"))
                            print '基因group{}导表成功！'.format(value_type)
                        if os.path.exists(gene_group_count_path):
                            gevent_list.append(gevent.spawn(self.add_express_group_detail, express_id, gene_group_count_path, "gene", "count", "rsem",
                                                          "group"))
                            print '基因group count导表成功！'
                        gevent_list.append(gevent.spawn(self.add_express_gragh, express_id,
                                               distribution_path_log2=distri_path + "/group_genes_fpkm_distribution/log2GroupGenes_distribution.xls",
                                               distribution_path_log10=distri_path + "/group_genes_fpkm_distribution/log10GroupGenes_distribution.xls",
                                               distribution_path=distri_path + "/group_genes_fpkm_distribution/GroupGenes_distribution.xls",
                                               sample_group="group", query_type="gene", value_type=value_type))
                        print '基因group graph{}导表成功！'.format(value_type)
                        gevent_list.append(gevent.spawn(self.add_express_gragh, express_id,
                                               distribution_path_log2=distri_path + "/group_genes_count_distribution/log2GroupGenes_distribution.xls",
                                               distribution_path_log10=distri_path + "/group_genes_count_distribution/log10GroupGenes_distribution.xls",
                                               distribution_path=distri_path + "/group_genes_count_distribution/GroupGenes_distribution.xls",
                                               sample_group="group", query_type="gene", value_type='count'))
                        print '基因group graph{}导表成功！'.format('count')
                        # #目前缺少count的graph数据
                        gevent_list.append(gevent.spawn(self.add_express_box, express_id, fpkm_path=gene_group_fpkm_path,
                                             sample_group="group", query_type="gene", value_type=value_type))
                        gevent_list.append(gevent.spawn(self.add_express_box, express_id, fpkm_path=gene_group_count_path, sample_group='group',
                                             query_type='gene', value_type='count'))
                elif re.search(r'^transcripts\.TMM', f):
                    fpkm_path = rsem_dir + "/" + f
                    count_path = rsem_dir + '/transcripts.counts.matrix'
                    print fpkm_path
                    print count_path

                    gevent_list.append(gevent.spawn(self.add_express_detail, express_id, fpkm_path, class_code, 'transcript', value_type, method,
                                            sample_group))
                    gevent_list.append(gevent.spawn(self.add_express_detail, express_id, count_path, class_code, 'transcript', 'count', method,
                                            sample_group))
                    gevent_list.append(gevent.spawn(self.add_express_gragh, express_id,
                                           distribution_path_log2=distri_path + "/log2transcript_distribution.xls",
                                           distribution_path_log10=distri_path + "/log10transcript_distribution.xls",
                                           distribution_path=distri_path + "/transcript_distribution.xls",
                                           sample_group="sample", query_type="transcript", value_type=value_type))
                    gevent_list.append(gevent.spawn(self.add_express_gragh, express_id,
                                           distribution_path_log2=distri_path + "/sample_count_distribution/log2transcript_distribution.xls",
                                           distribution_path_log10=distri_path + "/sample_count_distribution/log10transcript_distribution.xls",
                                           distribution_path=distri_path + "/sample_count_distribution/transcript_distribution.xls",
                                           sample_group="sample", query_type="transcript", value_type='count'))
                    gevent_list.append(gevent.spawn(self.add_express_box,express_id, fpkm_path=fpkm_path, sample_group="sample",
                                         query_type="transcript", value_type=value_type))
                    gevent_list.append(gevent.spawn(self.add_express_box,express_id, fpkm_path=count_path, sample_group="sample",
                                         query_type="transcript", value_type='count'))
                    if is_duplicate:
                        print value_type

                        if value_type == 'fpkm':
                            trans_group_fpkm_path = distri_path + "/Group.trans_transcripts.TMM.fpkm.matrix"
                            trans_group_count_path = distri_path + "/Group.trans_count_transcripts.TMM.fpkm.matrix"  # trans的count数据
                        if value_type == 'tpm':
                            trans_group_fpkm_path = distri_path + "/Group.trans_transcripts.TMM.EXPR.matrix"
                            trans_group_count_path = distri_path + "/Group.trans_count_transcripts.TMM.EXPR.matrix"  # trans的count数据
                        print trans_group_fpkm_path

                        if os.path.exists(trans_group_fpkm_path):
                            gevent_list.append(gevent.spawn(self.add_express_group_detail, express_id, trans_group_fpkm_path, "transcript", value_type,
                                                          "rsem", "group"))
                            gevent_list.append(gevent.spawn(self.add_express_group_detail, express_id, trans_group_fpkm_path, "transcript", 'count',
                                                          "rsem", "group"))
                        gevent_list.append(gevent.spawn(self.add_express_gragh, express_id,
                                               distribution_path_log2=distri_path + "/group_trans_fpkm_distribution/log2GroupTrans_distribution.xls",
                                               distribution_path_log10=distri_path + "/group_trans_fpkm_distribution/log10GroupTrans_distribution.xls",
                                               distribution_path=distri_path + "/group_trans_fpkm_distribution/GroupTrans_distribution.xls",
                                               sample_group="group", query_type="transcript", value_type=value_type))
                        gevent_list.append(gevent.spawn(self.add_express_gragh, express_id,
                                               distribution_path_log2=distri_path + "/group_trans_count_distribution/log2GroupTrans_distribution.xls",
                                               distribution_path_log10=distri_path + "/group_trans_count_distribution/log10GroupTrans_distribution.xls",
                                               distribution_path=distri_path + "/group_trans_count_distribution/GroupTrans_distribution.xls",
                                               sample_group="group", query_type="transcript", value_type='count'))
                        # 目前缺少trans count的graph数据
                        gevent_list.append(gevent.spawn(self.add_express_box, express_id,
                                             fpkm_path=trans_group_fpkm_path,
                                             sample_group="group", query_type="transcript", value_type=value_type))
                        gevent_list.append(gevent.spawn(self.add_express_box, express_id,
                                             fpkm_path=trans_group_count_path,
                                             sample_group="group", query_type="transcript", value_type='count'))
        gevent.joinall(gevent_list)
        return express_id

    # @report_check
    def add_express_feature(self, feature_dir=None, group_fpkm_path=None, is_duplicate=None, class_code=None,
                            samples=None,
                            params=None, name=None, express_diff_id=None, bam_path=None, major=True,
                            distri_path=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        value_type = params["type"]  # fpkm或者是tpm
        if params:
            params["submit_location"] = "express_feature"
        if 'express_method' not in params.keys():
            raise Exception("请在params中设置express_method参数!")
        express_method = params['express_method']
        re_name = "ExpStat_FeaCount_{}_".format(value_type.lower()) + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else re_name,
            'desc': '表达量计算主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params': (
            json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'specimen': samples,
            'status': 'end',
            'bam_path': bam_path,
            # 'transcript_fasta_path': transcript_fasta_path,
            'is_duplicate': is_duplicate
        }
        sample_group = "sample"
        if params:
            insert_data["genes"] = True
            insert_data["trans"] = False
            if params['group_detail']:
                if is_duplicate:
                    insert_data["group"] = params["group_detail"].keys()
        collection = self.db['sg_express']

        method = params["express_method"]
        # if params:
        #    params.update({"submit_location": "express_feature"})
        express_id = collection.insert_one(insert_data).inserted_id
        print "插入主表id是{}".format(express_id)

        # express_id=ObjectId("58f03a28a4e1af44d4139c79")
        if major:
            gevent_list = []
            if value_type == 'fpkm':
                fpkm_path = feature_dir + "/fpkm_tpm.fpkm.xls"
            if value_type == 'tpm':
                fpkm_path = feature_dir + "/fpkm_tpm.tpm.xls"
            count_path = feature_dir + "/count.xls"
            gevent_list.append(gevent.spawn(self.add_express_detail, express_id, fpkm_path, class_code, 'gene', value_type, method, sample_group))
            gevent_list.append(gevent.spawn(self.add_express_detail, express_id, count_path, class_code, 'gene', 'count', method, sample_group))
            gevent_list.append(gevent.spawn(self.add_express_gragh, express_id,
                                   distribution_path_log2=distri_path + "/{}/log2gene_distribution.xls".format(
                                       value_type),
                                   distribution_path_log10=distri_path + "/{}/log10gene_distribution.xls".format(
                                       value_type),
                                   distribution_path=distri_path + "/{}/gene_distribution.xls".format(value_type),
                                   sample_group="sample", query_type="gene", value_type=value_type))
            gevent_list.append(gevent.spawn(self.add_express_gragh, express_id,
                                   distribution_path_log2=distri_path + "/{}/log2gene_distribution.xls".format('count'),
                                   distribution_path_log10=distri_path + "/{}/log10gene_distribution.xls".format(
                                       "count"),
                                   distribution_path=distri_path + "/{}/gene_distribution.xls".format("count"),
                                   sample_group="sample", query_type="gene", value_type='count'))
            gevent_list.append(gevent.spawn(self.add_express_box, express_id, fpkm_path=fpkm_path, sample_group="sample", query_type="gene",
                                 value_type=value_type))
            #self.add_express_box(express_id, fpkm_path=count_path, sample_group="sample", query_type="gene",
            #                     value_type='count')
            if is_duplicate:
                if value_type == 'fpkm':
                    fpkm_group_path = group_fpkm_path + "/group.fpkm.xls"
                if value_type == 'tpm':
                    fpkm_group_path = group_fpkm_path + "/group.tpm.xls"
                count_group_path = distri_path + "/output/fpkm_count.xls"
                gevent_list.append(gevent.spawn(self.add_express_group_detail, express_id, fpkm_group_path, "gene", value_type, "featurecounts", "group"))
                gevent_list.append(gevent.spawn(self.add_express_group_detail, express_id, count_group_path, "gene", 'count', "featurecounts", "group"))
                gevent_list.append(gevent.spawn(self.add_express_gragh, express_id,
                                       distribution_path_log2=distri_path + "/group_fpkm_distribution/log2GroupGenes_distribution.xls",
                                       distribution_path_log10=distri_path + "/group_fpkm_distribution/log10GroupGenes_distribution.xls",
                                       distribution_path=distri_path + "/group_fpkm_distribution/GroupGenes_distribution.xls",
                                       sample_group="group", query_type="gene", value_type=value_type))
                gevent_list.append(gevent.spawn(self.add_express_gragh, express_id,
                                       distribution_path_log2=distri_path + "/group_count_distribution/log2GroupGenes_distribution.xls",
                                       distribution_path_log10=distri_path + "/group_count_distribution/log10GroupGenes_distribution.xls",
                                       distribution_path=distri_path + "/group_count_distribution/GroupGenes_distribution.xls",
                                       sample_group="group", query_type="gene", value_type='count'))
                gevent_list.append(gevent.spawn(self.add_express_box, express_id, fpkm_path=distri_path + "/group/group.{}.xls".format(value_type),
                                     sample_group="group", query_type="gene", value_type=value_type))
                gevent_list.append(gevent.spawn(self.add_express_box, express_id, fpkm_path=distri_path + "/output/fpkm_count.xls",
                                     sample_group="group", query_type="gene", value_type='count'))
        gevent.joinall(gevent_list)
        return express_id

    # @report_check
    def add_express_specimen_detail_feature(self, express_id, feature_result):
        """featurecounts单个样本的表达量信息"""
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        if not isinstance(express_id, ObjectId):
            if isinstance(express_id, types.StringTypes):
                express_id = ObjectId(express_id)
            else:
                raise Exception('express_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(feature_result):
            raise Exception('rsem_result所指定的路径：{}不存在，请检查！'.format(feature_result))
        with open(feature_result, 'r+') as f1:
            header = f1.readline().strip().split("\t")
            sample = header[6:]
            sample_num = len(sample)
            data_list = []
            for lines in f1:
                line = lines.strip().split("\t")
                insert_data = []
                for i in range(len(header)):
                    insert_data += [(header[i], line[i])]
                insert_data = SON(insert_data)
                data_list.append(insert_data)
        try:
            collection = self.db["sg_express_specimen_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入单样本表达量矩阵：%s信息出错:%s" % (feature_result, e))
        else:
            self.bind_object.logger.info("导入单样本表达量矩阵: %s信息成功" % (feature_result))

    # @report_check
    def add_express_group_detail(self, express_id, group_fpkm_path=None, query_type=None, value_type=None, method=None,
                                 sample_group=None):
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        from math import log
        if not isinstance(express_id, ObjectId):
            if isinstance(express_id, types.StringTypes):
                express_id = ObjectId(express_id)
            else:
                raise Exception('express_id必须为ObjectId对象或其对应的字符串！')
        with open(group_fpkm_path, 'r+') as f1:
            data_list = []
            group_name = f1.readline().strip().split("\t")
            group_name.sort()
            group_num = len(group_name)
            if not query_type:
                raise Exception("请设置query_type参数！")
            if not method:
                raise Exception("请设置method参数！")
            if not value_type:
                raise Exception("请设置value_type参数！")
            if not sample_group:
                raise Exception("请设置sample_group参数！")
            for lines in f1:
                line = lines.strip().split("\t")
                seq_id = line[0]
                insert_data = [
                    ("type", query_type),
                    ("value_type", value_type),
                    ("method", method),
                    ("sample_group", sample_group),
                    ("express_id", express_id),
                    ("seq_id", line[0])
                ]
                fpkm_data = line[1:]
                for i in range(len(fpkm_data)):
                    if float(fpkm_data[i]) >= 1e-02 and float(fpkm_data[i]) <= 1e+06:
                        # print '{}fpkm_data{}符合要求'.format(seq_id, str(float(fpkm_data[i])))
                        data_log2 = log(float(fpkm_data[i])) / log(2)
                        data_log10 = log(float(fpkm_data[i])) / log(10)
                        insert_data += [
                            ('{}_log2'.format(group_name[i]), float(data_log2)),
                            ('{}_log10'.format(group_name[i]), float(data_log10))
                        ]
                    else:
                        continue

                insert_data = SON(insert_data)
                data_list.append(insert_data)
        try:
            collection = self.db["sg_express_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入表达量矩阵group信息出错:%s" % e)
            # print ("导入表达量矩阵信息出错:%s" % e)
        else:
            # print ("导入表达量矩阵信息成功!")
            self.bind_object.logger.info("导入表达量矩阵group信息成功!")

    # @report_check
    def add_express_detail(self, express_id, path, class_code=None, query_type=None, value_type=None,
                                   method=None, sample_group=None, diff=True, workflow=True):
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        if not isinstance(express_id, ObjectId):
            if isinstance(express_id, types.StringTypes):
                express_id = ObjectId(express_id)
            else:
                raise Exception('express_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(path):
            raise Exception('path:{}所指定的路径不存在，请检查！'.format(path))

        # def class_code_get(class_code, query_type):
        #     if class_code:
        #         if os.path.exists(class_code):
        #             with open(class_code, 'r+') as cc:
        #                 class_code_dict = {}
        #                 for lines in cc:
        #                     line = lines.strip().split("\t")
        #                     if query_type == "gene":
        #                         if line[1] not in class_code_dict.keys():
        #                             class_code_dict[line[1]] = str(line[2])
        #                         else:
        #                             pass
        #                     if query_type == "transcript":
        #                         if line[0] not in class_code_dict.keys():
        #                             class_code_dict[line[0]] = str(line[2])
        #                         else:
        #                             pass
        #             return class_code_dict

        data_list = list()
        count_dict = {}
        sample_count = {}
        if class_code:
            class_code_info = self.get_gene_name(class_code, query_type, workflow=True)
            # class_code_get(class_code=class_code, query_type=query_type)
        min_log2 = float(math.log(1e-02) / math.log(2))
        max_log2 = float(math.log(1e+06) / math.log(2))
        min_log10 = float(math.log(1e-02) / math.log(10))
        max_log10 = float(math.log(1e+06) / math.log(10))
        with open(path, 'rb') as f:
            samples = f.readline().strip().split("\t")
            for l in f:
                l = l.strip().split('\t')
                sequence_id = l[0]
                # if re.search(r'(,)', seq_id):
                #     sequence_id = seq_id.split(",")[0]  # 以 ',' 为分隔符切割序列id和gene_name
                #     gene_name = seq_id.split(",")[1]
                # else:
                #     sequence_id = seq_id
                #     gene_name = None
                if class_code:
                    if sequence_id in class_code_info.keys():
                        _class_code = str(class_code_info[sequence_id]["class_code"])
                        if query_type == 'gene':
                            if _class_code == "u":
                                _class = True
                            else:
                                _class = False
                        if query_type == 'transcript':
                            _gene_id = class_code_info[sequence_id]['gene_id']
                            if _class_code == '=':
                                _class = False
                            else:
                                _class = True
                    else:
                        _class = None
                        if query_type == 'transcript':
                            _gene_id = '-'

                fpkm = l[1:]
                # if class_code:
                #     if sequence_id in class_code_info.keys():
                #         _class_code = class_code_info[sequence_id]
                #         if _class_code != "=":
                #             _class = True
                #         else:
                #             _class = False
                #     else:
                #         _class = None
                # else:
                #     _class = None

                data = [
                    ('seq_id', sequence_id),
                    ('type', query_type),
                    ('express_id', express_id),
                    ("value_type", value_type),
                    ("method", method),
                    # ("gene_name", gene_name),  # 添加gene_name信息
                    ("sample_group", sample_group)
                ]
                if class_code:
                    data += [
                        ("is_new", _class),
                        ("gene_name", class_code_info[sequence_id]["gene_name"])
                    ]
                if query_type == 'transcript':
                    data += [("gene_id", _gene_id)]
                for i in range(len(samples)):
                    log2_line_fpkm = math.log(float(fpkm[i]) + 1) / math.log(2)
                    log10_line_fpkm = math.log(float(fpkm[i]) + 1) / math.log(10)
                    data += [('{}_line_log2'.format(samples[i]), float(log2_line_fpkm))]
                    data += [('{}_line_log10'.format(samples[i]), float(log10_line_fpkm))]
                    if float(fpkm[i]) < (1e-02) or float(fpkm[i]) > (1e+06):
                        data += [('{}'.format(samples[i]), float(fpkm[i]))]
                        continue
                    else:
                        log2_fpkm = math.log(float(fpkm[i])) / math.log(2)
                        log10_fpkm = math.log(float(fpkm[i])) / math.log(10)
                        data += [('{}'.format(samples[i]), float(fpkm[i]))]
                        if float(log2_fpkm) >= min_log2 and float(log2_fpkm) <= max_log2:
                            data += [('{}_log2'.format(samples[i]), float(log2_fpkm))]
                        if float(log10_fpkm) >= min_log10 and float(log10_fpkm) <= max_log10:
                            data += [('{}_log10'.format(samples[i]), float(log10_fpkm))]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_express_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入表达量矩阵信息出错:%s" % e)
        else:
            self.bind_object.logger.info("导入表达量矩阵信息成功!")

    # @report_check
    def add_express_box(self, express_id, fpkm_path, sample_group, query_type=None, value_type=None):
        if not isinstance(express_id, ObjectId):
            if isinstance(express_id, types.StringTypes):
                express_id = ObjectId(express_id)
            else:
                raise Exception('express_id必须为ObjectId对象或其对应的字符串！')
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]

        min_log2 = float(math.log(1e-02) / math.log(2))
        max_log2 = float(math.log(1e+06) / math.log(2))
        min_log10 = float(math.log(1e-02) / math.log(10))
        max_log10 = float(math.log(1e+06) / math.log(10))

        def log_value(value, genes, log):
            """获取log值"""
            _id = value.apply(lambda x: x >= 1e-02 and x <= 1e+06)
            if log == 2:
                data = value[_id]
                new_genes = genes[_id]
                log2_data = np.log2(data.apply(lambda x: x))
                return log2_data, new_genes
            elif log == 10:
                data = value[_id]
                new_genes = genes[_id]
                log10_data = np.log10(data.apply(lambda x: x))
                return log10_data, new_genes
            else:
                return value, genes

        def box_info(fpkm, samples, log=None):
            box = {}
            gene_list = {}
            for sam in samples:
                gene_list[sam] = {}
                box[sam] = {}
                log_data, seq_id = log_value(fpkm[sam], fpkm[[0]], log)
                min = log_data.min()
                max = log_data.max()
                q1 = log_data.quantile(0.25)
                q3 = log_data.quantile(0.75)
                iqr = q3 - q1
                """过滤掉异常大的数值"""
                # max = q3 + 1.5 * iqr
                # min = q1 - 1.5 * iqr
                median = log_data.median()
                box[sam] = {"min": min, "max": max, 'q1': q1, 'q3': q3, 'median': median}
                min_q1 = seq_id[log_data.apply(lambda x: x >= min and x <= q1)].values
                gene_list[sam]['min-q1'] = [i[0] for i in min_q1]
                q1_median = seq_id[log_data.apply(lambda x: x > q1 and x <= median)].values
                gene_list[sam]['q1-median'] = [i[0] for i in q1_median]
                median_q3 = seq_id[log_data.apply(lambda x: x > median and x <= q3)].values
                gene_list[sam]['median-q3'] = [i[0] for i in median_q3]
                q3_max = seq_id[log_data.apply(lambda x: x > q3 and x <= max)].values
                gene_list[sam]['q3-max'] = [i[0] for i in q3_max]
            return box, gene_list

        express_info = self.db["sg_express"].find_one({"_id": ObjectId(express_id)})
        files = open(fpkm_path, 'r+')
        samples = files.readline().strip().split("\t")  # 此处已经有samples了
        files.close()
        fpkm = pd.read_table(fpkm_path)
        box = {}
        log2box = {}
        log10box = {}
        gene_list = {}
        log2gene_list = {}
        log10gene_list = {}

        box, gene_list = box_info(fpkm=fpkm, samples=samples)
        log2box, log2gene_list = box_info(fpkm=fpkm, log=2, samples=samples)
        log10box, log10gene_list = box_info(fpkm=fpkm, log=10, samples=samples)

        for sam in samples:
            data_log2 = [
                ("express_id", express_id),
                ("sample_group", sample_group),
                ("type", query_type),
                ("value_type", value_type)
            ]
            data_log2 += [
                ('{}_log2'.format(sam), log2box[sam])
            ]
            data_log2 += [
                ('{}_log10'.format(sam), log10box[sam])
            ]

            # data_log2=SON(data_log2)
            # log2_id = db['sg_express_box'].insert_one(data_log2).inserted_id  #每个样本的box值单独分开导表


            data_list_log2 = {}
            data_list_log10 = {}

            for keys, values in log2gene_list[sam].items():  # 每个样本的不同区段的gene_list分开导表
                # insert_data_log2 = [
                #     (keys,values),
                #     ("express_id",ObjectId(express_id)),
                #     ("box_id",ObjectId(log2_id))
                # ]
                # insert_data_log2 = SON(insert_data_log2)
                # data_list_log2.append(insert_data_log2)
                if values:
                    """如果values不为空"""
                    data_list_log2[keys] = ",".join(values)
                else:
                    data_list_log2[keys] = values
            for keys1, values1 in log10gene_list[sam].items():
                if values1:
                    """如果values1不为空"""
                    data_list_log10[keys1] = ",".join(values1)
                else:
                    data_list_log10[keys1] = values1
            data_log2 += [
                ('{}_log2_genelist'.format(sam), data_list_log2),
                ('{}_log10_genelist'.format(sam), data_list_log10)
            ]
            try:
                data_log2 = SON(data_log2)
                collection = self.db["sg_express_box"]
                box_id = collection.insert_one(data_log2).inserted_id
                print str(box_id)
            except Exception, e:
                self.bind_object.logger.error("导入盒形图错误：%s信息出错:%s" % (sam, e))
            else:
                self.bind_object.logger.info("导入盒形图: %s信息成功!" % sam)

    # @report_check
    def add_express_specimen_detail(self, express_id, rsem_result, rsem_type, sample=None):
        if not isinstance(express_id, ObjectId):
            if isinstance(express_id, types.StringTypes):
                express_id = ObjectId(express_id)
            else:
                raise Exception('express_id必须为ObjectId对象或其对应的字符串！')
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
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
            collection = self.db["sg_express_specimen_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入单样本表达量矩阵：%s信息出错:%s" % (rsem_result, e))
        else:
            self.bind_object.logger.info("导入单样本表达量矩阵: %s信息成功!" % rsem_result)

    # @report_check
    def get_diff_list(self, up_down_output, up_down=None):
        """
        :param diff_output_path, 差异分析生成的output文件夹，查找结尾是 'results_count'的文件
        :fc 筛选出fa的差异基因/转录本，导基因集的数据也是由此标准筛选之后生成的
        """
        import math
        if not os.path.exists(up_down_output):
            raise Exception("{}文件不存在，无法对up和down差异基因进行分类！".format(up_down_output))
        with open(up_down_output, 'r+') as f1:
            header = f1.readline()
            sequence = []

            for lines in f1:
                line = lines.strip().split("\t")
                seq_id = line[0]
                # print seq_id
                # print line[-2]
                regulate = line[-2]
                significant = line[-3]  # 显著性
                if significant == 'yes':
                    m_ = re.search(regulate, up_down)
                    if m_:
                        sequence.append(seq_id)
                    else:
                        pass
            if sequence:
                return sequence
            else:
                return None

    # @report_check

    def add_geneset(self, diff_stat_path, name=None, compare_name=None, ref_new=None, group_id=None,
                            express_method=None, type=None,
                            up_down=None, major=False):
        """
        添加sg_geneset主表, geneset的名字包括 up 和 down
        :param type: gene or transcript
        """
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(group_id, ObjectId):
            if isinstance(group_id, types.StringTypes):
                group_id = ObjectId(group_id)
            else:
                raise Exception('group_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(diff_stat_path):
            raise Exception('diff_stat_path所指的路径:{}不存在'.format(diff_stat_path))
        data_list_up = []
        data_list = []
        data_up = {
            'group_id': group_id,
            'task_id': task_id,
            'desc': '%s_vs_%s_差异基因集' % (name, compare_name),
            'project_sn': project_sn,
            'type': type
        }
        if not ref_new:
            raise Exception("请设置ref_new参数！")
        up_data = self.get_diff_list(up_down_output=diff_stat_path, up_down=up_down)
        if up_data:
            geneset_length = len(up_data)
            if geneset_length > 0:
                data_up['gene_length'] = int(geneset_length)
                if up_down == 'up' or up_down == 'down':
                    if type == 'gene':
                        data_up["name"] = '{}_vs_{}_{}_G_{}'.format(name, compare_name, up_down, ref_new)
                    if type == 'transcript':
                        data_up['name'] = '{}_vs_{}_{}_T_{}'.format(name, compare_name, up_down, ref_new)
                if up_down == 'up_down':
                    if type == 'gene':
                        data_up['name'] = '{}_vs_{}_G_{}'.format(name, compare_name, ref_new)
                    if type == 'transcript':
                        data_up['name'] = '{}_vs_{}_T_{}'.format(name, compare_name, ref_new)
                try:
                    collection = self.db["sg_geneset"]
                    print collection
                    geneset_up_id = collection.insert_one(data_up).inserted_id
                    print geneset_up_id
                    if major:
                        print("准备开始导入geneset_detail表！")
                        if geneset_up_id:
                            self.add_geneset_detail(geneset_up_id, diff_stat_path, up_data=up_data,
                                                    up_down=up_down)  # 直接导入detail表
                            print("导入geneset_detail表成功！")
                except Exception, e:
                    self.bind_object.logger.error("导入基因表达基因集：%s信息出错:%s" % (diff_stat_path, e))
                else:
                    self.bind_object.logger.info("导入基因表达基因集：%s信息成功!" % diff_stat_path)
                    return geneset_up_id
            else:
                print "{}对应{}调控基因集为空！".format(diff_stat_path, up_down)
        else:
            return None

    # @report_check
    def add_geneset_detail(self, geneset_id, diff_stat_path, up_data=None, up_down=None):
        """
        添加sg_geneset_detail表
        """
        if not isinstance(geneset_id, ObjectId):
            if isinstance(geneset_id, types.StringTypes):
                express_id = ObjectId(geneset_id)
            else:
                raise Exception('geneset_id必须为ObjectId对象或其对应的字符串！')
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        # geneset_id = str(geneset_id)
        data_list = []
        data_list_up = []
        if not up_data:
            up_data = self.get_diff_list(up_down_output=diff_stat_path, up_down=up_down)

        data = [
            ("geneset_id", ObjectId(geneset_id)),
            ("gene_list", up_data)
        ]
        data = SON(data)
        try:
            collection = self.db["sg_geneset_detail"]
            collection.insert_one(data)
        except Exception, e:
            self.bind_object.logger.error("导入基因集detail表: %s信息出错:%s" % (diff_stat_path, e))
        else:
            self.bind_object.logger.info("导入基因集detail表：%s信息成功!" % (diff_stat_path))

    # @report_check
    def add_express_diff(self, params, samples, compare_column, compare_column_specimen=None, ref_all=None,
                                workflow=True,is_duplicate=None, value_type="fpkm", express_method=None, diff_exp_dir=None,
                                class_code=None,query_type=None, express_id=None, name=None, group_id=None, group_detail=None,
                                control_id=None,major=True,pvalue_padjust = None,diff_method = None):

        # group_id, group_detail, control_id只供denovobase初始化时更新param使用
        """
        差异分析主表
        """
        if not isinstance(express_id, ObjectId):
            if isinstance(express_id, types.StringTypes):
                express_id = ObjectId(express_id)
            else:
                raise Exception('express_id必须为ObjectId对象或其对应的字符串！')
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        if params:
            params["submit_location"] = "expressdiff"
            params['task_id'] = task_id
            params['task_type'] = ''
        params['group_detail'] = group_detail_sort(params['group_detail'])
        if not express_method:
            raise Exception("add_express_diff函数需要设置express_method(选择表达量计算软件rsem或featurecounts)参数!")
        if not value_type:
            raise Exception("add_express_diff函数需要设置value_type(选择表达量水平fpkm或tpm)参数!")
        if "type" in params.keys() and "diff_method" in params.keys():
            diff_method = params['diff_method'].lower()
            re_name_info = {"gene": "G", "transcript": "T", "edger": "ER", "deseq2": "DS", "degseq": "DG",
                            "featurecounts": "FeaCount", "rsem": "RSEM"}
            re_name = 'DiffExp_{}_{}_{}_{}_'.format(re_name_info[query_type],
                                                    re_name_info[express_method.lower()], value_type.lower(),
                                                    re_name_info[diff_method]) + str(
                datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        else:
            raise Exception("params是字典格式，需要分别设置type和diff_method键值对!")
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else re_name,
            'desc': '表达量差异检测主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'specimen': sorted(samples),
            'status': 'end',
            'compare_column': compare_column,
            'group_detail': group_detail,
            'express_id': express_id,
            "is_duplicate": is_duplicate,
            "value_type": value_type
        }
        if params['type'] == 'gene':
            insert_data["genes"] = True
            insert_data["trans"] = False
        elif params['type'] == 'transcript':
            insert_data["genes"] = False
            insert_data["trans"] = True
        if compare_column_specimen:
            insert_data["compare_column_specimen"] = compare_column_specimen
        # if group_id == 'all':
        #     insert_data['group_detail'] = {'all': group_detail}
        collection = self.db['sg_express_diff']
        express_diff_id = collection.insert_one(insert_data).inserted_id
        if major:
            diff_exp_files = os.listdir(diff_exp_dir)
            for f in diff_exp_files:
                if re.search(r'_edgr_stat.xls$', f):
                    con_exp = f.split('_edgr_stat.xls')[0].split('_vs_')
                    name = con_exp[0]
                    compare_name = con_exp[1]
                    self.add_express_diff_detail(express_diff_id, name, compare_name, ref_all,
                                                 os.path.join(diff_exp_dir, f), workflow,
                                                 class_code, query_type, params["pvalue_padjust"])
        return express_diff_id

    def add_express_diff_detail(self, express_diff_id, name, compare_name, ref_all, diff_stat_path, workflow=False,
                            class_code=None, query_type=None, pvalue_padjust=None):
        """
        group:为两两比较的样本或分组名，列表
        query_type: gene/transcript
        diff_stat_path: 差异统计表
        workflow: 是否工作流导入文件
        ref_all: ref 还是 ref+new
        """
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        if not isinstance(express_diff_id, ObjectId):
            if isinstance(express_diff_id, types.StringTypes):
                express_diff_id = ObjectId(express_diff_id)
            else:
                raise Exception('express_diff_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(diff_stat_path):
            raise Exception('diff_stat_path所指定的路径:{}不存在，请检查！'.format(diff_stat_path))
        if class_code:
            if os.path.exists(class_code):
                name_seq_id = self.get_gene_name(class_code, query_type=query_type, workflow=workflow)
            else:
                raise Exception("{} not exist".format(class_code))
        if pvalue_padjust not in ['padjust', 'pvalue']:
            raise ValueError('pvalue_padjust must be padjust or pvalue')

        # dump data of diff_stat_path into mongodb
        diff_table = pd.read_table(diff_stat_path)
        # begin to filter repeated columns
        new_columns = list()
        target_columns = list()
        for col in diff_table.columns:
            judge1 = col in target_columns
            judge2 = col.endswith('.1') or col.endswith('.2') or col.endswith('.3')
            if (not judge1) and (not judge2):
                target_columns.append(col)
                new_columns.append(col)
            else:
                new_columns.append('not_need')
        diff_table.columns = new_columns
        diff_table = diff_table.loc[:, target_columns]
        # end of filtering repeated columns
        sig_status = list()
        sig_mark = diff_table['significant']
        reg_list = diff_table['regulate']
        if 'no' in list(sig_mark):
            sig_status.append('nosig')
        if 'yes' in list(sig_mark):
            if 'down' in list(reg_list[sig_mark == 'yes']):
                sig_status.append('down')
            if 'up' in list(reg_list[sig_mark == 'yes']):
                sig_status.append('up')
        sig_pvalues = diff_table[pvalue_padjust][diff_table['significant'] == "yes"]
        log10_pvalue_list = sorted([-log10(x) for x in sig_pvalues if x > 0])

        if len(log10_pvalue_list) == 0:
            log10_pvalue_cutoff = 320
        else:
            if len(sig_pvalues) > 2000:
                log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.85)]
            elif len(sig_pvalues) > 1000:
                log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.90)]
            elif len(sig_pvalues) > 500:
                log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.95)]
            elif len(sig_pvalues) > 250:
                log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.99)]
            elif len(sig_pvalues) == 0:
                tmp_list = sorted([-log10(x) for x in diff_table[pvalue_padjust] if x > 0])
                if len(tmp_list) == 0:
                    log10_pvalue_cutoff = 320
                else:
                    log10_pvalue_cutoff = tmp_list[int(len(tmp_list)*0.9)]
            else:
                log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.8)]

        self.bind_object.logger.info("log10_pvalue_cutoff: {}".format(log10_pvalue_cutoff))

        target_dict_list = list()
        marker_info = dict(name=name, ref_all=ref_all, compare_name=compare_name,
                           express_diff_id=express_diff_id)

        row_dict_list = diff_table.to_dict('records')
        for row_dict in row_dict_list:
            seq_id = row_dict['seq_id']
            if query_type == "gene":
                gene_name = name_seq_id[seq_id]['gene_name']
                gene_id = seq_id
            else:
                gene_name = name_seq_id[seq_id]["gene_name"]
                gene_id = name_seq_id[seq_id]["gene_id"]
            row_dict.update(dict(gene_name=gene_name, gene_id=gene_id))

            log2fc = row_dict['log2fc']
            row_dict['fc'] = round(2**log2fc, 3)

            pvalue = row_dict[pvalue_padjust]
            if pvalue <= 0:
                pvalue = 1e-320
            try:
                if -log10(pvalue) > log10_pvalue_cutoff:
                    log10_pvalue = log10_pvalue_cutoff
                else:
                    log10_pvalue = -log10(pvalue)
            except:
                print(type(pvalue))
                raise Exception('pvalue is {}'.format(pvalue))
            row_dict['log10_'+pvalue_padjust] = log10_pvalue

            row_dict.update(marker_info)

            target_dict_list.append(row_dict)

        try:
            collection = self.db["sg_express_diff_detail"]
            collection.insert_many(target_dict_list)
            con = self.db["sg_express_diff"]
            sig_status_name = name+'_vs_'+compare_name+'_'+ref_all+'_status'
            con.update({'_id': express_diff_id}, {"$set": {sig_status_name: sig_status}})
        except Exception, e:
            self.bind_object.logger.error("导入基因表达差异统计表：%s出错:%s" % (diff_stat_path, e))
        else:
            self.bind_object.logger.info("导入基因表达差异统计表：%s信息成功!" % diff_stat_path)

    def add_diff_summary_detail(self, diff_express_id, count_path, ref_all, query_type=None,
                                class_code=None, workflow=False):
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        if not isinstance(diff_express_id, ObjectId):
            if isinstance(diff_express_id, types.StringTypes):
                diff_express_id = ObjectId(diff_express_id)
            else:
                raise Exception('express_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(count_path):
            raise Exception('count_path:{}所指路径不存在'.format(count_path))
        if class_code:
            if os.path.exists(class_code):
                name_seq_id = self.get_gene_name(class_code, query_type, workflow=workflow)
        data_list = list()
        with open(count_path, 'rb') as f:
            sample = f.readline().strip().split('\t')
            lensam = len(sample)
            sample = sample[1:lensam]
            for line in f:
                l = line.strip().split('\t')
                gene_id = l[0]
                alen = len(l)
                blen = alen - 2
                alen = alen - 1
                fpkm = l[1:alen]
                sum_1 = l[alen]
                data = [
                    ("seq_id", gene_id),
                    ("express_diff_id", diff_express_id),
                    ('sum', int(sum_1)),
                    ('ref_all', ref_all)
                ]
                if class_code:
                    if name_seq_id:
                        if query_type == 'transcript':
                            if gene_id in name_seq_id.keys():
                                true_gene_id = name_seq_id[gene_id]['gene_id']
                                data.append(("gene_id", true_gene_id))
                                data.append(('gene_name', name_seq_id[gene_id]['gene_name']))
                        else:
                            data.append(('gene_name', name_seq_id[gene_id]['gene_name']))
                for j in range(blen):
                    data += [('{}_diff'.format(sample[j]), fpkm[j]), ]

                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_express_diff_summary"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入差异分析summary表出错:%s" % e)
        else:
            self.bind_object.logger.info("导入差异分析summary表成功!")


if __name__ == "__main__":
    ####################################################################################################################################
    ##################----------------------rsem fpkm 导表
    ####################################################################################################################################
    # rsem_dir = "/mnt/ilustre/users/sanger-dev/workspace/20170630/Single_rsem_stringtie_mouse_fpkm_2_diff_stat_2/Expresstest3/output/rsem"
    # rsem_dir = '/mnt/ilustre/users/sanger-dev/workspace/20170707/Single_rsem_stringtie_mouse_total_2/Expresstest2/MergeRsem1/output'
    # rsem_dir = '/mnt/ilustre/users/sanger-dev/workspace/20170707/Single_rsem_stringtie_mouse_total_2/Expresstest2/MergeRsem/output'
    # rsem_dir = '/mnt/ilustre/users/sanger-dev/workspace/20170702/Single_rsem_stringtie_mouse_total_1/Express/MergeRsem1/output'
    # rsem_dir = "/mnt/ilustre/users/sanger-dev/workspace/20170707/Single_rsem_stringtie_mouse_total_2/Expresstest2/MergeRsem1/output"
    # is_duplicate = True
    # # class_code = '/mnt/ilustre/users/sanger-dev/workspace/20170707/Single_rsem_stringtie_mouse_total_2/Expresstest2/MergeRsem/class_code'
    # # class_code = "/mnt/ilustre/users/sanger-dev/workspace/20170702/Single_rsem_stringtie_mouse_total_1/Express/MergeRsem1/class_code"
    # class_code = '/mnt/ilustre/users/sanger-dev/workspace/20170707/Single_rsem_stringtie_mouse_total_2/Expresstest2/MergeRsem1/class_code'
    # sample = "samples"
    # distri_path = '/mnt/ilustre/users/sanger-dev/workspace/20170707/Single_rsem_stringtie_mouse_total_2/Expresstest2/MergeRsem1/'
    # # distri_path = '/mnt/ilustre/users/sanger-dev/workspace/20170707/Single_rsem_stringtie_mouse_total_2/Expresstest2/MergeRsem'
    # # distri_path = "/mnt/ilustre/users/sanger-dev/workspace/20170707/Single_rsem_stringtie_mouse_total_2/Expresstest2/MergeRsem1"
    # # distri_path = '/mnt/ilustre/users/sanger-dev/workspace/20170702/Single_rsem_stringtie_mouse_total_1/Express/MergeRsem1/'
    # # distri_path = '/mnt/ilustre/users/sanger-dev/workspace/20170630/Single_rsem_stringtie_mouse_fpkm_2_diff_stat_2/Expresstest3/MergeRsem'
    # params = {}
    # params["group_id"] = "596452d7edcb255322d9e66e"
    # params['group_detail'] = {
    #     "A":['596452d7edcb255322d9dbe1','596452d7edcb255322d9dbdf','596452d7edcb255322d9dbe0'],
    #     "C":['596452d7edcb255322d9dbdd','596452d7edcb255322d9dbde','596452d7edcb255322d9dbdc'],
    #     "B":['596452d7edcb255322d9dbda','596452d7edcb255322d9dbdb','596452d7edcb255322d9dbd9']
    # }
    # # params['group_id'] = '5955f5e1edcb253a204f8988'
    # # params["group_detail"] = {
    # #     "X1": ["5955f5deedcb253a204f7ef5", "5955f5deedcb253a204f7ef4", "5955f5deedcb253a204f7ef3"],
    # #     "B1": ["5955f5deedcb253a204f7efb", "5955f5deedcb253a204f7ef9", "5955f5deedcb253a204f7efa"],
    # #     "Z1": ["5955f5deedcb253a204f7ef7", "5955f5deedcb253a204f7ef6", "5955f5deedcb253a204f7ef8"]}
    # params["express_method"] = "rsem"
    # params["type"] = "fpkm"
    # group_fpkm_path = distri_path + "/group"
    # samples = ['A_1', 'A_2', 'A_3', 'B_1', 'B_2', 'B_3', 'C_1', 'C_2', 'C_3']
    # # samples = ["X1_1","X1_2","X1_3","B1_1","B1_2","B1_3","Z1_1","Z1_2","Z1_3"]
    # a = RefrnaExpress2()
    # a.add_express(rsem_dir=rsem_dir, is_duplicate=True, group_fpkm_path=group_fpkm_path, samples=samples,
    #               class_code=class_code,
    #               params=params, major=True, distri_path=distri_path)
    # print 'end!'

    #####################################################################################################################################
    #################------------------------featurecounts 导表
    #####################################################################################################################################
    feature_dir = '/mnt/ilustre/users/sanger-dev/workspace/20170706/Single_feature_stringtie_mouse_1/Express/output/featurecounts'
    distri_path = '/mnt/ilustre/users/sanger-dev/workspace/20170706/Single_feature_stringtie_mouse_1/Express/Featurecounts'
    # feature_dir = "/mnt/ilustre/users/sanger-dev/workspace/20170629/Single_feature_stringtie_mouse_3/Express/output/featurecounts"
    # feature_dir = "/mnt/ilustre/users/sanger-dev/workspace/20170706/Single_feature_stringtie_mouse_1/Express/output/featurecounts"
    # distri_path = "/mnt/ilustre/users/sanger-dev/workspace/20170706/Single_feature_stringtie_mouse_1/Express/Featurecounts/"
    group_fpkm_path = distri_path + "/group"
    # group_fpkm_path = '/mnt/ilustre/users/sanger-dev/workspace/20170629/Single_feature_stringtie_mouse_3/Express/Featurecounts/group'
    is_duplicate = True
    samples = ['A_1', 'A_2', 'A_3', 'B_1', 'B_2', 'B_3', 'C_1', 'C_2', 'C_3']
    params = {}

    params["group_id"] = "596452d7edcb255322d9e66e"
    params['group_detail'] = {
        "A": ['596452d7edcb255322d9dbe1', '596452d7edcb255322d9dbdf', '596452d7edcb255322d9dbe0'],
        "C": ['596452d7edcb255322d9dbdd', '596452d7edcb255322d9dbde', '596452d7edcb255322d9dbdc'],
        "B": ['596452d7edcb255322d9dbda', '596452d7edcb255322d9dbdb', '596452d7edcb255322d9dbd9']
    }
    class_code = "/mnt/ilustre/users/sanger-dev/workspace/20170707/Single_rsem_stringtie_mouse_total_2/Expresstest2/MergeRsem/class_code"
    params["type"] = "fpkm"
    params["express_method"] = "featurecounts"

    a = RefrnaExpress2()
    a.add_express_feature(feature_dir=feature_dir, group_fpkm_path=group_fpkm_path, class_code=class_code,
                          is_duplicate=True, samples=samples,
                          params=params, major=True, distri_path=distri_path)
    print 'end!'
    ####################################################################################################


    ####################################################################################################################################
    ######################-------------gene set 导表 ref_new的参数有两种选择  ref和refandnew
    ####################################################################################################################################
    # a = RefrnaExpress2()
    # group_id = '596452d7edcb255322d9e66e'
    # # group_id = '5955f5e1edcb253a204f8988'
    # # path = "/mnt/ilustre/users/sanger-dev/workspace/20170527/DiffExpress_deno222_7450_1534/DiffExp/output"
    # # path = '/mnt/ilustre/users/sanger-dev/workspace/20170630/Single_rsem_stringtie_mouse_fpkm_2_diff_stat_2/Expresstest3/output/diff/genes_diff/diff_stat_dir'
    # # path = '/mnt/ilustre/users/sanger-dev/workspace/20170630/Single_rsem_stringtie_mouse_fpkm_2_diff_stat_2/Expresstest3/output/ref_diff/genes_ref_diff/diff_stat_dir'
    # # path = '/mnt/ilustre/users/sanger-dev/workspace/20170630/Single_rsem_stringtie_mouse_fpkm_2_diff_stat_2/Expresstest3/output/ref_diff/trans_ref_diff/diff_stat_dir'
    # # path = '/mnt/ilustre/users/sanger-dev/workspace/20170707/Single_rsem_stringtie_mouse_total_2/Expresstest2/output/diff/trans_diff/diff_stat_dir'
    # path = "/mnt/ilustre/users/sanger-dev/workspace/20170707/Single_rsem_stringtie_mouse_total_2/Expresstest2/output/ref_diff/genes_ref_diff/diff_stat_dir"
    # # path = '/mnt/ilustre/users/sanger-dev/workspace/20170630/Single_rsem_stringtie_mouse_fpkm_2_diff_stat_2/Expresstest3/output/diff/trans_diff/diff_stat_dir'
    # for files in os.listdir(path):
    #     if re.search(r'edgr_stat.xls', files):
    #         print files
    #         m_ = re.search(r'(\w+?)_vs_(\w+?).edgr_stat.xls', files)
    #         if m_:
    #             name = m_.group(1)
    #             compare_name = m_.group(2)
    #             up_down = a.add_geneset(diff_stat_path=path + "/" + files, group_id=group_id, name=name,
    #                                     compare_name=compare_name, ref_new='ref', express_method="rsem",
    #                                     type="gene", major=True, up_down='up_down')
    #             down_id = a.add_geneset(diff_stat_path=path + "/" + files, group_id=group_id, name=name, major=True,
    #                                     compare_name=compare_name, ref_new="ref", express_method="rsem",
    #                                     type="gene", up_down='down')
    #             up_id = a.add_geneset(diff_stat_path=path + "/" + files, group_id=group_id, name=name, major=True,
    #                                   compare_name=compare_name, ref_new='ref', express_method="rsem",
    #                                   type="gene", up_down='up')
    #
    #             print up_id
    #             print down_id
    #             print name, compare_name
    #         print 'end'

    ####################################################################################################################################
    ######################-------------差异分析 导表
    ####################################################################################################################################

    # path = "/mnt/ilustre/users/sanger-dev/workspace/20170524/Single_rsem_stringtie_fpkm_5/Express/output/diff/genes_diff/diff_stat_dir"
    # path = '/mnt/ilustre/users/sanger-dev/workspace/20170630/Single_rsem_stringtie_mouse_fpkm_2_diff_stat_2/Expresstest3/output/diff/genes_diff/diff_stat_dir'
    # path = '/mnt/ilustre/users/sanger-dev/workspace/20170630/Single_rsem_stringtie_mouse_fpkm_2_diff_stat_2/Expresstest3/output/diff/trans_diff/diff_stat_dir'
    # path = '/mnt/ilustre/users/sanger-dev/workspace/20170701/Single_rsem_stringtie_mouse_fpkm_2_diff_stat_2/Expresstest3/output/diff/trans_diff/diff_stat_dir'
    # is_duplicate = True
    # sample = ['X1_1', 'X1_2', 'X1_3', 'B1_1', 'B1_2', 'B1_3', 'Z1_1', 'Z1_2', 'Z1_3']
    # compare_column = ["X1|B1", "X1|Z1", "B1|Z1"]
    # params = {}
    # class_code = '/mnt/ilustre/users/sanger-dev/workspace/20170701/Single_rsem_stringtie_mouse_fpkm_2_diff_stat_2/Expresstest3/MergeRsem/class_code'
    # merge_path = '/mnt/ilustre/users/sanger-dev/workspace/20170701/Single_rsem_stringtie_mouse_fpkm_2_diff_stat_2/Expresstest3/output/diff/trans_diff/merge.xls'
    # params['control_id'] = '5955f821f2e3f7fddea08f6e'
    # params["group_id"] = "5955f5e1edcb253a204f8988"
    # params["group_detail"] = {
    #     "X1": ["5955f5deedcb253a204f7ef5", "5955f5deedcb253a204f7ef4", "5955f5deedcb253a204f7ef3"],
    #     "B1": ["5955f5deedcb253a204f7efb", "5955f5deedcb253a204f7ef9", "5955f5deedcb253a204f7efa"],
    #     "Z1": ["5955f5deedcb253a204f7ef7", "5955f5deedcb253a204f7ef6", "5955f5deedcb253a204f7ef8"]}
    # params['express_id'] = str("59560618a4e1af1ae5309bf3")
    # params['fc'] = 0.1
    # params['pvalue_padjust'] = 'padjust'
    # params['pvalue'] = 0.5
    # params['type'] = 'transcript'
    # params['diff_method'] = 'DESeq2'
    # compare_column_specimen = {"X1|B1": ['X1_1', 'X1_2', 'X1_3', 'B1_1', 'B1_2', 'B1_3'],
    #                            "X1|Z1": ['X1_1', 'X1_2', 'X1_3', 'Z1_1', 'Z1_2', 'Z1_3'],
    #                            "B1|Z1": ['B1_1', 'B1_2', 'B1_3', 'Z1_1', 'Z1_2', 'Z1_3']}
    # a = RefrnaExpress2()
    # diff_express_id = a.add_express_diff(params=params, samples=sample, compare_column=compare_column, ref_all='all',
    #                                      class_code=class_code, is_duplicate=True,
    #                                      diff_exp_dir=path, query_type="transcript",
    #                                      express_id=ObjectId("592e26fba4e1af397f263b38"),
    #                                      compare_column_specimen=compare_column_specimen,
    #                                      major=True, group_id=params["group_id"], workflow=True,
    #                                      pvalue_padjust='padjust',express_method = 'rsem',value_type='fpkm')
    # a.add_diff_summary_detail(diff_express_id, merge_path, ref_all='all', query_type='transcript',
    #                           class_code=class_code, workflow=True)
    # print 'end'
    # print "diff_express_id:ObjectId({})".format(str(diff_express_id))

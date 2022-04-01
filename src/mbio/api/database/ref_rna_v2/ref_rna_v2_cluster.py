# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# modify by khl

from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import json


class RefRnvV2Cluster(Base):
    def __init__(self, bind_object):
        super(RefRnvV2Cluster, self).__init__(bind_object)
        #self._db_name = Config().MONGODB + '_rna'
        #self._db_ref = Confif().MONGODB + "_ref_rna"
        self.denovo_db = Config().mongo_client[Config().MONGODB + "_rna"]
        # self.ref_db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        self._project_type = 'ref_rna_v2'  # modified by hongdong 20171102
        # self.ref_db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]

    def get_gene_name(self, class_code, query_type=None, workflow=False):
        """
        :params: 是否工作流根据class_code信息导入基因/转录本名称
        对转录本都加上相应的gene id信息
        """
        with open(class_code, 'r+') as f1:
            f1.readline()
            data = {}
            for lines in f1:
                line = lines.strip().split("\t")
                if workflow:
                    if query_type == 'transcript':
                        data[line[0]] = {"gene_name": line[3], "class_code": line[2], "gene_id": line[1]}
                    if query_type == 'gene':
                        data[line[1]] = {"gene_name": line[3], "class_code": line[2]}
                else:
                    if query_type == 'gene':
                        data[line[0]] = {"gene_name": line[1]}
                    if query_type == 'transcript':
                        data[line[0]] = {"gene_name": line[1], 'gene_id': line[3]}
            return data

    #@report_check
    def add_cluster(self, params, express_id, sample_tree=None, gene_tree=None, name=None, samples=None, genes=None, project='ref'):
        if not isinstance(express_id, ObjectId):
            if isinstance(express_id, types.StringTypes):
                express_id = ObjectId(express_id)
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        if gene_tree:
            with open(gene_tree, 'rb') as g:
                gene_tree = g.readlines()[0].strip('\n')
        if sample_tree:
            with open(sample_tree, 'rb') as s:
                sample_tree = s.readlines()[0].strip('\n')
        params['diff_fpkm'] = str(express_id)
        if 'express_method' in params.keys() and 'level' in params.keys() and 'type' in params.keys() and 'log' in params.keys():
            express_method = params['express_method']
            express_level = params['express_level']
            query_type = params['query_type']
            log_info = params['log']
            re_name = 'GSetCluster_{}_{}_{}_{}_'.format(express_method, express_level, query_type, log_info) + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        else:
            raise Exception("params需要有express_method、level、type、log字段!")
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else re_name,
            'desc': '基因集聚类分析主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'specimen': samples,
            'status': 'end',
            'express_id': express_id,
            'sample_tree': sample_tree,
            'gene_tree': gene_tree,
            'genes': genes,
        }
        if project == 'denovo':
            collection = self.denovo_db['sg_denovo_cluster']
        elif project == 'ref':
            collection = self.db["sg_geneset_cluster"]
        else:
            raise Exception("请设置project参数为denovo或ref")
        cluster_id = collection.insert_one(insert_data).inserted_id
        return cluster_id

    #@report_check
    def add_cluster_detail(self, cluster_id, sub, sub_path, class_code=None,project="ref",workflow=True,query_type=None):
        if not isinstance(cluster_id, ObjectId):
            if isinstance(cluster_id, types.StringTypes):
                cluster_id = ObjectId(cluster_id)
            else:
                raise Exception('cluster_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(sub_path):
            raise Exception('sub_path所指定的路径:{}不存在，请检查！'.format(sub_path))
        data_list = []
        # if query_type == 'transcript':
        #     if not os.path.exists(class_code):
        #         raise Exception("转录本的聚类需要输入class_code信息!")
        #     else:
        class_code_info = self.get_gene_name(class_code,query_type=query_type,workflow=workflow)
        with open(sub_path, 'rb') as f:
            head = f.readline().strip().split('\t')
            for line in f:
                line = line.strip().split('\t')
                data = [
                    ('sub_cluster', int(sub)),
                    ('cluster_id', cluster_id),
                    ('seq_id', line[0]),
                    ('gene_name', class_code_info[line[0]]['gene_name'])
                ]
                if query_type == 'transcript':
                    if line[0] in class_code_info.keys():
                        seq_id = class_code_info[line[0]]['gene_id']
                    else:
                        seq_id = '-'
                    data.append(('gene_id', seq_id))

                for i in range(len(head)):
                    data.append((head[i], round(float(line[i + 1]), 4)))
                data = SON(data)
                data_list.append(data)
        try:
            if project == "denovo":
                collection = self.denovo_db["sg_denovo_cluster_detail"]
            elif project == "ref":
                collection = self.db["sg_geneset_cluster_detail"]
            else:
                raise Exception("请设置project参数为denovo或ref!")
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入子聚类统计表：%s信息出错:%s" % (sub_path, e))
        else:
            self.bind_object.logger.info("导入子聚类统计表:%s信息成功!" % sub_path)

# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
import os
import re
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
import gridfs
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config


class AnnotationStat(Base):
    def __init__(self, bind_object):
        super(AnnotationStat, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + '_ref_rna'

    @report_check
    def add_annotation_stat(self, name=None, params=None, seq_type=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationStat_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'status': 'end',
            'desc': '注释统计主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'seq_type': seq_type
        }
        collection = self.db['sg_annotation_stat']
        stat_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_stat!")
        return stat_id

    @report_check
    def add_annotation_stat_detail(self, stat_id, stat_path, venn_path):
        """
        database: 进行统计的数据库
        stat_path: all_annotation_statistics.xls
        venn_path: venn图目录
        """
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(stat_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(stat_path))
        if not os.path.exists(venn_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(venn_path))
        data_list = []
        with open(stat_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('stat_id', stat_id),
                    ('type', line[0]),
                    ('transcript', int(line[1])),
                    ('gene', int(line[2])),
                    ('transcript_percent', round(float(line[3]), 4)),
                    ('gene_percent', round(float(line[4]), 4)),
                ]
                venn_list, gene_venn_list = None, None
                database = ["nr", "swissprot", "pfam", "kegg", "go", "string", "cog"]
                if line[0] in database:
                    venn = venn_path + "/" + line[0] + "_venn.txt"
                    gene_venn = venn_path + "/gene_" + line[0] + "_venn.txt"
                    if os.path.exists(venn) and os.path.exists(gene_venn):
                        with open(venn, "rb") as f:
                            venn_list = f.readline().strip('\n')
                            for line in f:
                                venn_list += ',{}'.format(line.strip('\n'))
                        with open(gene_venn, "rb") as f:
                            gene_venn_list = f.readline().strip('\n')
                            for line in f:
                                gene_venn_list += ',{}'.format(line.strip('\n'))
                        data.append(("gene_list", venn_list))
                        data.append(("transcript_list", gene_venn_list))
                    else:
                        raise Exception("{}对应的venn.txt文件不存在".format(line[0]))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_stat_detail']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入注释统计信息：%s出错!" % (stat_path))
        else:
            self.bind_object.logger.info("导入注释统计信息：%s成功!" % (stat_path))

    @report_check
    def add_stat_detail(self, old_stat_id, stat_id, nr_blast, gene_nr_blast, sw_blast, gene_sw_blast):
        """
        注释重运行时注释统计导表sg_annotation_stat_detail
        """
        if not isinstance(old_stat_id, ObjectId):
            if isinstance(old_stat_id, types.StringTypes):
                old_stat_id = ObjectId(old_stat_id)
            else:
                raise Exception('old_stat_id必须为ObjectId对象或其对应的字符串！')
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        collection = self.db["sg_annotation_stat_detail"]
        results = collection.find({"stat_id": old_stat_id})
        data_list, data = [], []
        total_anno_tran, total_anno_gene = [], []
        for result in results:
            db = result["type"]
            if db == "total":
                total_tran = result["transcript"]
                total_gene = result["gene"]
            data = [
                ('stat_id', stat_id),
                ('type', result["type"]),
                ('transcript', result["transcript"]),
                ('gene', result["gene"]),
                ('transcript_percent', result["transcript_percent"]),
                ('gene_percent', result["gene_percent"])
            ]
            if db == "total":
                data = SON(data)
                data_list.append(data)
            if db == "pfam":
                data.append(('gene_list', result["gene_list"]))
                data.append(('transcript_list', result["transcript_list"]))
                for g in result["gene_list"].split(","):
                    total_anno_gene.append(g)
                for t in result["transcript_list"].split(","):
                    total_anno_tran.append(t)
                data = SON(data)
                data_list.append(data)
        nr_ids = self.stat(stat_path=nr_blast)
        gene_nr_ids = self.stat(stat_path=gene_nr_blast)
        for t in nr_ids:
            total_anno_tran.append(t)
        for g in gene_nr_ids:
            total_anno_gene.append(g)
        data = [
            ('stat_id', stat_id),
            ('type', "nr"),
            ('transcript', len(nr_ids)),
            ('gene', len(gene_nr_ids)),
            ('transcript_percent', round(float(len(nr_ids))/total_tran, 4)),
            ('gene_percent', round(float(len(gene_nr_ids))/total_gene, 4)),
            ('gene_list', ",".join(gene_nr_ids)),
            ('transcript_list', ",".join(nr_ids))
        ]
        data = SON(data)
        data_list.append(data)
        sw_ids = self.stat(stat_path=sw_blast)
        gene_sw_ids = self.stat(stat_path=gene_sw_blast)
        for t in sw_ids:
            total_anno_tran.append(t)
        for g in gene_sw_ids:
            total_anno_gene.append(g)
        data = [
            ('stat_id', stat_id),
            ('type', "swissprot"),
            ('transcript', len(sw_ids)),
            ('gene', len(gene_sw_ids)),
            ('transcript_percent', round(float(len(sw_ids))/total_tran, 4)),
            ('gene_percent', round(float(len(gene_sw_ids))/total_gene, 4)),
            ('gene_list', ",".join(gene_sw_ids)),
            ('transcript_list', ",".join(sw_ids))
        ]
        data = SON(data)
        data_list.append(data)
        total_anno_tran = list(set(total_anno_tran))
        total_anno_gene = list(set(total_anno_gene))
        data = [
            ('stat_id', stat_id),
            ('type', "total_anno"),
            ('transcript', len(total_anno_tran)),
            ('gene', len(total_anno_gene)),
            ('transcript_percent', round(float(len(total_anno_tran))/total_tran, 4)),
            ('gene_percent', round(float(len(total_anno_gene))/total_gene, 4)),
            ('gene_list', ",".join(total_anno_gene)),
            ('transcript_list', ",".join(total_anno_tran))
        ]
        data = SON(data)
        data_list.append(data)
        # 细节表里添加一条记录，与流程保持一致,表示NR SWISSPROT PFAM三个库注释的并集 刘彬旭
        data = [
            ('stat_id', stat_id),
            ('type', "total_anno_nsp"),
            ('transcript', len(total_anno_tran)),
            ('gene', len(total_anno_gene)),
            ('transcript_percent', round(float(len(total_anno_tran))/total_tran, 4)),
            ('gene_percent', round(float(len(total_anno_gene))/total_gene, 4)),
            ('gene_list', ",".join(total_anno_gene)),
            ('transcript_list', ",".join(total_anno_tran))
        ]
        data = SON(data)
        data_list.append(data)

        try:
            collection = self.db['sg_annotation_stat_detail']
            collection.insert_many(data_list)
        except:
            self.bind_object.logger.error("导入注释统计信息出错")
        else:
            self.bind_object.logger.info("导入注释统计信息成功")

    def stat(self, stat_path):
        with open(stat_path, "rb") as f:
            id_list = []
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                q_id = line[5]
                id_list.append(q_id)
        id_list = list(set(id_list))
        return id_list

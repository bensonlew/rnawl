# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import re
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
import gridfs
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import pymongo


class AnnotationBlast(Base):
    def __init__(self, bind_object):
        super(AnnotationBlast, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + '_ref_rna'

    @report_check
    def add_annotation_blast(self, name=None, params=None, stat_id=None):
        """
        blast_path: 无参工作流中annotation输出目录
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationBlast_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'status': 'end',
            'desc': 'blast最佳比对结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_blast']
        blast_id = collection.insert_one(insert_data).inserted_id
        print("add ref_annotation_blast!")
        return blast_id

    @report_check
    def add_annotation_blast_detail(self, blast_id, seq_type, anno_type, database, blast_path):
        if not isinstance(blast_id, ObjectId):
            if isinstance(blast_id, types.StringTypes):
                blast_id = ObjectId(blast_id)
            else:
                raise Exception('blast_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(blast_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(blast_path))
        data_list = []
        with open(blast_path, 'r') as f:
            lines = f.readlines()
            best, filte = [], []
            for i in range(0, len(lines)):
                line = lines[i].strip().split("\t")
                hit = line[10]
                query = line[5]
                index = i
                for j in range(i+1, len(lines)):
                    second = lines[j].strip().split("\t")
                    second_hit = second[10]
                    second_query = second[5]
                    if query == second_query:
                        if second_hit != hit:
                            index = j
                            break
                    else:
                        index = j
                        break
                best.append(index)
            for x in best:
                if x not in filte:
                    filte.append(x)
            for i in filte:
                line = lines[i].strip().split('\t')
                data = [
                    ('blast_id', blast_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('database', database),
                    ('score', float(line[0])),
                    ('e_value', line[1]),
                    ('hsp_len', line[2]),
                    ('identity_rate', round(float(line[3]), 4)),
                    ('similarity_rate', round(float(line[4]), 4)),
                    ('query_id', line[5]),
                    ('q_len', line[6]),
                    ('q_begin', line[7]),
                    ('q_end', line[8]),
                    ('q_frame', line[9]),
                    ('hit_name', line[10]),
                    ('hit_len', line[11]),
                    ('hsp_begin', line[12]),
                    ('hsp_end', line[13]),
                    ('hsp_frame', line[14]),
                    ('description', line[15])
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_blast_detail']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入注释统计信息：%s出错!" % (blast_path))
        else:
            self.bind_object.logger.info("导入注释统计信息：%s成功!" % (blast_path))

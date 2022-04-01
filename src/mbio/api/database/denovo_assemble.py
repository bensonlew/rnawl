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


class DenovoAssemble(Base):
    def __init__(self, bind_object):
        super(DenovoAssemble, self).__init__(bind_object)
        self._db_name = Config().MONGODB + '_rna'

    @report_check
    def add_sequence(self, trinity_path, gene_path):
        if not os.path.exists(trinity_path) or not os.path.exists(gene_path):
            raise Exception('trinity_path or gene_path所指定的路径不存在，请检查！')
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': 'denovo组装拼接主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'trinity_path': trinity_path,
            'gene_path': gene_path,
            'status': 'end',
        }
        collection = self.db['sg_denovo_sequence']
        sequence_id = collection.insert_one(insert_data).inserted_id
        return sequence_id

    @report_check
    def add_sequence_detail(self, sequence_id, stat_path):
        if not isinstance(sequence_id, ObjectId):
            if isinstance(sequence_id, types.StringTypes):
                sequence_id = ObjectId(sequence_id)
            else:
                raise Exception('sequence_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(stat_path):
            raise Exception('stat_path所指定的路径不存在，请检查！')
        data_list = list()
        with open(stat_path, 'rb') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('stastistic_name', line[0]),
                    ('transcripts', float(line[2])),
                    ('genes', float(line[1])),
                    ('sequence_id', sequence_id)
                ]
                if line[3] != '--':
                    data += [('genes_count', int(line[3])), ('transcripts_count', int(line[4]))]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_denovo_sequence_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (stat_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % stat_path)

    @report_check
    def add_sequence_step(self, sequence_id, length_path, step):
        if not os.path.exists(length_path):
            raise Exception('length_path所指定的路径不存在，请检查！')
        if not isinstance(sequence_id, ObjectId):
            if isinstance(sequence_id, types.StringTypes):
                sequence_id = ObjectId(sequence_id)
            else:
                raise Exception('sequence_id必须为ObjectId对象或其对应的字符串！')
        data_list = list()
        with open(length_path, 'rb') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('length', line[0]),
                    ('transcripts_num', int(line[3])),
                    ('genes_num', int(line[1])),
                    ('transcripts_per', round(float(line[4]), 4)),
                    ('genes_per', round(float(line[2]), 4)),
                    ('sequence_id', sequence_id),
                    ('step', int(step))
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_denovo_sequence_step"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (length_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % length_path)

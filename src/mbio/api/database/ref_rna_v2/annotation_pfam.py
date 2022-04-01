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


class AnnotationPfam(Base):
    def __init__(self, bind_object):
        super(AnnotationPfam, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + '_ref_rna'

    @report_check
    def add_annotation_pfam(self, name=None, params=None, stat_id=None):
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
            'name': name if name else 'AnnotationPfam_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'status': 'end',
            'desc': 'pfam注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_pfam']
        pfam_id = collection.insert_one(insert_data).inserted_id
        print "add sg_annotation_pfam!"
        return pfam_id

    def add_annotation_pfam_detail(self, pfam_id, pfam_path, seq_type, anno_type):
        """
        pfam_path: pfam_domain
        """
        if not isinstance(pfam_id, ObjectId):
            if isinstance(pfam_id, types.StringTypes):
                pfam_id = ObjectId(pfam_id)
            else:
                raise Exception('pfam_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(pfam_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(pfam_path))
        data_list = []
        with open(pfam_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = [
                    ('pfam_id', pfam_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('seq_id', line[0]),
                    ('pfam', line[2]),
                    ('domain', line[3]),
                    ('description', line[4])
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_pfam_detail']
            collection.insert_many(data_list)
        except Exception, e:
            self.logger.info.error("导入pfam注释信息:%s失败！" % pfam_path)
        else:
            self.logger.info.info("导入pfam注释信息:%s成功" % pfam_path)

    @report_check
    def add_annotation_pfam_bar(self, pfam_id, pfam_path, seq_type, anno_type):
        pfam = []
        domain = {}
        with open(pfam_path, "rb") as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split("\t")
                if line[3] not in pfam:
                    pfam.append(line[3])
                    domain[line[3]] = 1
                else:
                    domain[line[3]] += 1
        if not isinstance(pfam_id, ObjectId):
            if isinstance(pfam_id, types.StringTypes):
                pfam_id = ObjectId(pfam_id)
            else:
                raise Exception('pfam_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(pfam_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(pfam_path))
        data_list = []
        for i in domain:
            data = [
                ('pfam_id', pfam_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('domain', i),
                ('num', domain[i])
            ]
            data = SON(data)
            data_list.append(data)
        try:
            collection = self.db['sg_annotation_pfam_bar']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入pfam注释信息:%失败！" % pfam_path)
        else:
            self.bind_object.logger.info("导入pfam注释信息:%成功！" % pfam_path)

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


class AnnotationNr(Base):
    def __init__(self, bind_object):
        super(AnnotationNr, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + '_ref_rna'

    @report_check
    def add_annotation_nr(self, name=None, params=None, stat_id=None):
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
            'name': name if name else 'AnnotationNr_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'status': 'end',
            'desc': 'nr注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_nr']
        nr_id = collection.insert_one(insert_data).inserted_id
        print "add sg_annotation_nr!"
        return nr_id

    @report_check
    def add_annotation_nr_pie(self, nr_id, evalue_path, similar_path, seq_type, anno_type):
        """
        """
        if not isinstance(nr_id, ObjectId):
            if isinstance(nr_id, types.StringTypes):
                nr_id = ObjectId(nr_id)
            else:
                raise Exception('nr_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(evalue_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(evalue_path))
        if not os.path.exists(similar_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(similar_path))
        evalue, evalue_list, similar, similar_list = [], [], [], []
        data_list = []
        with open(evalue_path, "r") as f1, open(similar_path, "r") as f2:
            lines1 = f1.readlines()
            lines2 = f2.readlines()
            for line1 in lines1[1:]:
                line1 = line1.strip().split('\t')
                value = {"key": line1[0], "value": int(line1[1])}
                try:
                    value_list = {"key": line1[0], "value": line1[2]}
                except:
                    value_list = {"key": line1[0], "value": None}
                evalue.append(value)
                evalue_list.append(value_list)
            for line2 in lines2[1:]:
                line2 = line2.strip().split('\t')
                similarity = {"key": line2[0], "value": int(line2[1])}
                try:
                    similarity_list = {"key": line2[0], "value": line2[2]}
                except:
                    similarity_list = {"key": line2[0], "value": None}
                similar.append(similarity)
                similar_list.append(similarity_list)
        data = [
            ('nr_id', nr_id),
            ('seq_type', seq_type),
            ('anno_type', anno_type),
            ('e_value', evalue),
            ('similar', similar),
            ('evalue_list', evalue_list),
            ('similar_list', similar_list),
        ]
        data = SON(data)
        data_list.append(data)
        try:
            collection = self.db['sg_annotation_nr_pie']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入nr库注释作图信息evalue,similar：%s、%s出错!" % (evalue_path, similar_path))
        else:
            self.bind_object.logger.info("导入nr库注释作图信息evalue,similar：%s、%s成功!" % (evalue_path, similar_path))

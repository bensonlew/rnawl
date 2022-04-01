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


class AnnotationCog(Base):
    def __init__(self, bind_object):
        super(AnnotationCog, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + '_ref_rna'

    @report_check
    def add_annotation_cog(self, name=None, params=None, cog_sum=None):
        """
        cog注释导表函数
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        summary = None
        with open(cog_sum, "rb") as f:
            header = f.readline()
            total = re.match(r".*:(.*)$", header)
            if total:
                summary = total.group(1)
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationCog_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'status': 'end',
            'desc': 'cog注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'total_seq': summary
        }
        collection = self.db['sg_annotation_cog']
        cog_id = collection.insert_one(insert_data).inserted_id
        print "add ref_annotation_cog!"
        return cog_id

    @report_check
    def add_annotation_cog_detail(self, cog_id, cog_path, seq_type, anno_type):
        '''
        cog_path: cog_summary.xls
        seq_type: ref/new
        anno_type: transcript/gene
        '''
        if not isinstance(cog_id, ObjectId):
            if isinstance(cog_id, types.StringTypes):
                cog_id = ObjectId(cog_id)
            else:
                raise Exception('cog_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(cog_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(cog_path))
        data_list = list()
        with open(cog_path, 'r') as f:
            lines = f.readlines()
            for line in lines[2:]:
                line = line.strip().split('\t')
                data = [
                    ('cog_id', cog_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('type', line[0]),
                    ('function_categories', line[1]),
                    ('cog', int(line[2])),
                    ('nog', int(line[3])),
                    ('kog', int(line[4]))
                ]
                try:
                    data.append(('cog_list', line[5]))
                except:
                    data.append(('cog_list', None))
                try:
                    data.append(('nog_list', line[6]))
                except:
                    data.append(('nog_list', None))
                try:
                    data.append(('kog_list', line[7]))
                except:
                    data.append(('kog_list', None))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_cog_detail']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入cog注释信息：%s出错!" % (cog_path))
        else:
            self.bind_object.logger.info("导入cog注释信息：%s成功!" % (cog_path))

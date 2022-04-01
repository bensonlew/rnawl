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


class AnnotationKegg(Base):
    def __init__(self, bind_object):
        super(AnnotationKegg, self).__init__(bind_object)
        self._project_type = 'ref_rna'

        #self._db_name = Config().MONGODB + '_ref_rna'

    @report_check
    def add_annotation_kegg(self, name=None, params=None):
        """
        kegg注释导表函数
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationKegg_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'status': 'end',
            'desc': 'kegg注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        collection = self.db['sg_annotation_kegg']
        kegg_id = collection.insert_one(insert_data).inserted_id
        print "add ref_annotation_kegg!"
        return kegg_id

    @report_check
    def add_annotation_kegg_categories(self, kegg_id, seq_type, anno_type, categories_path):
        """
        categories_path:kegg_layer.xls
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(categories_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(categories_path))
        data_list = list()
        with open(categories_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                data = [
                    ('kegg_id', kegg_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('first_catergory', line[0]),
                    ('second_catergory', line[1]),
                    ('num', int(line[2])),
                    ('seq_list', line[3]),
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_kegg_categories']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入kegg注释分类信息：%s出错!" % (categories_path))
        else:
            self.bind_object.logger.info("导入kegg注释分类信息：%s 成功!" % categories_path)

    @report_check
    def add_annotation_kegg_level(self, kegg_id, seq_type, anno_type, level_path, png_dir):
        """
        level_path: pathway_table.xls
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(level_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(level_path))
        if not os.path.exists(png_dir):
            raise Exception('{}所指定的路径不存在，请检查！'.format(png_dir))
        data_list = []
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                fs = gridfs.GridFS(self.db)
                pid = re.sub('path:', '', line[0])
                pngid = fs.put(open(png_dir + '/' + pid + '.pdf', 'rb'))
                insert_data = {
                    'kegg_id': kegg_id,
                    'seq_type': seq_type,
                    'anno_type': anno_type,
                    'pathway_id': line[0],
                    'first_category': line[1],
                    'second_category': line[2],
                    'pathway_definition': line[3],
                    'number_of_seqs': int(line[4]),
                    'seq_list': line[5],
                    'graph_id': pngid,
                }
                data_list.append(insert_data)
        try:
            collection = self.db['sg_annotation_kegg_level']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入kegg注释层级信息：%s、%s出错!" % (level_path, png_dir))
        else:
            self.bind_object.logger.info("导入kegg注释层级信息：%s、%s 成功!" % (level_path, png_dir))

    @report_check
    def add_annotation_kegg_table(self, kegg_id, seq_type, anno_type, table_path):
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(table_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(table_path))
        with open(table_path, 'rb') as r:
            data_list = []
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                insert_data = {
                    'kegg_id': kegg_id,
                    'seq_type': seq_type,
                    'anno_type': anno_type,
                    'query_id': line[0],
                    'ko_id': line[1],
                    'ko_name': line[2],
                    'hyperlink': line[3],
                    'paths': line[4],
                }
                data_list.append(insert_data)
        try:
            collection = self.db['sg_annotation_kegg_table']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入kegg注释table信息：%s出错!" % (table_path))
        else:
            self.bind_object.logger.info("导入kegg注释table信息：%s成功!" % (table_path))

# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from collections import OrderedDict
import os
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import json
from api_base import ApiBase
import pandas as pd
import sqlite3


class Ngloc(ApiBase):
    def __init__(self, bind_object):
        super(Ngloc, self).__init__(bind_object)

    @report_check
    def add_annotation_subloc(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="52501046")
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationSubloc_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir + "/subloc",
            'status': 'start',
            'desc': 'subloc注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_subloc']
        subloc_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_annotation_subloc!")
        return subloc_id

    def get_description(self, acc_id):
        self.cursor.execute('select * from acc2des where uniacc="{}"'.format(acc_id))
        try:
            desc = self.cursor.fetchall()[0][2]
            description = desc
        except:
            description = ""
        return description

    def add_annotation_subloc_detail(self, subloc_id, subloc_path, seq_type=None, anno_type=None):
        """
        subloc_path: subloc_domain
        """
        self.bind_object.logger.info("start import")
        nracc2des = Config().SOFTWARE_DIR + '/database/Annotation/all/Uniprot/version_202009/uniprot_acc2des.db'
        self.bind_object.logger.info("start import {}".format(nracc2des))
        self.conn = sqlite3.connect(nracc2des)
        self.cursor = self.conn.cursor()
        if not isinstance(subloc_id, ObjectId):
            if isinstance(subloc_id, types.StringTypes):
                subloc_id = ObjectId(subloc_id)
            else:
                self.bind_object.set_error('subloc_id必须为ObjectId对象或其对应的字符串！', code="52501047")
        if not os.path.exists(subloc_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(subloc_path), code="52501048")
        data_list = []
        with open(subloc_path, "r") as f:
            lines = f.readlines()
            for line in lines[0:]:
                line = line.strip().split("\t")
                data = [
                    ('ngloc_id', subloc_id),
                    ('accession_id', line[0]),
                    ('description', self.get_description(line[0])),
                    ('subloc1', line[1].split(':')[0]),
                    ('subloc2', line[2].split(':')[0]),
                    ('subloc3', line[3].split(':')[0]),
                    ('subloc1_prob', float(line[1].split(':')[1])),
                    ('subloc2_prob', float(line[2].split(':')[1])),
                    ('subloc3_prob', float(line[3].split(':')[1])),
                ]
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['ngloc_detail']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入subloc注释信息:%s失败！" , variables=( subloc_path), code="52501049")
            else:
                self.bind_object.logger.info("导入subloc注释信息:%s成功" % subloc_path)

    @report_check
    def add_annotation_subloc_bar(self, subloc_id, subloc_path, seq_type=None, anno_type=None):
        subloc = []
        domain = {}
        if not isinstance(subloc_id, ObjectId):
            if isinstance(subloc_id, types.StringTypes):
                subloc_id = ObjectId(subloc_id)
            else:
                self.bind_object.set_error('subloc_id必须为ObjectId对象或其对应的字符串！', code="52501050")
        if not os.path.exists(subloc_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(subloc_path), code="52501051")
        data_list = []
        with open(subloc_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = [
                    ('ngloc_id', subloc_id),
                    ('name', line[0]),
                    ('value', int(line[1])),
                    ('type', "column")
                ]
                data = SON(data)
                data_list.append(data)

        if data_list:
            try:
                collection = self.db['ngloc_bar']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.set_error("导入subloc注释信息:%s失败！" , variables=( subloc_path), code="52501052")
            else:
                self.bind_object.logger.info("导入subloc注释信息:%s成功！" % subloc_path)



    def create_db_table(self, table_name, content_dict_list, tag_dict=None):
        '''
        Create main/detail table in database system.
        :param table_name: table name
        :param content_dict_list: list with dict as elements
        :param tag_dict: a dict to be added into each record in content_dict_list.
        :return: None or main table id
        '''
        table_id = None
        conn = self.db[table_name]
        if tag_dict:
            for row_dict in content_dict_list:
                row_dict.update(tag_dict)
        record_num = len(content_dict_list)
        try:
            if record_num > 5000:
                for i in range(0, record_num, 3000):
                    tmp_list = content_dict_list[i: i + 3000]
                    conn.insert_many(tmp_list)
            else:
                if record_num >= 2:
                    conn.insert_many(content_dict_list)
                else:
                    table_id = conn.insert_one(content_dict_list[0]).inserted_id
        except Exception as e:
            if record_num >= 2:
                self.bind_object.logger.warn('fail to insert records into table {} -> ({})'.format(table_name, e))
            else:
                self.bind_object.logger.warn('fail to insert record into table {} -> ({})'.format(table_name, e))
        else:
            if record_num >= 2:
                self.bind_object.logger.info('succeed in inserting records into table {}'.format(table_name))
            else:
                self.bind_object.logger.info('succeed in inserting record into table {}'.format(table_name))
            return table_id


if __name__ == '__main__':
    from mbio.workflows.denovo_rna_v2.denovo_test_api import DenovoTestApiWorkflow
    from biocluster.wsheet import Sheet
    import random


    data = {
        "id": "denovo_rna_v2_upgrade",
        #+ str(random.randint(1,10000)),
        #"id": "denovo_rna_v2",
        "project_sn": "denovo_rna_v2_upgrade",
        #+ str(random.randint(1,10000)),
        "type": "workflow",
        "name": "denovo_rna_v2.denovo_test_api",
        "options": {
        },
    }
    wsheet = Sheet(data=data)
    wf = DenovoTestApiWorkflow(wsheet)

    wf.IMPORT_REPORT_DATA = True
    wf.IMPORT_REPORT_AFTER_END = False
    wf.test_api = wf.api.api("denovo_rna_v2.cdslen")

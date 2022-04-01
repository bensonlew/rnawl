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


class OrfPredict(ApiBase):
    def __init__(self, bind_object):
        super(OrfPredict, self).__init__(bind_object)

    @report_check
    def add_orf_detail(self, cds_id, cds_detail):
        cds_pd = pd.read_table(cds_detail, header=0)
        cds_pd['cds_id'] = ObjectId(cds_id)
        result_dict_list = cds_pd.to_dict("records")
        self.create_db_table('orf_predict_detail', result_dict_list)

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

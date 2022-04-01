# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.dia.api_base import ApiBase
import re
import pandas as pd
from bson.objectid import ObjectId

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class Express(ApiBase):
    def __init__(self, bind_object):
        super(Express, self).__init__(bind_object)
        self._project_type = 'dia'

    #@report_check
    def add_express(self, params=None, project_sn=None,
                    task_id=None, express_exp=None,
                    main_id=None,name=None):
        if params is None:
            params = {"software": "peaks 8.5"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'params': params if params else "",
            'status': 'start',
            'desc': 'exp main table',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            'type': 'raw',
            'version':'v3',
            'name':name,
            }
        express_id = self.create_db_table('sg_express', [insert_data])
        self.add_express_detail(express_exp, express_id)
        df = pd.read_table(express_exp, sep ="\t", header = 0)
        self.update_db_record('sg_express', express_id,
                              status="end", main_id=express_id, samples=list(df.columns[1:]))
        return express_id

    #@report_check
    def add_express_detail(self, express_exp, express_id):
        df = pd.read_table(express_exp, sep="\t", header=0, dtype={0: str})
        df.rename(columns={df.columns[0]: "accession_id"}, inplace=True)
        df = df.round(6)
        df['express_id'] = ObjectId(express_id)
        data_list = df.to_dict('records')
        self.create_db_table('sg_express_detail', data_list)

if __name__ == '__main__':
    anno = Express(None)
    express_exp = "/mnt/ilustre/users/sanger-dev/workspace/20201130/Diav3_202011301316/raw_treat_ref"
    project_sn="test_api_database___project_sn"
    task_id="test_api_database___task_id"
    anno.add_express(params=None, main_id=None,
                                 project_sn=project_sn, task_id=task_id,
                                 express_exp=express_exp,name="Raw")

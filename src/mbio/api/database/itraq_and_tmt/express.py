# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.itraq_and_tmt.api_base import ApiBase
import re
import pandas as pd
from bson.objectid import ObjectId

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class Express(ApiBase):
    def __init__(self, bind_object):
        super(Express, self).__init__(bind_object)
        self._project_type = 'itraq_tmt'

    #@report_check
    def add_express(self, params=None, project_sn=None,
                    task_id=None, type=None, express_exp=None,
                    main_id=None):
        if params is None:
            params = {"software": "Proteome Discoverer"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'params': params if params else "",
            'status': 'start',
            'desc': 'express结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            'type': type,
            }
        express_id = self.create_db_table('sg_express', [insert_data])
        self.add_express_detail(express_exp, express_id)
        df = pd.read_table(express_exp, sep ="\t", header = 0)
        self.update_db_record('sg_express', express_id,
                              status="end", main_id=express_id, samples=list(df.columns[1:]))

    #@report_check
    def add_express_detail(self, express_exp, express_id):
        df = pd.read_table(express_exp, sep ="\t", header = 0, dtype = {0:str})
        df.rename(columns={df.columns[0]: "accession_id"}, inplace=True)
        df = df.round(6)
        df['express_id'] = ObjectId(express_id)
        data_list = df.to_dict('records')
        self.create_db_table('sg_express_detail', data_list)

if __name__ == '__main__':
    anno = Express(None)
    express_exp = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_itraq_and_tmt/data4/result/ratio.exp.txt'
    type="ratio"
    task_id="tsg_28978"
    project_sn="188_5ab1e54fa8afd"
    anno.add_express(express_exp=express_exp, type=type,project_sn=project_sn, task_id=task_id)

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
class Searchdb(ApiBase):
    def __init__(self, bind_object):
        super(Searchdb, self).__init__(bind_object)
        self._project_type = 'dia'

    #@report_check
    def add_searchdb(self, params=None, project_sn='dia',
                task_id='dia', name=None, protein_selected = None):
        if params is None:
            params = {"software": "Spectronaut"}
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'params': params,
            'status': 'start',
            'desc': 'search_db结果表',
            'name': name if name else 'search_db_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            'version':'v3'
            }
        searchdb_id = self.create_db_table('sg_search_db', [insert_data])
        print("9999999999999")
        print(searchdb_id)
        print("666666666666")
        self.add_searchdb_detail(protein_selected, searchdb_id)
        self.update_db_record('sg_search_db', searchdb_id,
                              status="end", main_id=searchdb_id)
        # 更新的表虽然也更新了时间，但是sort时候还是会排在后面
        '''
        self.update_sgtask_record(task_id=task_id,
                                  protein_sliced=protein_selected,
                                  created_ts=datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'))
        # self.update_db_record('sg_task', task_id, protein_selected=protein_selected, main_id=searchdb_id)
        print("更新sg_task表成功")
        '''
    #@report_check
    def add_searchdb_detail(self, protein_selected, searchdb_id):
        df = pd.read_table(protein_selected, sep="\t", header=0, dtype={0:str}).fillna("_")
        df = df.round(6)
        df.columns = ["accession_id", "description", "peptide", "aas", "mw", "calc_pi"]
        df['searchdb_id'] = ObjectId(searchdb_id)
        data_list = df.to_dict('records')
        self.create_db_table('sg_search_db_detail', data_list)

if __name__ == '__main__':
    anno = Searchdb(None)
    task_id="sanger_146293"
    project_sn="25727_5c12060282185"
    protein_selected = '/mnt/lustre/users/sanger/workspace/20181214/Itraqtmt_sanger_146293/Searchdb/protein_sliced.xls'
    anno.add_searchdb(protein_selected=protein_selected, task_id=task_id, project_sn=project_sn)


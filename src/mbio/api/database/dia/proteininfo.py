# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.dia.api_base import ApiBase
import pandas as pd
import json
import shutil
import os

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class Proteininfo(ApiBase):
    def __init__(self, bind_object):
        super(Proteininfo, self).__init__(bind_object)
        self._project_type = 'dia'

    #@report_check
    def add_proteininfo(self, params=None, project_sn='dia',
                    task_id='dia', main_id=None, proteininfo_exp=None):

        df = pd.read_table(proteininfo_exp, sep ="\t", header = 0)
        # 感觉有一次修改之后，列名就对不上了，没办法就只能暴力解决了
        old_columns = ['Total Spectrum', 'Identified Spectrum', 'Peptide number',
                       'Protein number', 'Protein group number',
                       'Peptide-Spectrum Matches', 'Peptide sequences', 'Proteins', 'Protein groups']
        new_columns = ['total_spectrum', 'identified_spectrum','peptide_num', 'protein_num', 'protein_group_num',
                       'identified_spectrum', 'peptide_num', 'protein_num', 'protein_group_num']
        columns=zip(old_columns, new_columns)
        df_col = df.columns.tolist()
        # df.rename(columns=dict(columns),inplace=True)
        for o, n in columns:
            if o in df_col:
                df[n] = df[o]
        # 咱也不知道前端那边咋弄的啊，还是多加几列吧
        try:
            df['Peptide sequence'] = df['peptide_num']
            df['Peptide sequences'] = df['peptide_num']
            df['Protein groups'] = df['protein_group_num']
            df['Proteins'] = df['protein_num']
        except:
            pass
        data_list = df.to_dict('records')
        info_id = self.create_db_table('sg_protein_info', data_list)
        if params is None:
            params = {"software": "peaks 8.5"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'params': params if params else "",
            'status': 'start',
            'desc': 'proteininfo结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            'version':'v3'
            }
        self.update_db_record_2('sg_protein_info', info_id, insert_data,
                                main_id=info_id)
        self.update_db_record('sg_protein_info', info_id, status='end')

        qc_dir = os.path.join(self.bind_object.work_dir, 'qc')
        if os.path.exists(qc_dir):
            shutil.copy(proteininfo_exp, os.path.join(qc_dir, 'Protein_infomation.xls'))

        return data_list


if __name__ == '__main__':
    anno = Proteininfo(None)
    proteininfo_exp = '/mnt/ilustre/users/sanger-dev/sg-users/litangjian' \
                '/protein_dev/Protein_information.xls'
    anno.add_proteininfo(proteininfo_exp=proteininfo_exp)


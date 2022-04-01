# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'

import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.ref_rna_v2.api_base import ApiBase
import os
import pandas as pd
from bson.objectid import ObjectId
from collections import OrderedDict
from biocluster.config import Config
from bson.objectid import ObjectId

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class GenesetSelf(ApiBase):
    def __init__(self, bind_object):
        super(GenesetSelf, self).__init__(bind_object)
        self._project_type = 'ref_rna_v2'

    #@report_check
    def add_geneset(self, geneset_output_dir, name=None, type='transcript_id',
                    project_sn='ref_rna_v2', task_id='ref_rna_v2', main_id=None):
        if main_id is None:
            # prepare main table info
            time=datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            desc="the_geneset_upload_at_%s_by_the_customer"%time
            if name is None:
                name=desc

            main_info = dict(
                name=name,
                task_id=task_id,
                project_sn=project_sn,
                group_id="",
                type=type,
                desc=desc,
            )
            main_id = self.create_db_table('sg_geneset', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        # prepare detail table info
        # target_file = os.path.join(geneset_output_dir, 'GenesetSelf.txt')
        with open(geneset_output_dir, 'r') as f:
            f.readline()
            seq_list = [x.strip() for x in f]
            gene_length=len(seq_list)
        if gene_length < 5:
            self.db['sg_geneset'].update({'_id': main_id}, {'$set': {'params': None, 'status': 'failed'}})
            self.bind_object.set_error("预上传基因集中基因与该项目中基因匹配数目过少（＜5个），不予上传，请核查ID书写是否规范或者从交互页面直接创建")
        row_dict_list = [{"seq_list": seq_list}]
        # insert detail
        tag_dict = dict(geneset_id=main_id)
        self.create_db_table('sg_geneset_detail', row_dict_list, tag_dict=tag_dict)
        self.update_db_record('sg_geneset', main_id, status="end", main_id=main_id, gene_length=gene_length)
        # status_record_id = self.db['sg_status'].find_one({'table_id': main_id})['_id']
        # self.update_db_record('sg_status', status_record_id, status="end")

if __name__ == '__main__':
    pass
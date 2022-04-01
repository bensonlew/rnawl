# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify:20161220
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
from biocluster.config import Config
from bson.objectid import ObjectId
import re


class DenovoGoPvalSort(Base):
    def __init__(self, bind_object):
        super(DenovoGoPvalSort, self).__init__(bind_object)
        self._db_name = Config().MONGODB + '_rna'

    @report_check
    def add_pvalue(self, name=None, params=None, pval_sort=None, express_diff_id=None):
        '''
        pval_sort：go显著富集分析结果表
        express_diff_id: 用于工作流初始化代码截停时，更新params
        '''
        if express_diff_id:
            params['express_diff_id'] = str(express_diff_id)
        project_sn = self.bind_object.sheet.project_sn
        task_id = self.bind_object.sheet.id
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'GoPvalueSort_' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'status': 'end',
            'desc': 'go富集调控pval排序',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db['sg_denovo_go_pvalue']
        sort_id = collection.insert_one(insert_data).inserted_id
        if os.path.exists(pval_sort):
            self.add_pvalue_detail(sort_id, pval_sort)
        return sort_id
        self.bind_object.logger.info("add sg_denovo_go_pvalue sucess!")

    @report_check
    def add_pvalue_detail(self, sort_id, pval_sort):
        if not isinstance(sort_id, ObjectId):
            if isinstance(sort_id, types.StringTypes):
                sort_id = ObjectId(sort_id)
            else:
                raise Exception('sort_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(pval_sort):
            raise Exception('pval_sort所指定的路径:{}不存在，请检查！'.format(pval_sort))
        data_list = []
        with open(pval_sort, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                insert_data = {
                    'go_pvalue_id': sort_id,
                    'term_type': line[0],
                    'term': line[1],
                    'go_id': line[2],
                    'pvalue': round(float(line[8]), 6),
                    'corrected_pvalue': round(float(line[7]), 6),
                    'up_num': int(line[3]),
                    'down_num': int(line[5]),
                    'up_percent': round(float(line[4]), 6),
                    'down_percent': round(float(line[6]), 6)
                }
                data_list.append(insert_data)
        try:
            collection = self.db["sg_denovo_go_pvalue_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入go富集调控排序表：%s出错:%s" % (pval_sort, e))
        else:
            self.bind_object.logger.info("导入go富集调控排序表:%s成功" % pval_sort)

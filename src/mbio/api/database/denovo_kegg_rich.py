# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# last_modify:20161124
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
from biocluster.config import Config
from bson.objectid import ObjectId
import re
import json


class DenovoKeggRich(Base):
    def __init__(self, bind_object):
        super(DenovoKeggRich, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + '_rna'

    @report_check
    def add_kegg_rich(self, name=None, params=None, kegg_enrich_table=None, express_diff_id=None):
        '''
        kegg_enrich_table：kegg富集的分析结果表
        express_diff_id: 用于工作流初始化代码截停时，更新params
        '''
        if express_diff_id:
            params['express_diff_id'] = str(express_diff_id)
        project_sn = self.bind_object.sheet.project_sn
        task_id = self.bind_object.sheet.id
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'KeggEnrich_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'status': 'end',
            'desc': 'kegg富集分析',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db['sg_denovo_kegg_enrich']
        enrich_id = collection.insert_one(insert_data).inserted_id
        if kegg_enrich_table:
            self.add_kegg_enrich_detail(enrich_id, kegg_enrich_table)
        self.bind_object.logger.info("add sg_denovo_kegg_enrich sucess!")
        return enrich_id

    def add_kegg_enrich_detail(self, enrich_id, kegg_enrich_table):
        if not isinstance(enrich_id, ObjectId):
            if isinstance(enrich_id, types.StringTypes):
                enrich_id = ObjectId(enrich_id)
            else:
                raise Exception('kegg_enrich_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(kegg_enrich_table):
            raise Exception('kegg_enrich_table所指定的路径:{}不存在，请检查！'.format(kegg_enrich_table))
        data_list = []
        with open(kegg_enrich_table, 'rb') as r:
            for line in r:
                if re.match(r'\w', line):
                    line = line.strip('\n').split('\t')
                    insert_data = {
                        'kegg_enrich_id': enrich_id,
                        'term': line[0],
                        'database': line[1],
                        'id': line[2],
                        'study_number': int(line[3]),
                        'backgroud_number': int(line[4]),
                        'pvalue': round(float(line[5]), 4),
                        'corrected_pvalue': round(float(line[6]), 4),
                        'gene_lists': line[7],
                        'hyperlink': line[8]
                    }
                    data_list.append(insert_data)
            if data_list:
                try:
                    collection = self.db['sg_denovo_kegg_enrich_detail']
                    collection.insert_many(data_list)
                except Exception, e:
                    self.bind_object.logger.error("导入kegg富集统计表：%s信息出错:%s" % (kegg_enrich_table, e))
                else:
                    self.bind_object.logger.info("导入kegg富集统计表:%s信息成功!" % kegg_enrich_table)
            else:
                coll = self.db['sg_denovo_kegg_enrich']
                coll.update({'_id': enrich_id}, {'$set': {'desc': 'no_result'}})
                self.bind_object.logger.info("kegg富集统计表没结果：%s" % kegg_enrich_table)

# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# last_modify:20161124
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
from biocluster.config import Config
from bson.objectid import ObjectId
import gridfs
import json


class DenovoKeggRegulate(Base):
    def __init__(self, bind_object):
        super(DenovoKeggRegulate, self).__init__(bind_object)
        self._db_name = Config().MONGODB + '_rna'

    @report_check
    def add_kegg_regulate(self, name=None, params=None, kegg_regulate_table=None, pathways_dir=None, express_diff_id=None):
        '''
        kegg_regulate_table：kegg调控的分析结果表
        express_diff_id: 用于工作流初始化代码截停时，更新params
        '''
        if express_diff_id:
            params['express_diff_id'] = str(express_diff_id)
        project_sn = self.bind_object.sheet.project_sn
        task_id = self.bind_object.sheet.id
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'KeggRegulate_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'status': 'end',
            'desc': 'kegg调控统计分析',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db['sg_denovo_kegg_regulate']
        regulate_id = collection.insert_one(insert_data).inserted_id
        if kegg_regulate_table:
            self.add_kegg_regulate_detail(regulate_id, kegg_regulate_table)
        if pathways_dir:
            self.add_kegg_regulate_pathway(regulate_id, pathways_dir)
        self.bind_object.logger.info("add sg_denovo_kegg_regulate sucess!")
        return regulate_id

    def add_kegg_regulate_detail(self, regulate_id, kegg_regulate_table):
        if not isinstance(regulate_id, ObjectId):
            if isinstance(regulate_id, types.StringTypes):
                regulate_id = ObjectId(regulate_id)
            else:
                raise Exception('kegg_regulate_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(kegg_regulate_table):
            raise Exception('kegg_regulate_table所指定的路径:{}不存在，请检查！'.format(kegg_regulate_table))
        data_list = []
        with open(kegg_regulate_table, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                insert_data = {
                    'kegg_regulate_id': regulate_id,
                    'pathway_id': line[0],
                    'ko_ids': line[1],
                    'up_numbers': int(line[2]),
                    'down_numbers': int(line[3]),
                    'up_genes': line[4],
                    'down_genes': line[5],
                }
                data_list.append(insert_data)
            try:
                collection = self.db['sg_denovo_kegg_regulate_detail']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.logger.error("导入kegg调控统计表：%s信息出错:%s" % (kegg_regulate_table, e))
            else:
                self.bind_object.logger.info("导入kegg调控统计表:%s信息成功!" % kegg_regulate_table)

    def add_kegg_regulate_pathway(self, regulate_id, pathway_dir):
        if not isinstance(regulate_id, ObjectId):
            if isinstance(regulate_id, types.StringTypes):
                regulate_id = ObjectId(regulate_id)
            else:
                raise Exception('kegg_regulate_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(pathway_dir):
            raise Exception('pathway_dir所指定的路径:{}不存在，请检查！'.format(pathway_dir))
        data_list = []
        files = os.listdir(pathway_dir)
        fs = gridfs.GridFS(self.db)
        for f in files:
            png_id = fs.put(open(os.path.join(pathway_dir, f), 'rb'))
            insert_data = {
                'kegg_regulate_id': regulate_id,
                'pathway_png': png_id,
                'pathway_id': f.split('.pdf')[0]
            }
            data_list.append(insert_data)
        try:
            collection = self.db['sg_denovo_kegg_regulate_pathway']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入kegg调控pathway：%s信息出错:%s" % (pathway_dir, e))
        else:
            self.bind_object.logger.info("导入kegg调控pathway:%s信息成功!" % pathway_dir)

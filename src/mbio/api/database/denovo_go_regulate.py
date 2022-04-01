# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify:20161205
import os
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import json


class DenovoGoRegulate(Base):
    def __init__(self, bind_object):
        super(DenovoGoRegulate, self).__init__(bind_object)
        #self._db_name = Config().MONGODB + '_rna'
        self._project_type = 'ref_rna'

    @report_check
    def add_go_regulate(self, name=None, params=None, go_regulate_dir=None, express_diff_id=None):
        '''
        go_regulate_dir：go分类统计的分析结果表
        express_diff_id: 用于工作流初始化代码截停时，更新params
        '''
        if express_diff_id:
            params['express_diff_id'] = str(express_diff_id)
        project_sn = self.bind_object.sheet.project_sn
        task_id = self.bind_object.sheet.id
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'GoRegulate_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'status': 'end',
            'desc': 'go调控分析主表',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db['sg_denovo_go_regulate']
        go_regulate_id = collection.insert_one(insert_data).inserted_id
        if os.path.exists(go_regulate_dir):
            self.add_go_regulate_detail(go_regulate_id, go_regulate_dir)
        return go_regulate_id

    @report_check
    def add_go_regulate_detail(self, go_regulate_id, go_regulate_dir):
        if not isinstance(go_regulate_id, ObjectId):
            if isinstance(go_regulate_id, types.StringTypes):
                go_regulate_id = ObjectId(go_regulate_id)
            else:
                raise Exception('go_enrich_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(go_regulate_dir):
            raise Exception('{}所指定的路径不存在，请检查！'.format(go_regulate_dir))
        data_list = []
        with open(go_regulate_dir, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                line[3] = int(line[3])
                line[4] = round(float(line[4]), 6)
                line[5] = int(line[5])
                line[6] = round(float(line[6]), 6)
                data = [
                    ('go_regulate_id', go_regulate_id),
                    ('go_type', line[0]),
                    ('go', line[1]),
                    ('go_id', line[2]),
                    ('up_num', line[3]),
                    ('up_percent', line[4]),
                    ('down_num', line[5]),
                    ('down_percent', line[6]),
                ]
                try:
                    data += [('up_genes', line[7])]
                except:
                    data += [('up_genes', '')]
                try:
                    data += [('down_genes', line[8])]
                except:
                    data += [('down_genes', '')]
                data = SON(data)
                data_list.append(data)
            try:
                collection = self.db['sg_denovo_go_regulate_detail']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入go调控信息：%s出错:%s" % (go_regulate_dir, e))
            else:
                self.bind_object.logger.info("导入go调控信息：%s成功!" % (go_regulate_dir))

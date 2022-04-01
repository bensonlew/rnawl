# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify:20161205
import os
import re
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
import gridfs
import json
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config


class DenovoGoEnrich(Base):
    def __init__(self, bind_object):
        super(DenovoGoEnrich, self).__init__(bind_object)
        self._project_type = 'ref_rna'

    @report_check
    def add_go_enrich(self, name=None, params=None, go_graph_dir=None, go_enrich_dir=None, express_diff_id=None):
        '''
        go_graph_dir: go有向无环图的路径
        go_enrich_dir：go富集的分析结果表
        express_diff_id: 用于工作流初始化代码截停时，更新params
        '''
        if express_diff_id:
            params['express_diff_id'] = str(express_diff_id)
        project_sn = self.bind_object.sheet.project_sn
        task_id = self.bind_object.sheet.id
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'GOEnrich_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'status': 'end',
            'desc': 'go富集分析主表',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db['sg_denovo_go_enrich']
        go_enrich_id = collection.insert_one(insert_data).inserted_id
        if go_graph_dir:
            fs = gridfs.GridFS(self.db)
            gra = fs.put(open(go_graph_dir, 'rb'))
            try:
                collection.update({"_id": ObjectId(go_enrich_id)}, {"$set": {'go_directed_graph': gra}})
            except Exception, e:
                self.bind_object.set_error("导入%s信息出错：%s" % (go_graph_dir, e))
            else:
                self.bind_object.logger.info("导入%s信息成功！" % (go_graph_dir))
        if os.path.exists(go_enrich_dir):
            self.add_go_enrich_detail(go_enrich_id, go_enrich_dir)
        self.bind_object.logger.info("add sg_denovo_go_enrich sucess!")
        return go_enrich_id

    @report_check
    def update_directed_graph(self, go_enrich_id, go_graph_dir):
        collection = self.db['sg_denovo_go_enrich']
        if go_graph_dir:
            fs = gridfs.GridFS(self.db)
            gra = fs.put(open(go_graph_dir, 'rb'))
            try:
                collection.update({"_id": ObjectId(go_enrich_id)}, {"$set": {'go_directed_graph': gra}})
            except Exception, e:
                self.bind_object.set_error("导入%s信息出错：%s" % (go_graph_dir, e))
            else:
                self.bind_object.logger.info("导入%s信息成功！" % (go_graph_dir))

    @report_check
    def add_go_enrich_detail(self, go_enrich_id, go_enrich_dir):
        if not isinstance(go_enrich_id, ObjectId):
            if isinstance(go_enrich_id, types.StringTypes):
                go_enrich_id = ObjectId(go_enrich_id)
            else:
                raise Exception('go_enrich_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(go_enrich_dir):
            raise Exception('{}所指定的路径不存在。请检查！'.format(go_enrich_dir))
        data_list = []
        with open(go_enrich_dir, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                if float(line[8]):
                    m = re.match(r"(.+)/(.+)", line[5])
                    pop_count = int(m.group(1))
                    line[6] = round(float(line[6]), 6)
                    line[7] = int(line[7])
                    line[8] = int(line[8])
                    data = [
                        ('go_enrich_id', go_enrich_id),
                        ('go_id', line[0]),
                        ('go_type', line[1]),
                        ('enrichment', line[2]),
                        ('discription', line[3]),
                        ('ratio_in_study', line[4]),
                        ('ratio_in_pop', line[5]),
                        ('p_uncorrected', line[6]),
                        ('depth', line[7]),
                        ('study_count', line[8]),
                        ('pop_count', pop_count),
                        ('diff_genes', line[-1]),
                    ]
                    try:
                        data += [('p_fdr', round(float(line[-2]), 6))]
                    except:
                        data += [('p_fdr', None)]
                    data = SON(data)
                    data_list.append(data)
            try:
                collection = self.db['sg_denovo_go_enrich_detail']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入go富集信息：%s出错:%s" % (go_enrich_dir, e))
            else:
                self.bind_object.logger.info("导入go富集信息：%s成功!" % (go_enrich_dir))

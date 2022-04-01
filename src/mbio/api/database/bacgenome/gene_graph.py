# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20180313
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId


class GeneGraph(Base):
    def __init__(self, bind_object):
        super(GeneGraph, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_gene_graph(self, fasta_path,kegg_path,cog_path, main=False, task_id=None, main_id=None,project_sn=None,
                    params=None, name=None,specimen_id =None):
        collection = self.db['gene_graph']
        if main:
            if task_id is None:
                task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn if project_sn is None else project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '基因组图谱',
                'created_ts': created_ts,
                'name': name if name else 'GeneGraph_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end'
            }
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
            collection.update({'_id': main_id},{'$set':{'main_id':main_id}}) #guanqing 20180813
        else:
            if main_id is None:
                #raise Exception("main为False时需提供main_id!")
                self.bind_object.set_error("main为False时需提供main_id!", code="51402601")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        collection_detail = self.db["gene_graph_detail"]
        data_list = []
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                data = [('gene_graph_id',main_id),('specimen_id',specimen_id)]
                if line.startswith('>'):
                    line = line.strip().split(' ')
                    location = line[3].split('_')[0]
                    data.extend(
                        [("location", location), ("gene_id", line[0][1:]), ("from", int(line[4])), ("to", int(line[5]))])
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入序列信息出错！")
            self.bind_object.set_error("导入序列信息出错！", code="51402602")
        else:
            self.bind_object.logger.info("导入序列信息成功!")
        with open(kegg_path, 'r') as f1:
            lines = f1.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                if line[1] != '-':
                    collection_detail.update({'gene_graph_id':main_id,'gene_id':line[0],'specimen_id':specimen_id},{'$set':{'gene_name':line[1]}})
        with open(cog_path, 'r') as f2:
            lines = f2.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                collection_detail.update({'gene_graph_id':main_id,'gene_id':line[0],'specimen_id':specimen_id},{'$set':{'cog_type':line[3][0]}})
        return main_id

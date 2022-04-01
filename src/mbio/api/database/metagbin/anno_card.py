# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId


class AnnoCard(Base):
    def __init__(self, bind_object):
        super(AnnoCard, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_anno_card(self, anno_card_path, main=False, task_id=None, project_sn=None, main_id=None,
                      params=None, name=None, update_id=None, genome_id=None, database_version="new"):
        if main:
            if task_id is None:
                task_id = self.bind_object.sheet.id
            project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': 'CARD注释',
                'created_ts': created_ts,
                'name': name if name else 'AnnoCard_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end'
            }
            collection = self.db['anno_card']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
            collection.update({'_id': main_id},{'$set':{'main_id':main_id}})
        else:
            if main_id is None:
                raise Exception("main为False时需提供main_id!")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        self.genome_id = genome_id
        if os.path.exists(anno_card_path + '/' + genome_id + '_card_anno.xls'):
            self.add_anno_card_type(anno_card_path + '/' + genome_id + '_card_category.xls', main_id, genome_id=genome_id)
            self.add_anno_card_detail(anno_card_path + '/' + genome_id + '_card_anno.xls', main_id, update_id=update_id,
                                  genome_id=genome_id, database_version=database_version)
        return main_id

    def add_anno_card_type(self, type_file, main_id=None, genome_id=None):
        data_list = []
        pie_list = []
        collection_type = self.db['anno_card_type']
        collection_pie = self.db['anno_card_pie']
        type_num = {}
        pie_num = {}
        with open(type_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                data = [("anno_card_id", main_id)]
                line = line.strip().split('\t')
                type_num[line[1]] = int(line[2])
                data.extend(
                    [("genome_id", genome_id), ("aro", line[1]),
                     ("anno_num", int(line[2]))])
                data_son = SON(data)
                data_list.append(data_son)
            sort_num = sorted(type_num.items(), key=lambda e: e[1], reverse=True)
            if len(type_num) > 11:
                n = 1
                for i in sort_num:
                    if n < 11:
                        pie_num[i[0]] = type_num[i[0]]
                    else:
                        pie_num['other'] = pie_num['other'] + type_num[i[0]] if pie_num.has_key('other') else type_num[
                            i[0]]
                    n = n + 1
            else:
                pie_num = type_num
            pie_data = [("anno_card_id", main_id)]
            for key in pie_num:
                pie_data.extend(
                    [("genome_id", genome_id), ("aro", key),
                     ("anno_num", pie_num[key])])
                data_son = SON(pie_data)
                pie_list.append(data_son)
        try:
            collection_type.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入anno_card_type%s信息出错:%s" % (type_file, e))
        else:
            self.bind_object.logger.info("导入anno_card_type%s信息成功!" % type_file)
        try:
            collection_pie.insert_many(pie_list)
        except Exception as e:
            self.bind_object.logger.error("导入anno_card_pie%s信息出错:%s" % (type_file, e))
        else:
            self.bind_object.logger.info("导入anno_card_pie%s信息成功!" % type_file)

    def add_anno_card_detail(self, detail_file, main_id=None, update_id=None, genome_id=None, database_version="new"):
        data_list = []
        collection_detail = self.db['anno_card_detail']
        summary_collection = self.db['anno_summary_detail']
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                data = [("anno_card_id", main_id)]
                line = line.strip().split('\t')
                data.extend(
                    [("gene_id", line[0]), ("location", line[1]), ("genome_id", genome_id),
                     ("aro_name", line[3]),
                     ("accession", line[4]), ("des", line[5])])
                if database_version  in ['new']:
                    data.extend([("category", line[8]), ("identity", float(line[7])),("evalue", float(line[6])),
                                 ("resis_mechanism", line[9])])
                else:
                    data.extend([("category", line[6]), ("identity", float(line[8])),("evalue", float(line[7]))])
                data_son = SON(data)
                data_list.append(data_son)
                if update_id:
                    if not isinstance(update_id, ObjectId):
                        update_id = ObjectId(update_id)
                    summary_collection.update_one(
                        {'summary_id': update_id, 'gene_id': line[0], 'genome_id': genome_id},
                        {'$set': {'aro_name': line[3], 'aro_category': line[6]}},
                        upsert=False)
        try:
            collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入anno_card_detail%s信息出错:%s" % (detail_file, e))
        else:
            self.bind_object.logger.info("导入anno_card_detail%s信息成功!" % detail_file)

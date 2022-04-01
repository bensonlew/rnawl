# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20180310
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
        self._project_type = "bacgenome"

    @report_check
    def add_anno_card(self, anno_card_path, main=False, task_id=None, project_sn=None, main_id=None,
                      params=None, name=None, update_id=None, specimen_id=None):
        self.specimen_id = specimen_id
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
                'status': 'end',
                "version": "3.0",
                "settled_params": json.dumps({"version": "card_v3.0.9"})
            }
            collection = self.db['anno_card']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
            collection.update({'_id': main_id},{'$set':{'main_id':main_id}}) #guanqing 20180813
        else:
            if main_id is None:
                raise Exception("main为False时需提供main_id!")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        self.add_anno_card_type(anno_card_path + '/' + specimen_id + '_card_category.xls', main_id)
        self.add_anno_card_detail(anno_card_path + '/' + specimen_id + '_card_anno.xls', main_id, update_id=update_id)
        self.add_anno_card_stat(anno_card_path + '/sample_stat.xls',main_id)
        return main_id

    def add_anno_card_type(self, type_file, main_id=None):
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
                    [("specimen_id", self.specimen_id), ("drug_class", line[1]),
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
                    [("specimen_id", self.specimen_id), ("drug_class", key),
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

    def add_anno_card_detail(self, detail_file, main_id=None, update_id=None):
        data_list = []
        collection_detail = self.db['anno_card_detail']
        summary_collection = self.db['anno_summary_detail']
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                data = [("anno_card_id", main_id)]
                line = line.strip().split('\t')
                data.extend(
                    [("gene_id", line[0]), ("location", line[1]), ("specimen_id", self.specimen_id),
                     ("aro_name", line[3]),
                     ("accession", line[4]), ("des", line[5]), ("identity", float(line[7])),
                     ("evalue", float(line[6]))])   #guanqing.zou 20180911 互换 identity 和 evalue的值
                      ## qingchen.zhang 去掉category
                if len(line) > 10:
                    data.extend([("drug_class", line[8]),("resis_mechanism", line[9]),("score", float(line[10])),("coverage", float(line[11]))])
                data_son = SON(data)
                data_list.append(data_son)
                if update_id:
                    pass
                    # if not isinstance(update_id, ObjectId):
                    #     update_id = ObjectId(update_id)
                    # summary_collection.update_one(
                    #     {'summary_id': update_id, 'gene_id': line[0], 'specimen_id': self.specimen_id},
                    #     {'$set': {'aro_name': line[3], 'aro_category': line[6]}},
                    #     upsert=False)
        try:
            collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入anno_card_detail%s信息出错:%s" % (detail_file, e))
        else:
            self.bind_object.logger.info("导入anno_card_detail%s信息成功!" % detail_file)


    def add_anno_card_stat(self,stat_file,main_id):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        with open(stat_file) as f:
            f.readline()
            line2 = f.readline()
            sp_line = line2.split('\t')
            sample = sp_line[0]
            gene_num = sp_line[1]

        insert_data = {
            "specimen_id" : sample,
            "gene_num" : int(gene_num),
            "anno_card_id" : main_id
        }
        collection_stat = self.db['anno_card_stat']
        try:
            collection_stat.insert_one(insert_data)
        except Exception as e:
            self.bind_object.logger.error('导anno_card_stat表失败')
        else:
            self.bind_object.logger.info("导入anno_card_stat表成功")


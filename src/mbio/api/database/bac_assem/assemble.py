# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
# last_modify:20190417
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId
import pandas as pd


class Assemble(Base):
    def __init__(self, bind_object):
        super(Assemble, self).__init__(bind_object)
        self._project_type = "bac_assem"

    @report_check
    def add_draft(self, params, task_id=None, project_sn=None, name=None):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': 'draft_assem',
            'created_ts': created_ts,
            'name': name if name else 'draft_assem',
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end'
        }
        collection = self.db['draft']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': main_id},{'$set':{'main_id':main_id}})
        return main_id

    def check_id(self, main_id):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！')
        return main_id

    @report_check
    def add_draft_stat_detail(self, main_id, type, stat_file, samp, seq_path):
        main_id = self.check_id(main_id)
        data = pd.read_table(stat_file, header=0, names=["scaf_no", "scaf_base", "scaf_no_large", "scaf_base_large",
                                                         "scaf_len_max", "scaf_n50", "scaf_n90", "gc_rate", "n_rate",
                                                         "contig_no", "contig_base", "contig_no_large",
                                                         "contig_base_large", "contig_len_max", "contig_n50", "contig_n90"])
        data["draft_id"] = main_id
        data["sof_type"] =type
        data["samp"] = samp
        # data["coverage"] = 1  # 目前表格里没有这个数据
        data["seq_path"] = seq_path
        data_list = data.to_dict(orient="records")
        try:
            self.db["draft_stat_detail"].insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("draft_stat_detail export %s error: %s" % (stat_file, e))
        else:
            self.bind_object.logger.info("导入draft_stat_detail")

    @report_check
    def add_draft_seq_detail(self, main_id, detail_path, samp, type, sof_type, status_dict=None):
        data_list =[]
        with open(detail_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split()
                data = {
                    "draft_id": ObjectId(main_id),
                    "seq_id": line[0],
                    "len": int(line[1]),
                    "gc_rate": float(line[2]),
                    "base_stat": line[3],
                    "samp": samp,
                    "type": type,
                    "sof_type": sof_type,
                }
                if status_dict:
                    data["status"] = status_dict[line[0]]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            self.db["draft_seq_detail"].insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("draft_seq_detail export %s error: %s" % (detail_path, e))
        else:
            self.bind_object.logger.info("导入draft_seq_detail信息成功")

    @report_check
    def add_draft_seq_cov(self, main_id, detail_path, samp, list_seqid):
        main_id = self.check_id(main_id)
        names = ["seq_id", "region", "cov"]
        data = pd.read_table(detail_path, header=None, names=names)
        data["cov_id"] = main_id
        data["samp"] = samp
        for i in data.seq_id.drop_duplicates():
            new_data = data.loc[data.seq_id==i,]
            new_data["sort"] = range(1, 101)
            data_list = new_data.to_dict(orient="records")
            try:
                self.db["seq_cov_detail"].insert_many(data_list)
            except Exception as e:
                self.bind_object.set_error("seq_cov_detail export %s error: %s: %s" % (detail_path, i, e))
        collection = self.db['seq_cov']
        collection.update({'_id': main_id}, {'$set': {'scaffolds': list_seqid}})
        self.bind_object.logger.info("导入seq_cov_detail信息成功")

    @report_check
    def add_draft_seq_bar(self, main_id, bar_file, samp, type, step):
        main_id = self.check_id(main_id)
        data = pd.read_table(bar_file, header=0, names=["len", "num", "prop_num", "sort"])
        data["draft_id"] = main_id
        data["step"] = step
        data["samp"] = samp
        data["type"] = type
        data_list = data.to_dict(orient="records")
        try:
            self.db["draft_seq_bar"].insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("draft_seq_bar export %s error: %s" % (bar_file, e))
        else:
            self.bind_object.logger.info("导入draft_seq_bar信息成功")

    @report_check
    def add_complete_seq_detail(self):
        pass

    @report_check
    def add_complete_stat_detail(self):

        pass


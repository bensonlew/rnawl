# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20180308
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId


class AnnoTcdb(Base):
    def __init__(self, bind_object):
        super(AnnoTcdb, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_anno_tcdb(self, anno_tcdb_path, main=False, task_id=None, project_sn=None, main_id=None,
                      params=None, name=None, update_id=None, specimen_id=None,analysis=None):
        self.specimen_id = specimen_id
        if main:
            if task_id is None:
                task_id = self.bind_object.sheet.id
            project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': 'TCDB注释',
                'created_ts': created_ts,
                'name': name if name else 'AnnoTcdb_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end',
                'version':'3.0',
                'settled_params': json.dumps({"version": "tcdb_v20200917"})  # by zhaozhigang 20200929
            }
            collection = self.db['anno_tcdb']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
            collection.update({'_id': main_id},{'$set':{'main_id':main_id}}) #guanqing 20180813
        else:
            if main_id is None:
                #raise Exception("main为False时需提供main_id!")
                self.bind_object.set_error("main为False时需提供main_id!", code="51401401")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        if analysis == 'complete':
            if specimen_id + '_whole_genome_tcdb_anno.xls' in os.listdir(anno_tcdb_path):    #guanqing.zou 20180912
                anno_file = anno_tcdb_path + '/' + specimen_id + '_whole_genome_tcdb_anno.xls'
            else:
                self.bind_object.set_error("找不到TCDB注释文件", code="51401402")
        else:
            if specimen_id + '_tcdb_anno.xls' in os.listdir(anno_tcdb_path):
                anno_file = anno_tcdb_path + '/' + specimen_id + '_tcdb_anno.xls'
            else:
                self.bind_object.set_error("找不到TCDB注释文件", code="51401402")

        self.add_anno_tcdb_detail(anno_file, main_id, update_id=update_id)
        return main_id

    def add_anno_tcdb_detail(self, detail_file, main_id=None, update_id=None):
        data_list = []
        collection_detail = self.db['anno_tcdb_detail']
        summary_collection = self.db['anno_summary_detail']
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                data = [("anno_tcdb_id", main_id)]
                line = line.strip().split('\t')
                data.extend(
                    [("gene_id", line[0]), ("location", line[1]), ("specimen_id", self.specimen_id), ("tcdb", line[3]),
                     ("des", line[4]), ("family", line[5]), ("subclass", line[6]), ("class", line[7]),
                     ("identity", float(line[8])), ("evalue", float(line[9])), ("score", float(line[10])),("coverage",float(line[12])*100)])
                data_son = SON(data)
                data_list.append(data_son)
                if update_id:
                    pass
                    # if not isinstance(update_id, ObjectId):
                    #     update_id = ObjectId(update_id)
                    # summary_collection.update_one(
                    #     {'summary_id': update_id, 'gene_id': line[0], 'specimen_id': self.specimen_id},
                    #     {'$set': {'tcdb_class': line[7]}},
                    #     upsert=False)
        try:
            collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入anno_tcdb_detail%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入anno_tcdb_detail%s信息出错:%s" , variables=(detail_file, e), code="51401403")
        else:
            self.bind_object.logger.info("导入anno_tcdb_detail%s信息成功!" % detail_file)
        collection_stat = self.db['anno_tcdb_stat']
        insert_data = {"anno_tcdb_id": main_id, "specimen_id": self.specimen_id, "anno_num": len(lines[1:])}
        try:
            collection_stat.insert_one(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入anno_tcdb_stat%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入anno_tcdb_stat%s信息出错:%s" , variables=(detail_file, e), code="51401404")
        else:
            self.bind_object.logger.info("导入anno_tcdb_stat%s信息成功!" % detail_file)

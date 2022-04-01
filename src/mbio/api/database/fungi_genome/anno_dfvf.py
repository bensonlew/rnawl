# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modify:20180526
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId


class AnnoDfvf(Base):
    def __init__(self, bind_object):
        super(AnnoDfvf, self).__init__(bind_object)
        self._project_type = "fungigenome"

    @report_check
    def add_anno_dfvf(self, anno_dfvf_path, main=False, task_id=None, project_sn=None, main_id=None,
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
                'desc': 'DFVF 注释',
                'created_ts': created_ts,
                'name': name if name else 'AnnoDfvf_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end'
            }
            collection = self.db['anno_dfvf']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
        else:
            if main_id is None:
                #raise Exception("main为False时需提供main_id!")
                self.bind_object.set_error("main为False时需提供main_id!", code="52100501")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        self.add_anno_dfvf_detail(anno_dfvf_path, main_id=main_id, update_id=update_id,specimen_id=specimen_id)
        return main_id

    def add_anno_dfvf_detail(self, detail_file_path, main_id=None, update_id=None, specimen_id=None):
        data_list = []
        collection_detail = self.db['anno_dfvf_detail']

        summary_collection = self.db['anno_summary_detail']
        detail_file = os.path.join(detail_file_path,"{}_dfvf.anno.xls".format(specimen_id))
        # gene_info = os.path.join(detail_file_path,"gene_info.dic")
        # fr = open(gene_info)
        # genedic = eval(fr.readline())
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                data = [("anno_dfvf_id", main_id)]
                line = line.strip().split('\t')
                tmp = line[0].split("__")
                scf = tmp[3].split("_")[0]
                data.extend(
                    [("gene_id", tmp[0]),("description",line[4]),("host",line[5]),("host_key",line[6]),
                     ("disease",line[7]),("disease_key",line[8]),("symbol",line[2]),("dfvf_name",line[1]),
                     ("identity", float(line[9])),("evalue",float(line[10])), ("score", float(line[11])),
                     ("specimen_id", self.specimen_id),("location",scf.capitalize())])
                     #("location",line[0].split("_")[0])])

                data_son = SON(data)
                data_list.append(data_son)
                if update_id:
                    if not isinstance(update_id, ObjectId):
                        update_id = ObjectId(update_id)
                    summary_collection.update_one(
                        {'summary_id': update_id, 'gene_id': tmp[0], 'specimen_id': self.specimen_id},
                        {'$set': {'dfvf_id': line[1]}},
                        upsert=False)
        try:
            collection_detail.insert_many(data_list)
            main_collection = self.db['anno_dfvf']
            main_collection.update({"_id": ObjectId(main_id)}, {"$set": {"status": 'end',"origin_id": ObjectId(main_id)}})

        except Exception as e:
            self.bind_object.logger.error("导入anno_dfvf_detail%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入anno_dfvf_detail%s信息出错:%s" , variables=(detail_file, e), code="52100502")
        else:
            self.bind_object.logger.info("导入anno_dfvf_detail%s信息成功!" % detail_file)
        collection_stat = self.db['anno_dfvf_stat']
        insert_data = {"anno_dfvf_id": main_id, "specimen_id": self.specimen_id, "anno_num": len(lines[1:])}
        try:
            collection_stat.insert_one(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入anno_dfvf_stat%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入anno_dfvf_stat%s信息出错:%s" , variables=(detail_file, e), code="52100503")
        else:
            self.bind_object.logger.info("导入anno_dfvf_stat%s信息成功!" % detail_file)


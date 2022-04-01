# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20181214


from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
import pandas as pd


class Assembly(Base):
    def __init__(self, bind_object):
        super(Assembly, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_assembly(self, genome_id, params=None, name=None):
        old_task_id = str(self.bind_object.sheet.id).split('_')[0]
        num = str(self.bind_object.sheet.id).split('_')[1]
        task_id = old_task_id + "_" + num
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "assembly主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "genome_id":genome_id,# G_bin_01
            "name": name if name else "Assembly_Origin",
            #"genome_path":"s3://metag_bin",
        }
        collection = self.db["assembly"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_assembly_detail(self, inserted_id, genome_id, insert_size, max_length, file_path=None, bin_cover_path1=None, assem_cover_path2=None,
                            genome_path=None, table_name=None, name=None, params=None, stat=None):
        if not isinstance(inserted_id, ObjectId):
            inserted_id = ObjectId(inserted_id)
        bin_cover = ''
        assem_cover = ''
        if bin_cover_path1:
            with open(bin_cover_path1, "r") as f_bin:
                lines = f_bin.readlines()
                bin_cover = int(lines[1].strip('\r\n').split('\t')[0])
        if assem_cover_path2:
            with open(assem_cover_path2, "r") as f_geno:
                lines = f_geno.readlines()
                assem_cover = int(lines[1].strip('\r\n').split('\t')[0])
        old_bin_id = genome_id.split('_')[1:]
        bin_id = ''.join(old_bin_id)
        data_list = []
        taxon = '-'
        if stat != None:
            with open(stat, "r") as m:
                line = m.readlines()
                taxon = line[1].strip('\r\n').split('\t')[-1]
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip('\r\n').split("\t")
                self.bind_object.logger.info("line长度%s".format(len(line)))
                data = {
                    "assemble_id": inserted_id,
                    "genome_id": genome_id,
                    "bin_id": bin_id,
                    "scaf_num":line[0],
                    "scaf_base": line[1],
                    "scaf_g_num": line[2],
                    "scaf_g_base": line[3],
                    "scaf_len_max": line[4],
                    "scaf_n50": line[5],
                    "scaf_n90": line[6],
                    "gc_rate": line[7],
                    "n_rate": line[8],
                    "binning_cover": bin_cover,
                    "genome_cover": assem_cover,
                    "taxon": taxon
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["assembly_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["assembly"]
            main_collection.update({"_id": ObjectId(inserted_id)},
                                   {"$set": {"status": 'end',
                                             'main_id': ObjectId(inserted_id),
                                             "insert": int(insert_size),
                                             "max": int(max_length),
                                             "genome_id": genome_id,
                                             "genome_path": genome_path,
                                             }})
            if name != None:
                main_collection = self.db["assembly"]
                main_collection.update({"_id": ObjectId(inserted_id)}, {"$set": {"name": name,
                                                                                 "status": 'end',
                                                                                 "params": json.dumps(params, sort_keys=True, separators=(',', ':'))}})
            if table_name != None:
                sg_status_collection = self.db["sg_status"]
                sg_status_collection.update({"table_id": ObjectId(inserted_id)},{"$set": {"status": 'end',
                                                                                          "desc": "Job has been finished",
                                                                                          "genome_id": genome_id}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

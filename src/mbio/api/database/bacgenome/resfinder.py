# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20201031
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from Bio import SeqIO
from bson.son import SON
from bson.objectid import ObjectId


class Resfinder(Base):
    def __init__(self, bind_object):
        super(Resfinder, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_resfinder(self, task_id=None, project_sn=None,params=None, name=None, specimen_id=None):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': '耐药基因ResFinder预测',
            'created_ts': created_ts,
            'name': name if name else 'Resfinder_Origin',
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end',
            "software" : "ResFinder_v4.1.0"
            # "specimen_id": specimen_id if specimen_id else ""
        }
        collection = self.db['resfinder']
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    @report_check
    def add_resfinder_dir(self, dir,sample, main_id=None):
        """
        导入详情表
        :param dir:文件夹路径
        :param sample: 样品名
        :param main_id: 主表id
        :return:
        """
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        else:
            main_id = main_id

        sample_stat = os.path.join(dir, sample + ".stat.xls")
        if os.path.exists(sample_stat):
            self.add_resfinder_stat(sample_stat, main_id=main_id)
        stat_file = os.path.join(dir, sample + "_resfinder.class.xls")
        if os.path.exists(stat_file):
            self.add_resfinder_type(stat_file, sample, "resfinder", main_id=main_id)
        stat_file2 = os.path.join(dir, sample + "_disinfinder.class.xls")
        if os.path.exists(stat_file2):
            self.add_resfinder_type(stat_file2, sample, "disinfinder", main_id=main_id)
        detail_file = os.path.join(dir, sample + "_resfinder.detail.xls")
        if os.path.exists(detail_file):
            self.add_resfiner_detail(detail_file, sample, "resfinder", main_id=main_id)
        detail_file2 = os.path.join(dir, sample + "_disinfinder.detail.xls")
        if os.path.exists(detail_file2):
            self.add_resfiner_detail(detail_file2, sample, "disinfinder", main_id=main_id)

    @report_check
    def add_resfinder_stat(self, detail_file, main_id=None):
        data_list = []
        collection_detail = self.db['resfinder_stat']
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = [("res_id", main_id),
                        ("specimen", line[0]),
                        ("gene_num", int(line[1]))]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入resfinder_stat%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入resfinder_stat%s信息出错" )
        else:
            self.bind_object.logger.info("导入resfinder_stat%s信息成功!" % detail_file)

    @report_check
    def add_resfinder_type(self, detail_file, sample,type, main_id=None):
        data_list = []
        collection_detail = self.db['resfinder_type']
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [("res_id", main_id),
                        ("class", line[1].strip()),
                        ("specimen", sample),
                        ("type", type),
                        ("gene_num", int(line[2]))]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            if len(data_list) > 0:
                collection_detail.insert_many(data_list)
                s_collection = self.db['resfinder']
                s_collection.update({'_id': main_id}, {'$set': {'main_id': main_id}})
        except Exception as e:
            self.bind_object.logger.error("导入resfinder_type%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入resfinder_type%s信息出错" )
        else:
            self.bind_object.logger.info("导入resfinder_type%s信息成功!" % detail_file)

    @report_check
    def add_resfiner_detail(self, detail_file, sample,type, main_id=None):
        data_list = []
        collection_detail = self.db['resfinder_detail']
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            data_count = len(lines) -1
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = [("res_id", main_id),
                        ("specimen", sample),
                        ("location", line[1]),
                        ("gene", line[0]),
                        ("resis_gene", line[3]),
                        ("class", line[4]),
                        ("pheno", line[5]),
                        ("resis", line[6]),
                        ("g_start", int(line[7])),
                        ("g_end", int(line[8])),
                        ("h_start", int(line[9])),
                        ("h_end", int(line[10])),
                        ("iden", float(line[11])),
                        ("cov", float(line[12])),
                        ("accession", line[13])]
                if type in ['resfinder']:
                    data.append(("database", "ResFinder"))
                else:
                    data.append(("database", "Disinfinder"))
                data_son = SON(data)
                data_list.append(data_son)
        try:
            if len(data_list) > 0:
                collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入resfinder_detail%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入resfinder_detail%s信息出错" )
        else:
            self.bind_object.logger.info("导入resfinder_detail%s信息成功!" % detail_file)

        try:
            collection = self.db["resfinder"]
            result = collection.find_one({"_id": main_id})
            if result:
                if 'specimen_id' in result:
                    specimen_id_list = str(result['specimen_id']).split(",")
                    if sample not in specimen_id_list:
                        specimen_id_list.append(sample)
                    new_specimen_id = ",".join(specimen_id_list)
                else:
                    new_specimen_id = sample
                if 'database_type' in result:
                    database_type_list = str(result['database_type']).split(",")
                    if "" in database_type_list:
                        database_type_list.remove("")
                    if type not in database_type_list:
                        if data_count != 0:
                            database_type_list.append(type)
                    database_type = ",".join(database_type_list)
                else:
                    if data_count != 0:
                        database_type = type
                    else:
                        database_type = ""
            else:
                new_specimen_id = sample
                if data_count != 0:
                    database_type = type
                else:
                    database_type = ""
            collection.update({"_id": main_id}, {"$set": {"specimen_id": new_specimen_id,
                                                          "database_type": database_type}})
        except Exception as e:
            self.bind_object.logger.error("更新主表resfinder主表信息出错:%s" % (e))
            self.bind_object.set_error("更新主表resfinder主表信息出错:%s" )
        else:
            self.bind_object.logger.info("更新主表resfinder主表信息成功!")








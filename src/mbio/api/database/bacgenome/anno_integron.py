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


class AnnoIntegron(Base):
    def __init__(self, bind_object):
        super(AnnoIntegron, self).__init__(bind_object)
        self._project_type = "bacgenome"
        self.seq_convert = {
            "Chromosome" : 'Chr',
            "Chromosome1" : "Chr1",
            "Chromosome2" : "Chr2",
            "Chromosome3" : "Chr3",
            "PlasmidA" :'pA',
            "PlasmidB" : 'pB',
            "PlasmidC" :'pC',
            "PlasmidD" : 'pD',
            "PlasmidE" :'pE',
            "PlasmidF" : 'pF',
            "PlasmidG" :'pG',
            "PlasmidH" : 'pH',
            "PlasmidI" : 'pI',
            "PlasmidJ" :'pJ',
            "PlasmidK" : 'pK',
            "PlasmidL" :'pL',
            "PlasmidM" : 'pM',
            "PlasmidN" : 'pN',
            "PlasmidO" :'pO',
            "PlasmidP" : 'pP',
        }

    @report_check
    def add_anno_integron(self, task_id=None, project_sn=None,params=None, name=None, specimen_id=None):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': '整合子分析',
            'created_ts': created_ts,
            'name': name if name else 'Integron_Origin',
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end',
            # 'specimen_id': specimen_id,
            "software" : "Integron_Finder_v2.0",
        }
        collection = self.db['integron']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': main_id},{'$set':{'main_id':main_id}})
        return main_id

    @report_check
    def add_anno_integron_stat(self, detail_file, main_id=None):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        data_list = []
        collection_detail = self.db['integron_stat']
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = [("integron_id", main_id),
                        ("sample", line[0]),
                        ("integron_num", int(line[1]))]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入integron_stat%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入integron_stat%s信息出错" )
        else:
            self.bind_object.logger.info("导入integron_stat%s信息成功!" % detail_file)

    @report_check
    def add_anno_integron_detail(self, detail_file, seq_file, main_id=None):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        data_list = []
        collection_detail = self.db['integron_detail']
        seq_dict = {}
        for seq_record in SeqIO.parse(seq_file, 'fasta'):
            id = seq_record.id
            seq = str(seq_record.seq)
            if id not in seq_dict:
                seq_dict[id] = seq
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [("integron_id", main_id),
                        ("sample", line[2]),
                        ("integron", line[0]),
                        ("strand", line[3]),
                        ("start", int(line[4])),
                        ("end", int(line[5])),
                        ("len", int(line[6])),
                        ("type", line[7]),
                        ("cds_num", int(line[8])),
                        ]
                if line[1] in self.seq_convert:
                    data.append(("location", line[1]))
                else:
                    data.append(("location", line[1]))
                if line[0] in seq_dict:
                    data.append(("seq", seq_dict[line[0]]))
                else:
                    data.append(("seq", "-"))
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection_detail.insert_many(data_list)
            s_collection = self.db['integron']
            s_collection.update({'_id': main_id}, {'$set': {'main_id': main_id}})
        except Exception as e:
            self.bind_object.logger.error("导入integron_detail%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入integron_detail%s信息出错" )
        else:
            self.bind_object.logger.info("导入integron_detail%s信息成功!" % detail_file)
        ## 下面代码更新circos表的字段，供前端判断使用
        task_id = "_".join(self.bind_object.sheet.id.split("_")[0:2])
        detail_table = self.db['circos_table_detail']
        try:
            type_main = self.db['integron']
            type_main_id = type_main.find_one({"_id": main_id})['_id']
            type_main_detail = self.db['integron_detail']
            main_table = self.db['circos_table']
            circos_table_id = main_table.find_one({"task_id": task_id})["_id"]
            if not isinstance(circos_table_id, ObjectId):
                circos_table_id = ObjectId(circos_table_id)
            for one in type_main_detail.find({"integron_id": type_main_id}):
                if one['location'].startswith('Scaffold'):
                    detail_table.update_one(
                        {"circos_table_id": circos_table_id, "specimen_id": one['sample']},
                        {'$set': {'Integron': 1}}, upsert=True)
                else:
                    detail_table.update_one(
                        {"circos_table_id": circos_table_id, "location": one['location'], "specimen_id": one['sample']},
                        {'$set': {'Integron': 1}}, upsert=True)
        except Exception as e:
            self.bind_object.logger.info('没有整合子信息')
        else:
            self.bind_object.logger.info('导入整合子信息成功')

    @report_check
    def add_anno_integron_summary(self, detail_file, seq_file, sample, main_id=None):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        seq_dict = {}
        for seq_record in SeqIO.parse(seq_file, 'fasta'):
            id = seq_record.id
            seq = str(seq_record.seq)
            if id not in seq_dict:
                seq_dict[id] = seq
        data_list = []
        collection_detail = self.db['integron_component']
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                if "ID_integron" in line:
                    continue
                data = [("integron_id", main_id),
                        ("sample", sample),
                        ("integron", "In" + line[0].split("_")[1]),
                        ("component", line[2]),
                        ("anno", line[8]),
                        ("start", int(line[3])),
                        ("end", int(line[4])),
                        ("evalue", line[6]),
                        ("model", line[9]),
                        ("distance", line[12]),
                        ("seq", line[1])]
                if line[1] in self.seq_convert:
                    data.append(("location", line[1]))
                else:
                    data.append(("location", line[1]))
                if line[5] in [1, "1"]:
                    data.append(("strand", '+'))
                    length = abs(int(int(line[4]) - int(line[3])))
                    data.append(("length", length))
                else:
                    data.append(("strand", '-'))
                    length = abs(int(int(line[3]) - int(line[4])))
                    data.append(("length", length))
                if line[2] in seq_dict:
                    data.append(("seq", seq_dict[line[2]]))
                else:
                    data.append(("seq", "-"))
                if line[8] in ['attC']:
                    data.append(("anno_type", "attC"))
                elif line[8] in ['protein']:
                    data.append(("anno_type", "tetm"))
                elif line[8] in ['intI']:
                    data.append(("anno_type", "Transposase"))
                else:
                    data.append(("anno_type", "tesC"))

                data_son = SON(data)
                data_list.append(data_son)

        try:
            collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入integron_component%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入integron_component%s信息出错" )
        else:
            self.bind_object.logger.info("导入integron_component%s信息成功!" % detail_file)

        try:
            collection = self.db["integron"]
            result = collection.find_one({"_id": main_id})
            if result:
                if 'specimen_id' in result:
                    specimen_id_list = str(result['specimen_id']).split(",")
                    if sample not in specimen_id_list:
                        specimen_id_list.append(sample)
                    new_specimen_id = ",".join(specimen_id_list)
                else:
                    new_specimen_id = sample
            else:
                new_specimen_id = sample
            collection.update({"_id": main_id}, {"$set": {"specimen_id": new_specimen_id}})
        except Exception as e:
            self.bind_object.logger.error("更新主表integron主表信息出错:%s" % (e))
            self.bind_object.set_error("更新主表integron主表信息出错:%s" )
        else:
            self.bind_object.logger.info("更新主表integron主表信息成功!")

    def add_main_id(self, main_id=None):
        try:
            s_collection = self.db['integron']
            s_collection.update({'_id': main_id}, {'$set': {'main_id': main_id}})
        except Exception as e:
            self.bind_object.logger.error("导入integron main_id出错:%s" % (e))
            self.bind_object.set_error("导入integron main_id出错" )





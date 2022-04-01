# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'
# last modify: 2020.10.26

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
from bson.objectid import ObjectId
from pymongo import MongoClient

class GenomeQuality(Base):
    def __init__(self, bind_object):
        super(GenomeQuality, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_genome_quality(self, main = None,params=None, name=None, task_id = None, project_sn = None, main_id=None):
        if main:
            if task_id is None:
                task_id = self.bind_object.sheet.id
            project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '基因组评估主表',
                'created_ts': created_ts,
                'name': name if name else 'GenomeQuality',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'start'
            }
            collection = self.db['genome_quality']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
            collection.update({'_id': main_id}, {'$set': {'main_id': main_id}})
        else:
            if main_id is None:
                # raise Exception("main为False时需提供main_id!")
                self.bind_object.set_error("main为False时需提供main_id!", code="51401101")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        return main_id

    def add_genome_quality_detail(self, task_id, type, project_sn = None, main_id=None, params=None, name=None):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        params = {"task_id":task_id}
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '基因组评估主表',
                'created_ts': created_ts,
                'name': name if name else 'GenomeQuality',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'start'
        }
        collection = self.db['genome_quality']
        # 将主表名称写在这里
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': main_id}, {'$set': {'main_id': main_id}})

        data_list = []
        collection_detail = self.db['genome_quality_detail']
        if type == "uncomplete":
            genome_size_main = self.db.datastat
            genome_size_detail = self.db.datastat_summary
            scaffold_num_main = self.db.assemble
            scaffold_num_detail = self.db.assemble_detail
            gene_num_main = self.db.gene_predict
            gene_num_detail = self.db.gene_predict_specimen
            genome_size_main_data = []
            genome_size_detail_data = []
            scaffold_num_main_data = []
            scaffold_num_detail_data = []
            gene_num_main_data = []
            gene_num_detail_data = []
            sample_id = []
            information = []

            # genome_size 提取
            genome_main_data = genome_size_main.find({"task_id": task_id})
            for i in genome_main_data:
                genome_size_main_data.append(i)
            genome_detail = genome_size_detail.find({"datastat_id": genome_size_main_data[0]["main_id"]})
            for i in genome_detail:
                genome_size_detail_data.append(i)
            if len(genome_size_detail_data) == 0:
                size = 0
            else:
                size = 1
            # scaffold数量提取
            scaffold_main = scaffold_num_main.find({"task_id": task_id})
            for i in scaffold_main:
                scaffold_num_main_data.append(i)
            scaffold_detail = scaffold_num_detail.find({"assemble_id": scaffold_num_main_data[-1]["main_id"]})
            for i in scaffold_detail:
                scaffold_num_detail_data.append(i)
            for i in range(len(scaffold_num_detail_data)):
                if scaffold_num_detail_data[i]["specimen_id"] not in sample_id:
                    sample_id.append(scaffold_num_detail_data[i]["specimen_id"])
            # gene 数量提取
            gene_main = gene_num_main.find({"task_id": task_id})
            for i in gene_main:
                gene_num_main_data.append(i)
            gene_detail_data = gene_num_detail.find({"predict_id": gene_num_main_data[-1]["main_id"]})
            for i in gene_detail_data:
                gene_num_detail_data.append(i)
            for i in sample_id:
                if size:
                    for a in range(len(genome_size_detail_data)):
                        if genome_size_detail_data[a]["specimen_id"] == i:
                            information.append(i)
                            information.append(str(float(genome_size_detail_data[a]["g_size"])/(1000*1000)))
                            break
                else:
                    information.append(i)
                    information.append("-")
                for b in range(len(scaffold_num_detail_data)):
                    if scaffold_num_detail_data[b]["specimen_id"] == i:
                        information.append(scaffold_num_detail_data[b]["scaf_no"])
                        break
                for c in range(len(gene_num_detail_data)):
                    if gene_num_detail_data[c]["specimen_id"] == i:
                        information.append(gene_num_detail_data[c]["num"])
                        information.append(gene_num_detail_data[c]["aver_base"])
                        break
                if (1 <= float(information[1]) <= 10) and (1 <= int(information[2]) <= 100) and (
                        1000 <= float(information[3]) <= 10000) and (800 <= float(information[4]) <= 1200):
                    information.append("yes")
                else:
                    information.append("no")
                data = {
                    "quality_id": main_id,
                    "specimen_id": information[0],
                    "genome_size": information[1],
                    "scaffold_num": information[2],
                    "cds_num": information[3],
                    "cds_average_len": information[4],
                    "normal": information[5]
                    }
                information = []
                data_son = SON(data)
                data_list.append(data_son)
        else:
            genome_size_main = self.db.datastat
            genome_size_detail = self.db.datastat_summary
            chr_num_main = self.db.assemble
            chr_num_detail = self.db.assemble_detail
            gene_num_main = self.db.gene_predict
            gene_num_detail = self.db.gene_predict_specimen
            gene0001_main = self.db.anno_summary
            gene0001_detail = self.db.anno_summary_detail
            genome_size_main_data = []
            genome_size_detail_data = []
            chr_num_main_data = []
            chr_num_detail_data = []
            gene_num_main_data = []
            gene_num_detail_data = []
            sample_id = []
            gene0001_data = []
            print(task_id)
            # genome_size 提取
            genome_main_data = genome_size_main.find({"task_id": task_id})
            for i in genome_main_data:
                genome_size_main_data.append(i)
            if genome_size_main_data:
                genome_detail = genome_size_detail.find({"datastat_id": genome_size_main_data[0]["main_id"]})
                for i in genome_detail:
                    genome_size_detail_data.append(i)
            if len(genome_size_detail_data) == 0:
                size = 0
            else:
                size = 1
            # chr和plas数量提取
            chr_main = chr_num_main.find({"task_id": task_id})
            for i in chr_main:
                chr_num_main_data.append(i)
            chr_detail = chr_num_detail.find({"assemble_id": chr_num_main_data[-1]["main_id"]})
            for i in chr_detail:
                chr_num_detail_data.append(i)
            for i in range(len(chr_num_detail_data)):
                if chr_num_detail_data[i]["specimen_id"] not in sample_id:
                    sample_id.append(chr_num_detail_data[i]["specimen_id"])
            # gene 数量提取
            gene_main = gene_num_main.find({"task_id": task_id})
            for i in gene_main:
                gene_num_main_data.append(i)
            gene_detail_data = gene_num_detail.find({"predict_id": gene_num_main_data[-1]["main_id"]})
            for i in gene_detail_data:
                gene_num_detail_data.append(i)
            for i in sample_id:
                information = []
                if size:
                    for a in range(len(genome_size_detail_data)):
                        if genome_size_detail_data[a]["specimen_id"] == i:
                            information.append(i)
                            information.append(str(float(genome_size_detail_data[a]["g_size"])/(1000*1000)))
                            break
                else:
                    information.append(i)
                    information.append("-")
                for b in range(len(chr_num_detail_data)):
                    if chr_num_detail_data[b]["specimen_id"] == i:
                        information.append(chr_num_detail_data[b]["chr_no"])
                        information.append(chr_num_detail_data[b]["pla_no"])
                        break
                for c in range(len(gene_num_detail_data)):
                    if gene_num_detail_data[c]["specimen_id"] == i:
                        information.append(gene_num_detail_data[c]["num"])
                        information.append(gene_num_detail_data[c]["aver_base"])
                        break
                result = list(gene0001_main.find({"task_id": task_id}))
                main_id1 = result[0]['main_id']
                gene0001tmp = list(gene0001_detail.find(
                    {"summary_id": main_id1, "gene_id": "gene0001", "specimen_id": i}))
                if len(gene0001tmp) >0:
                    if gene0001tmp[0]["gene_name"] == "dnaA":
                        information.append("yes")
                    else:
                        information.append("no")
                else:
                    information.append("no")
                if size:
                    if 1 <= float(information[1]) <= 10 and (0 <= int(information[2]) <= 3) and (
                            0 <= int(information[3]) <= 10) and (1000 <= float(information[4]) <= 10000) and (
                            800 <= float(information[5]) <= 1200) and information[6] == "yes":
                        information.append("yes")
                    else:
                        information.append("no")
                else:
                    if (str(information[1]) == "-" and (0 <= int(information[2]) <= 3) and (
                            0 <= int(information[3]) <= 10) and (1000 <= float(information[4]) <= 10000) and (
                            800 <= float(information[5]) <= 1200) and information[6] == "yes"):
                        information.append("yes")
                    else:
                        information.append("no")
                data = {
                    "quality_id": main_id,
                    "specimen_id": information[0],
                    "genome_size": information[1],
                    "chromosome_num": information[2],
                    "plasmid_num": information[3],
                    "cds_num": information[4],
                    "cds_average_len": information[5],
                    "first_gene": information[6],
                    "normal": information[7]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection_detail.insert_many(data_list)
            s_collection = self.db['genome_quality']
            s_collection.update({'_id': main_id}, {'$set': {'main_id': main_id, 'status': 'end'}})
        except Exception as e:
            print("error")
            self.bind_object.logger.error("导入genome_quality_detail%s信息出错:%s" % (task_id, e))
        else:
            print("ok")
            self.bind_object.logger.info("导入genome_quality_detail%s信息成功!" % task_id)

if __name__ == "__main__":
    a = GenomeQuality(None)
    #task_id = "majorbio_236401"
    #project_sn = "test01"
    #a.add_genome_quality(main="123",project_sn="234",task_id="tsg_36243",params="123456")
    #a.add_genome_quality_detail("tsg_37312", "uncomplete")
    #a.add_genome_quality_detail("tsg_34765", "complete",main_id= insert_id)

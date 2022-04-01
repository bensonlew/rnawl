# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
#20190902

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON


class GeneStat(Base):
    def __init__(self, bind_object):
        super(GeneStat, self).__init__(bind_object)
        self._project_type = "bac_comparative"

    @report_check
    def add_gene_stat(self, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "基因预测主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": "GeneStat_Origin",
        }
        collection = self.db["gene_stat"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    def add_gene_stat_detail(self, inserted_id, specimen_id, stat):
        with open(stat, "r") as f:
            lines = f.readlines()
            data_list = []
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "stat_id": ObjectId(inserted_id),
                    "specimen_id": specimen_id,
                    "cds": int(line[1]),
                    "rrna": int(line[2]),
                    "s5": int(line[5]),
                    "s23": int(line[4]),
                    "s16": int(line[3]),
                    "trna": int(line[6])
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["gene_stat_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (stat, e))

    def add_genomes_detail(self, specimen_id, stat):
        """
        用于更新基因组基本特征详情表
        inserted_id：data表主表id；
        specimen_id：样品名称；
        stat：每个样本对应的统计表
        """
        task_id = self.bind_object.sheet.id
        data = self.db['data']
        data_id = data.find_one({"task_id": task_id})['_id']
        self.bind_object.logger.info(data_id)
        with open(stat, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                cds = int(line[1])
                rrna = int(line[2])
                trna = int(line[6])
                s16 = int(line[3])
                try:
                    collection = self.db["genomes_detail"]
                    collection.update({"data_id": ObjectId(data_id), "specimen_id": specimen_id}, {'$set': {"cds": cds, "rrna": rrna, "trna": trna, "s16_rna": s16}})
                except Exception, e:
                    self.bind_object.logger.info("更新%s结果表出错:%s" % (specimen_id, e))
                else:
                    self.bind_object.logger.info("更新%s结果表成功!" % specimen_id)

    def add_data_detail(self, stat):
        """
        用于更新基因组基本特征详情表
        stat：每个样本对应的统计表
        Assembly_accession	Total_length	Genome_length	GC	Chrom_num
        Plasmid_num	Scaffold_num	Contig_num	Scaffold_n50	Contig_n50	Seq_status	Species	Gene_prefix
        """
        data_list = []
        task_id = self.bind_object.sheet.id
        data = self.db['data']
        data_id = data.find({"task_id": task_id})['_id']
        self.bind_object.logger.info(data_id)
        with open(stat, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                assembly = line[0]
                g_size = str(line[1])
                location = line[2]
                g_location = line[3]
                seq_size = str(line[4])
                gc = str(line[5])
                chr_no = str(line[6])
                pla_no = str(line[7])
                scf_no = str(line[8])
                con_no = str(line[9])
                scf_n50 = str(line[10])
                con_n50 = str(line[11])
                seq_status = line[12]
                species = line[13]
                gene_name = line[14]
                data = {
                    "data_id": ObjectId(data_id),
                    "ass_assession": assembly,
                    "specimen_id": assembly,
                    "g_size": g_size,
                    "location": location,
                    "g_location": g_location,
                    "seq_size": seq_size,
                    "gc": gc,
                    "chr_no": chr_no,
                    "pla_no": pla_no,
                    "scf_no": scf_no,
                    "con_no": con_no,
                    "scf_n50": scf_n50,
                    "con_n50": con_n50,
                    "seq_status": seq_status,
                    "species": species,
                    "gene_name": gene_name,
                    "from_refdb": "F"
                }
                data_son = SON(data)
                data_list.append(data_son)
            try:
                collection = self.db["data_detail"]
                collection.insert_many(data_list)
            except:
                self.bind_object.logger.info("更新data_detail结果表出错!")
            else:
                self.bind_object.logger.info("更新data_detail结果表成功!")


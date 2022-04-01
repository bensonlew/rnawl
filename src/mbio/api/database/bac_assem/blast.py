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
import random
from Bio import SeqIO


class Blast(Base):
    def __init__(self, bind_object):
        super(Blast, self).__init__(bind_object)
        self._project_type = "bac_assem"

    @report_check
    def add_blast_nt(self, params, task_id=None, project_sn=None, name=None):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': 'blast_nt',
            'created_ts': created_ts,
            'name': name if name else 'origin',
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end'
        }
        collection = self.db['blast_nt']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': main_id},{'$set': {'main_id': main_id}})
        return main_id

    def check_id(self, main_id):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！')
        return main_id

    def export_db(self, db_name, data_list, table):
        try:
            self.db[db_name].insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("%s export %s error %s" % (db_name, table, e))
        else:
            self.bind_object.logger.info("%s export success" % db_name)

    @report_check
    def add_blast_nt_detail(self, table_name, main_id, table, list_seqid, genome_id=None):
        main_id = self.check_id(main_id)
        if genome_id:
            genome_id = self.check_id(genome_id)
        index = ["samp", "seq_id", "aln_len", "q_len", "q_start", "q_end", "q_direction", "sbj_id", "s_len", "s_start", "s_end", "s_direction", "desc", "ident", "cov", "evalue", "score"]
        data = pd.read_table(table, header=0, names=index)
        data["blast_id"] = main_id
        data_list = data.to_dict(orient="records")
        if not data.empty:
            self.export_db(table_name, data_list, table)
        if genome_id:
            self.db["sg_status"].update({"table_id": main_id}, {"$set": {"genome_id": genome_id}})
        collection = self.db['blast_nt']
        collection.update({'_id': main_id}, {'$set': {'scaffolds': list_seqid}})

    @report_check
    def add_compblast_nt_detail(self, table_name, main_id, table, list_seqid):
        main_id = self.check_id(main_id)
        index = ["samp", "seq_id", "aln_len", "q_len", "q_start", "q_end", "q_direction", "sbj_id", "s_len", "s_start",
                 "s_end", "s_direction", "desc", "ident", "cov", "evalue", "score"]
        data = pd.read_table(table, header=0, names=index)
        data["blast_id"] = main_id
        data_list = data.to_dict(orient="records")
        if not data.empty:
            self.export_db(table_name, data_list, table)
        collection = self.db['complete_blastnt']
        collection.update({'_id': main_id}, {'$set': {'scaffolds': list_seqid}})

    @report_check
    def add_gap_fill(self, params, task_id=None, project_sn=None, name=None, assem_para=None):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': 'gap_fill',
            'created_ts': created_ts,
            'name': name if name else "origin",
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end',
        }
        collection = self.db['gap_fill']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': main_id},{'$set':{'main_id':main_id}})
        return main_id

    @report_check
    def add_first_gap_fill_detail(self, main_id, sample, data, seq_path, assem_para):
        main_id = self.check_id(main_id)
        collection = self.db["gap_fill"]
        collection.update({"_id": main_id}, {"$set": {"sample": sample, "seq_path": seq_path, "assem_para": assem_para}})
        data = data.reindex(columns=["samp", "seq_id", "q_len", "chrome_info"])
        data.columns = ["samp", "scaf_name", "length", "scaf_status"]
        data["gap_id"] = main_id
        data["run_log"] = data["scaf_name"]
        data["location"] = "-"
        data_list = data.to_dict(orient="records")
        self.export_db("gap_fill_detail", data_list, "gap_fill_detail")

    @report_check
    def add_gap_fill_detail(self, main_id, sample, table, seq_path, sw=None, status_dict=None):
        main_id = self.check_id(main_id)
        collection = self.db["gap_fill"]
        collection.update({"_id": main_id}, {"$set": {"sample": sample, "seq_path": seq_path}})
        if status_dict:
            data = pd.read_table(table, header=None, names=["scaf_name", "length", "gc", "base"])
            data = data.loc[:,["scaf_name","length"]]
            data["run_log"] = data["scaf_name"]
            data.loc[data["length"] >= 1000000, "location"] = "Chromosome"
            data.loc[data["length"] < 1000000, "location"] = "Plasmid"
        else:
            data = pd.read_table(table, header=None, names=["location", "scaf_name", "run_log", "scaf_status", "length"])
        data["samp"] = sample
        data["gap_id"] = main_id
        if sw: # 交互的情况
            data["run_log"] = data["run_log"].str.replace(",", "\n") + "\n" + sw
        data_list = data.to_dict(orient="records")
        self.export_db("gap_fill_detail", data_list, table)

    @report_check
    def add_gap_fill2_detail(self, main_id, sample, table, seq_path, status_dict=None):
        main_id = self.check_id(main_id)
        collection = self.db["gap_fill"]
        collection.update({"_id": main_id}, {"$set": {"sample": sample, "seq_path": seq_path}})
        list = 'abcdefghijklmnopurstuvwxyz'
        list2 = '0abcdefghijklmnopurstuvwxyz'
        dict = {}
        n = 1
        for i in list2:
            for j in list:
                if i in ['0']:
                    dict[n] = j.upper()
                else:
                    dict[n] = i.upper() + j.upper()
                n += 1
        data_list = []
        i =1
        j =1
        with open(table, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split()
                data = {
                    "gap_id": ObjectId(main_id),
                    "scaf_name": line[0],
                    "length": int(line[1]),
                    "samp": sample,
                    "run_log": line[0],
                }
                if status_dict:
                    data["scaf_status"] = status_dict[line[0]]
                if int(line[1]) >= 1000000:
                    data['location'] = "Chromosome" + str(j)
                    j += 1
                else:
                    data['location'] = "Plasmid" + dict[i]
                    i += 1
                data_son = SON(data)
                data_list.append(data_son)
        try:
            self.db["gap_fill_detail"].insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("gap_fill_detail export %s error: %s" % (table, e))
        else:
            self.bind_object.logger.info("导入gap_fill_detail信息成功")

    @report_check
    def add_complete_stat_detail(self, main_id, genome_id, samp, seq_path, table1, table2):
        main_id = self.check_id(main_id)
        genome_id = self.check_id(genome_id)
        data = pd.read_table(table2, header=None, names=["seq_id", "len", "gc_rate", "seq_from"])
        data["complete_id"] = main_id
        data["samp"] = samp
        data.loc[data["seq_id"].str.contains("Chromosome"),"type"] = "Chromosome"
        data.loc[data["seq_id"].str.contains("Plasmid"),"type"] = "Plasmid"
        data_list = data.reindex(columns=["complete_id", "samp", "seq_id", "type", "len", "gc_rate", "seq_from"]).to_dict(orient="records")
        with open(table1, "r") as file:
            line = file.readlines()
            print(line[1])
            line = line[1].strip().split("\t")
            print(line[0])
            main_data = {
                "chr_no": int(line[0]),
                "pls_no": int(line[1]),
                "genome_size": int(line[2]),
                "gc_rate": float(line[3]),
                "seq_path": seq_path
            }
        self.db["complete_stat"].update({"_id": main_id}, {"$set": main_data})
        self.export_db("complete_stat_detail", data_list, table2)
        self.db["sg_status"].update({"table_id": main_id}, {"$set": {"genome_id": genome_id}})

    @report_check
    def add_kmer_pca_detail(self, table_name, main_id,samp, table, table_importance, genome_id=None):
        main_id = self.check_id(main_id)
        if genome_id:
            genome_id = self.check_id(genome_id)
        data = pd.read_table(table, header=0)
        data = data.reindex(columns=["Sample_ID", "PC1", "PC2"])
        data.columns = ["label", "pc1", "pc2"]
        data["samp"] = samp
        data["kmer_pca_id"] = main_id
        data["group"] = data.apply(lambda x: x["label"].split(".")[0], axis=1)
        data_list = data.to_dict(orient="records")
        self.export_db(table_name, data_list, table)
        data2 = pd.read_table(table_importance, header=0)
        importance = [data2.iloc[0,1], data2.iloc[1,1]]
        if genome_id:
            self.db["complete_kmer_pca"].update({"_id": main_id}, {'$set': {"group": data["group"].drop_duplicates().tolist(), "importance": importance}})
            self.db["sg_status"].update({"table_id": main_id}, {"$set": {"genome_id": genome_id}})
        else:
            self.db["kmer_pca"].update({"_id": main_id}, {'$set': {"group": data["group"].drop_duplicates().tolist(), "importance": importance}})

    @report_check
    def add_complete_organism_hk(self, table_name, main_id, genome_id, samp, table, listname):
        main_id = self.check_id(main_id)
        genome_id = self.check_id(genome_id)
        data = pd.read_table(table, header=0, names=["house", "org", "ident", "cov", "evalue", "score"])
        if data.empty:
            return
        data["samp"] = samp
        data["organism_id"] = main_id
        data_list = data.to_dict(orient="records")
        self.export_db(table_name, data_list, "hgenetable")
        collection = self.db['complete_organism']
        collection.update({'_id': main_id}, {'$set': {'names': listname}})
        self.db["sg_status"].update({"table_id": main_id}, {"$set": {"genome_id": genome_id}})



    @report_check
    def add_complete_organism_16s(self, table_name, main_id, genome_id, samp, table):
        main_id = self.check_id(main_id)
        genome_id = self.check_id(genome_id)
        data = pd.read_table(table, header=0, names=["org", "ident", "cov", "evalue", "score"])
        if data.empty:
            return
        data["samp"] = samp
        data["organism_id"] = main_id
        data_list = data.to_dict(orient="records")
        self.export_db(table_name, data_list, "16stable")
        self.db["sg_status"].update({"table_id": main_id}, {"$set": {"genome_id": genome_id}})

    @report_check
    def add_kmer_dist(self, params, task_id=None, project_sn=None, name=None, main_id=None):
        main_id = self.check_id(main_id)
        if task_id is None:
            task_id = self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': 'kmer_dist',
            'created_ts': created_ts,
            'name': name if name else "Origin",
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end',
            "main_id": main_id
        }
        collection = self.db['kmer_dist']
        collection.insert_one(insert_data)

    @report_check
    def add_gc_depth(self, params, task_id=None, project_sn=None, name=None, main_id=None):
        main_id = self.check_id(main_id)
        if task_id is None:
            task_id = self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': 'gc_depth',
            'created_ts': created_ts,
            'name': name if name else "Origin",
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end',
            "main_id": main_id
        }
        collection = self.db['gc_depth']
        collection.insert_one(insert_data)

    def update_task_id(self):
        for coll in ["gap_fill", "complete_stat", "kmer_dist", "gc_depth"]:
            self.db[coll].update({"task_id": "1", "project_sn": "1"}, {'$set': {"fake": "1", "task_id": self.bind_object.sheet.id, "project_sn": self.bind_object.sheet.project_sn}})

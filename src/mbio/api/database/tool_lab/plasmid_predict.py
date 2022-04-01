# -*- coding: utf-8 -*-

from bson.objectid import ObjectId
import datetime
import os, re
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase

class PlasmidPredict(ApiBase):
    def __init__(self, bind_object):
        super(PlasmidPredict, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_plasmid_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "plasmid_predict"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='plasmid_predict',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('plasmid_predict', [main_info])
        else:
            main_id = ObjectId(main_id)
            table1_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name","type": "string"},
                    {"field": "contig_num", "filter": "false", "sort": "false", "title": "Contig num","type": "int"},
                    {"field": "chromosome_num", "filter": "false", "sort": "false", "title": "Chromosome num","type": "int"},
                    {"field": "plasmid_num", "filter": "false", "sort": "false", "title": "Plasmid num","type": "int"},
                    {"field": "unclassified_num", "filter": "false", "sort": "false", "title": "Unclassified num","type": "int"}],
                "condition": {}}
            table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
            table2_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
                    {"field": "contig_name", "filter": "false", "sort": "false", "title": "Contig name", "type": "string"},
                    {"field": "plasmid_type", "filter": "false", "sort": "false", "title": "Plasmid type", "type": "string"},
                    {"field": "probability", "filter": "false", "sort": "false", "title": "Probability(%)", "type": "float"},
                    {"field": "length", "filter": "false", "sort": "false", "title": "Contig length (bp)", "type": "int"},
                    {"field": "gc_content", "filter": "false", "sort": "false", "title": "GC content (%)", "type": "float"},
                    {"field": "acc_num", "filter": "false", "sort": "false", "title": "Acc num", "type": "string"},
                    {"field": "acc_name", "filter": "false", "sort": "false", "title": "Acc plasmid name", "type": "string"},
                    {"field": "identity", "filter": "false", "sort": "false", "title": "Identity (%)", "type": "float"},
                    {"field": "evalue", "filter": "false", "sort": "false", "title": "E-value", "type": "float"},
                    {"field": "taxon", "filter": "false", "sort": "false", "title": "Taxon name", "type": "string"}],
                "condition": {}}
            table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
            try:
                self.update_db_record('plasmid_predict', main_id, status="end", main_id=main_id,
                                      stat_table=table1_info, detail_table=table2_info)
            except Exception as e:
                self.bind_object.logger.error("导入plasmid_predict数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入plasmid_predict数据成功")
        return main_id

    def add_plasmid_detail(self, main_id, file_path, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "prokka"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='prokka',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('prokka', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data1 = [];
        insert_data2 = [];
        all_plas = {};
        base = 0;
        gc = 0;
        blast_result = {};
        plsdb = {}
        self.get_name(file_path)
        with open(file_path + '/' + self.sample_name +"_plasflow.log") as f:
            data1 = f.readlines()
            for i in range(len(data1)):
                if data1[i].startswith("Resulting plasmid sequences"):
                    #if ("chromosome" in data1[i].strip().split()) and ("unclassified" in data1[i].strip().split()) and ("plasmid" in data1[i].strip().split()):
                    tmp_dir ={}
                    for kj in range(len(data1[i+1].strip().split())):
                        tmp_dir[data1[i+1].strip().split()[kj]] = data1[i + 2].strip().split()[kj + 1]
                    if "chromosome" in tmp_dir:
                        chromosome = int(tmp_dir["chromosome"])
                    else:
                        chromosome = 0
                    if "plasmid" in tmp_dir:
                        plasmid = int(tmp_dir["plasmid"])
                    else:
                        plasmid = 0
                    if "unclassified" in tmp_dir:
                        unclassified = int(tmp_dir["unclassified"])
                    else:
                        unclassified = 0
                    contig_num = chromosome + plasmid + unclassified
                    insert_data1.append({
                        "contig_num": contig_num,
                        "chromosome_num": chromosome,
                        "plasmid_num": plasmid,
                        "unclassified_num": unclassified,
                        "sample_name": self.sample_name,
                        "plasmid_id": main_id
                    })
                    break
        file2 = file_path + '/' + self.sample_name + "_plasflow_predictions.tsv_plasmids.fasta"
        if os.path.getsize(file2) != 0:
            with open(file2) as jk:
                data2 = jk.read()
                for i in range(len(data2.split(">"))):
                    if len(data2.split(">")[i].split("\n")) > 1:
                        all_plas[data2.split(">")[i].split("\n")[0].split()[0]] = []
                        for s in data2.split(">")[i].split("\n")[1:]:
                            base += len(s)
                            gc += (s.count('g') + s.count('G') + s.count('c') + s.count('C'))
                        gc_ratio = float(gc) / base
                        all_plas[data2.split(">")[i].split("\n")[0].split()[0]].append(base)
                        all_plas[data2.split(">")[i].split("\n")[0].split()[0]].append(gc_ratio)
                        gc = 0
                        base = 0
                self.bind_object.logger.info(all_plas)
            file1 = file_path + '/' + self.sample_name + "_plasflow_predictions.tsv"
            with open(file1) as lg:
                data1 = lg.readlines()
                for n in data1[1:]:
                    if n.split("\t")[2] in all_plas.keys():
                        # all_plas[n.split("\t")[2]].append(n.split("\t")[5].split(".")[1])
                        for z in range(len(data1[0].strip().split("\t"))):
                            if n.split("\t")[5] == data1[0].strip().split("\t")[z]:
                                all_plas[n.split("\t")[2]].append(n.split("\t")[z+1])
                        if len(all_plas[n.split("\t")[2]]) <3:
                            all_plas[n.split("\t")[2]].append("-")
                self.bind_object.logger.info(all_plas)
            file3 = file_path + '/' + self.sample_name + "_plasmid_blast.txt"
            file4 = file_path + '/' + self.sample_name + "_plasmid_plsdb.tsv"
            if os.path.getsize(file3) != 0:
                with open(file3) as ea, open(file4) as ca:
                    data3 = ea.readlines()
                    data4 = ca.readlines()
                    for b in data4:
                        plsdb[b.split("\t")[1]] = b
                    for i in data3:
                        if i.split("\t")[0] in blast_result.keys():
                            pass
                        else:
                            blast_result[i.split("\t")[0]] = [i.split("\t")[1]]
                            blast_result[i.split("\t")[0]].append(i.split("\t")[2])
                            blast_result[i.split("\t")[0]].append(i.split("\t")[10])
                    for key in all_plas.keys():
                        if key in blast_result.keys():
                            tmp_data = plsdb[blast_result[key][0]]
                            all_plas[key].append(blast_result[key][0])
                            all_plas[key].append(blast_result[key][1])
                            all_plas[key].append(blast_result[key][2])
                            if "plasmid" in tmp_data.split("\t")[2]:
                                all_plas[key].append(tmp_data.split("\t")[2].split("plasmid")[1].strip().split(",")[0])
                            else:
                                all_plas[key].append("-")
                            if tmp_data.split("\t")[24]:
                                all_plas[key].append(tmp_data.split("\t")[24])
                            else:
                                all_plas[key].append("-")
                            if tmp_data.split("\t")[49]:
                                all_plas[key].append(tmp_data.split("\t")[49].split(",")[0].strip())
                            else:
                                all_plas[key].append("-")
                        else:
                            all_plas[key].append("-")
                            all_plas[key].append("-")
                            all_plas[key].append("-")
                            all_plas[key].append("-")
                            all_plas[key].append("-")
                            all_plas[key].append("-")
            else:
                for key in all_plas.keys():
                    all_plas[key].append("-")
                    all_plas[key].append("-")
                    all_plas[key].append("-")
                    all_plas[key].append("-")
                    all_plas[key].append("-")
                    all_plas[key].append("-")
            for bj in all_plas.keys():
                insert_data2.append({
                    "sample_name": self.sample_name,
                    "contig_name": bj,
                    "plasmid_type": all_plas[bj][8],
                    "probability": "-" if all_plas[bj][2] == "-" else round(float(all_plas[bj][2])*100, 2),
                    "length": int(all_plas[bj][0]),
                    "gc_content": round(float(all_plas[bj][1])*100, 2),
                    "acc_num": all_plas[bj][3],
                    "acc_name": all_plas[bj][6],
                    "identity": "-" if all_plas[bj][4] == "-" else round(float(all_plas[bj][4]), 3),
                    "evalue": all_plas[bj][5],
                    "taxon": all_plas[bj][7],
                    "plasmid_id": main_id,
                })
        #self.bind_object.logger.info(insert_data1)
        #self.bind_object.logger.info(insert_data2)
        try:
            self.create_db_table('plasmid_predict_stat', insert_data1)
            self.create_db_table('plasmid_predict_detail', insert_data2)
        except Exception as e:
            self.bind_object.logger.error("导入plasmid_predict数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入plasmid_predict数据成功")

    def get_name(self,path):
        for file in os.listdir(path):
            if file.endswith("plasflow_predictions.tsv"):
                self.sample_name = file.split("_plasflow_predictions.tsv")[0]
# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'
# last modify: 2020.10.22

from biocluster.api.database.base import Base, report_check
import datetime
import json
import re
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId


class Repeatmasker(Base):
    def __init__(self, bind_object):
        super(Repeatmasker, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_repeat_predict(self, params=None):
        task_id = self.bind_object.sheet.id
        #print task_id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "start",
            "desc": "散在重复序列预测主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": "Repeatmasker",
            "software": "RepeatMasker_4.0.7"
        }
        collection = self.db["interpersed_repeat"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    @report_check
    def get_id(self, task_id):
        collection = self.db["interpersed_repeat"]
        id = list(collection.find({"task_id":task_id},sort= [('id',-1)]))[-1]['_id']
        collection.update({'_id': id}, {'$set': {'main_id': id}})
        return id

    @report_check
    def add_repeat_predict_detail(self, inserted_id, specimen_id, predict_gff ,fa_path):
        gff = pd.read_table(predict_gff, sep='\t', header=0)
        h = open(fa_path)
        seq_all = h.read()
        seq_data = seq_all.split(">")
        if len(gff) < 1:
            return
        num = len(gff)
        total_base = (gff["End"].sum() - gff["Begin"].sum() + num -1)
        main_collection = self.db["interpersed_repeat"]
        task_id = main_collection.find_one({"_id": ObjectId(inserted_id)})['task_id']
        main_gene = self.db["gene_predict"]
        #print task_id
        #cds_id = main_gene.find_one({"task_id": task_id})['main_id']
        pre_all = main_gene.find({"task_id": task_id})
        detail_gene = self.db["gene_predict_specimen"]
        genome_total_base = ''
        for g_p in pre_all:
            cds_id = g_p['main_id']
            g_p_detail = detail_gene.find_one({"specimen_id": specimen_id, "predict_id": cds_id})
            print ObjectId(cds_id)
            if g_p_detail:
                genome_total_base = g_p_detail['total_base']
                break
            print genome_total_base
        # self.bind_object.logger.info(genome_total_base)
        if genome_total_base not in ['','-']:
            in_genome = float(total_base) / float(genome_total_base)
        else:
            in_genome = '-'
        gff["location"] = gff["Sequence Name"]
        gff["match"] = (gff["Attributes"].str.split(";", expand=True)[1]).str.split(" ", expand=True)[0]
        gff["TE"] = (gff["Attributes"].str.split(";", expand=True)[0]).str.split("=", expand=True)[1]
        gff["class"] = (gff["Attributes"].str.split(";", expand=True)[2]).str.split("/", expand=True)[0]
        gff["family"] = (gff["Attributes"].str.split(";", expand=True)[2]).str.split("/", expand=True)[1]
        gff["repeat_start"] = (gff["Attributes"].str.split(";", expand=True)[1]).str.split(" ", expand=True)[1]
        gff["repeat_end"] = (gff["Attributes"].str.split(";", expand=True)[1]).str.split(" ", expand=True)[2]
        data_list = []
        sine = 0
        line = 0
        ltr = 0
        dna = 0
        other = 0
        if seq_data[1].replace("\n","")[0:8] == "Scaffold":
            analysis = "uncomplete"
            for i in range(len(gff)):
                for x in range(len(seq_data)-1):
                    seq_length = len(gff["Sequence Name"][i])
                    if gff["Sequence Name"][i] + "\n" == seq_data[x][0:seq_length+1]:
                        sequence = seq_data[x].replace("\n","")[(gff["Begin"][i] + seq_length -1):(gff["End"][i] + seq_length)]
                        data_detail = {
                            "predict_id": ObjectId(inserted_id),
                            "specimen_id": specimen_id,
                            "location": gff["Sequence Name"][i],
                            "repeat": gff["Sequence Name"][i] + '_' + gff["TE"][i],
                            "start": gff["Begin"][i],
                            "end": gff["End"][i],
                            "len": gff["End"][i] - gff["Begin"][i] + 1,
                            "match": gff["match"][i].split("=")[1],
                            "class": gff["class"][i].strip("?")[6:],
                            "family": gff["family"][i] if gff["family"][i] else "-",
                            "repeat_start": int(gff["repeat_start"][i]),
                            "repeat_end": int(gff["repeat_end"][i]),
                            "repeat_size": int(gff["repeat_end"][i]) - int(gff["repeat_start"][i]) + 1,
                            "seq": sequence,
                        }
                        data_son = SON(data_detail)
                        data_list.append(data_son)
                if "LINE" in gff["class"][i][6:]:
                    line += 1
                elif "SINE" in gff["class"][i][6:]:
                    sine += 1
                elif "LTR" in gff["class"][i][6:]:
                    ltr += 1
                elif "DNA" in gff["class"][i][6:]:
                    dna += 1
                else:
                    other += 1
        elif seq_data[1].replace("\n","")[0:9] == "Chromosom" or seq_data[1].replace("\n","")[0:7] == "Plasmid":
            for i in range(len(gff)):
                for x in range(len(seq_data)):
                    if seq_data[x].replace("\n","")[0:len(gff["Sequence Name"][i])] == gff["Sequence Name"][i]:
                        sequence = seq_data[x].replace("\n", "")[(gff["Begin"][i] + (len(gff["Sequence Name"][i])-1)):(gff["End"][i] + len(gff["Sequence Name"][i]))]
                        data_detail = {
                            "predict_id": ObjectId(inserted_id),
                            "specimen_id": specimen_id,
                            "location": gff["Sequence Name"][i],
                            "repeat": gff["Sequence Name"][i] + '_' + gff["TE"][i],
                            "start": gff["Begin"][i],
                            "end": gff["End"][i],
                            "len": gff["End"][i] - gff["Begin"][i] + 1,
                            "match": gff["match"][i].split("=")[1],
                            "class": gff["class"][i][6:],
                            "family": gff["family"][i],
                            "repeat_start": int(gff["repeat_start"][i]),
                            "repeat_end": int(gff["repeat_end"][i]),
                            "repeat_size": int(gff["repeat_end"][i]) - int(gff["repeat_start"][i]) + 1,
                            "seq": sequence,
                        }
                        data_son = SON(data_detail)
                        data_list.append(data_son)
                        if gff["class"][i][6:] == "LINE":
                            line += 1
                        elif gff["class"][i][6:] == "SINE":
                            sine += 1
                        elif gff["class"][i][6:] == "LTR":
                            ltr += 1
                        elif gff["class"][i][6:] == "DNA":
                            dna += 1
                        else:
                            other += 1
        data = {
            "predict_id": ObjectId(inserted_id),
            "specimen_id": specimen_id,
            "sine": sine,
            "line": line,
            "ltr": ltr,
            "transposon": dna,
            "unclassified": other,
            "in_genome": in_genome,
            "total_base": total_base,
        }

        try:
            collection = self.db["interpersed_repeat_detail"]
            collection.insert_many(data_list)
            collection2 = self.db["interpersed_repeat_specimen"]
            collection2.insert_one(data)
            collection3 = self.db["interpersed_repeat"]
            collection3.update({'_id': inserted_id}, {'$set': {'main_id': inserted_id,"status": "end"}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (predict_gff, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % predict_gff)
        #return line,sine,ltr,dna,other,genome_total_base,total_base

    @report_check
    def add_repeat_predict_specimen(self, inserted_id, each, line_all,sine_all,ltr_all,dna_all,other_all,genome_total_base,total_base_all):
        if genome_total_base not in ['','-']:
            in_genome = float(total_base_all) / float(genome_total_base) * 100
        else:
            in_genome = '-'
        data = {
            "predict_id": ObjectId(inserted_id),
            "specimen_id": each,
            "sine": sine_all,
            "line": line_all,
            "ltr": ltr_all,
            "transposon": dna_all,
            "unclassified": other_all,
            "in_genome": round(in_genome,2),
            "total_base": total_base_all,
        }
        collection2 = self.db["interpersed_repeat_specimen"]
        collection2.insert_one(data)

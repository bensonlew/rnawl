# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId


class AnnoAntismash(Base):
    def __init__(self, bind_object):
        super(AnnoAntismash, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_anno_antismash(self, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "Antismash注释",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "version": "3.0",
            "name": "AnnoAntiSMASH_Origin",
            "software":"Antismash 5.1.2"
        }
        collection = self.db["anno_antismash"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_anno_antismash_detail(self, inserted_id, specimen_id, anno):  #not used
        data_list = []
        ann = pd.read_table(anno, sep='\t', header=0)
        ann["Location"] = ann["Cluster ID"].str.rsplit("_",n=1, expand=True)[0]
        ann["cluster"] = ann["Cluster ID"].str.rsplit("_", n=1, expand=True)[1]
        for i in range(len(ann)):
            loc = str(ann["Location"][i])
            if loc.startswith("P"):
                loc = "p" + loc.lstrip("P")
            data = {
                "antismash_id": inserted_id,
                "specimen_id": specimen_id,
                "gene_id":  ann["Gene ID"][i],
                "location": loc,
                "cluster_id": ann["cluster"][i],
                "pfam_id": ann["Pfam_id"][i],
                "gene_des": ann["NR Description"][i],

            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_antismash_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)

    @report_check
    def add_anno_antismash_stat(self, inserted_id, specimen_id, stat):
        data_list = []
        ann = pd.read_table(stat, sep='\t', header=0)
        if len(ann) < 1:
            return
        ann["Location"] = ann["Cluster ID"].str.rsplit("_", n=1, expand=True)[0]
        ann["cluster"] = ann["Cluster ID"].str.rsplit("_", n=1 , expand=True)[1]
        for i in range(len(ann)):
            data = {
                "antismash_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "location": ann["Location"][i],
                "cluster_id": ann["cluster"][i],
                "type": ann["Type"][i],
                "start": ann["Start"][i],
                "end": ann["End"][i],
                "gene_num": ann["Gene No."][i],
                "mskc": ann["Most Similar Cluster"][i],
                "similarity": ann["Similarity"][i],
                "accession": ann["MIBiG accession"][i],
                "structure": ann["predicted_structure"][i]
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_antismash_stat"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (stat, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % stat)
        ann = ann.ix[:, ["Location", "cluster", "Genes"]]
        antis = ann.drop("Genes", axis=1).join(ann["Genes"].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename("Gene ID"))
        antis.index = antis["Gene ID"]
        task_id = self.bind_object.sheet.id
        summary = self.db["anno_summary"]
        summary_coll = summary.find_one({"task_id": task_id})
        summary_id = summary_coll["_id"]
        summary_de = self.db["anno_summary_detail"]
        data_list2 = []
        for i in range(len(antis)):
            #print antis["Gene ID"][i]
            res = summary_de.find_one({"summary_id": summary_id, "specimen_id": specimen_id, "gene_id": antis["Gene ID"][i]})
            if not res:
                self.bind_object.logger.info("summary_id"+str(summary_id)+"specimen_id"+ str(specimen_id)+"gene_id" + str(antis["Gene ID"][i]))
                continue
            data = {
                "antismash_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id":  antis["Gene ID"][i],
                "location": antis["Location"][i],
                "cluster_id": antis["cluster"][i],
                "pfam_id": res["pfam"],
                "gene_des": res["gene_des"]
            }
            data_son = SON(data)
            data_list2.append(data_son)
        try:
            collection = self.db["anno_antismash_detail"]
            if len(data_list2) != 0:
                collection.insert_many(data_list2)
            else:
                self.bind_object.logger.info("导入anno_antismash_detail的数据为空")
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错" % ("anno_antismash_detail", e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % "anno_antismash_detail")


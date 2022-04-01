# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId
import re


class AnnoPfam(Base):
    def __init__(self, bind_object):
        super(AnnoPfam, self).__init__(bind_object)
        self._project_type = "fungigenome"

    @report_check
    def add_anno_pfam(self, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "Pfam注释",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": "AnnoPfam_Origin",
            "settled_params": json.dumps({"version": "pfam_v33.1"})
        }
        collection = self.db["anno_pfam"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_anno_pfam_detail(self, inserted_id, specimen_id, anno):
        data_list = []
        data_gene_list = []
        ann = pd.read_table(anno, sep='\t', header=0)
        if len(ann) < 1:
            return
        ann_gene = ann[["Gene ID", "Location"]].drop_duplicates(subset=["Gene ID", "Location"], keep='first')
        gene_li = list(set(ann_gene['Gene ID']))
        gene_li.sort(key=list(ann_gene['Gene ID']).index)   #guanqing.zou 20181112
        pfam = []
        domain = []
        nu = []
        for ge in gene_li:
            pfam_g = ";".join(ann[ann["Gene ID"] == ge]["Pfam_id"])
            domain_g = ";".join(ann[ann["Gene ID"] == ge]["Domain"])
            nu_g = len(ann[ann["Gene ID"] == ge]["Pfam_id"])
            pfam.append(pfam_g)
            domain.append(domain_g)
            nu.append(nu_g)
        ann_gene["nu"] = nu
        ann_gene["Pfam_id"] = pfam
        ann_gene["Domain"] = domain
        ann_gene.index = ann_gene["Gene ID"]
        pat = re.compile('\.[0-9]*')
        for i in range(len(gene_li)):
            data = {
                "pfam_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id":  ann_gene["Gene ID"][i],
                "location": ann_gene["Location"][i],
                "nu": ann_gene["nu"][i],
                "pfam": ann_gene["Pfam_id"][i],
                "domain": ann_gene["Domain"][i],
            }
            data['pfam'] = re.sub(pat,'',data['pfam'])
            data_son = SON(data)
            data_gene_list.append(data_son)
        pat = re.compile('\.[0-9]*')  #zouguanqing 20190103
        for i in range(len(ann)):
            pfam_g = re.sub(pat, '', ann["Pfam_id"][i])  #zouguanqing 20190103
            data = {
                "pfam_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id":  ann["Gene ID"][i],
                "location": ann["Location"][i],
                "pfam": pfam_g,
                "domain": ann["Domain"][i],
                "domain_des": ann["DomainDescription"][i],
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_pfam_detail"]
            collection.insert_many(data_list)
            collection2 = self.db["anno_pfam_stat"]
            collection2.insert_many(data_gene_list)
            main_collection = self.db["anno_pfam"]
            main_collection.update({"_id": ObjectId(inserted_id)},
                                   {"$set": {"status": 'end',"origin_id": ObjectId(inserted_id)}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)

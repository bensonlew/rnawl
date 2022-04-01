# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json,re
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId


class AnnoPfam(Base):
    def __init__(self, bind_object):
        super(AnnoPfam, self).__init__(bind_object)
        self._project_type = "bacgenome"

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
            "version" : '3.1',
            "settled_params": json.dumps({"version": "pfam_v33.1"})
        }
        collection = self.db["anno_pfam"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_anno_pfam_detail(self, inserted_id, specimen_id, anno):
        data_list = []
        data_gene_list = []
        ann = pd.read_table(anno, sep='\t', header=0)
        ll = re.compile('\.[0-9]*')
        ann['Pfam_id'] = ann['Pfam_id'].replace(ll, '')
        if len(ann) < 1:
            return
        ann_gene = ann[["Gene ID", "Location"]].drop_duplicates(subset=["Gene ID", "Location"], keep='first')
        gene_li = list(set(ann_gene['Gene ID']))
        gene_li.sort(key=list(ann_gene['Gene ID']).index)   #guanqing.zou 20181112
        pfam = []
        domain = []
        nu = []
        domain_des = []
        for ge in gene_li:
            pfam_g = ";".join(set(ann[ann["Gene ID"] == ge]["Pfam_id"]))## qingchen.zhang 20201207
            domain_g = ";".join(set(ann[ann["Gene ID"] == ge]["Domain"]))## qingchen.zhang 20201207
            domain_d_g =  ";".join(set(ann[ann["Gene ID"] == ge]["DomainDescription"]))## qingchen.zhang 20201207
            nu_g = len(set(ann[ann["Gene ID"] == ge]["Pfam_id"])) ## qingchen.zhang 20201207
            pfam.append(pfam_g)
            domain.append(domain_g)
            domain_des.append(domain_d_g)
            nu.append(nu_g)
        ann_gene["nu"] = nu
        ann_gene["Pfam_id"] = pfam
        ann_gene["Domain"] = domain
        ann_gene["Domain_desc"] = domain_des
        ann_gene['Gene ID'] = gene_li
        ann_gene.index = ann_gene["Gene ID"]
        for i in range(len(gene_li)):
            data = {
                "pfam_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id":  ann_gene["Gene ID"][i],
                "location": ann_gene["Location"][i],
                "nu": ann_gene["nu"][i],
                "pfam": ann_gene["Pfam_id"][i],
                "domain": ann_gene["Domain"][i],
                "domain_desc" : ann_gene["Domain_desc"][i]
            }
            data_son = SON(data)
            data_gene_list.append(data_son)
        for i in range(len(ann)):
            data = {
                "pfam_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id":  ann["Gene ID"][i],
                "location": ann["Location"][i],
                "pfam": ann["Pfam_id"][i],
                "domain": ann["Domain"][i],
                "domain_des": ann["DomainDescription"][i],
                "evalue": ann["DomainE-Value"][i],
                "score": ann["score"][i] # 增加score 值 zouguanqing 20190401
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_pfam_detail"]
            collection.insert_many(data_list)
            collection2 = self.db["anno_pfam_stat"]
            collection2.insert_many(data_gene_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)

# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json,re
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId


class AnnoSummary(Base):
    def __init__(self, bind_object):
        super(AnnoSummary, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_anno_summary(self, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "注释汇总",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": "AnnoSummary_Origin",
            "ref": [],
            "version" : "3.1",

        }
        collection = self.db["anno_summary"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_anno_summary_detail(self, inserted_id, specimen_id, anno, dict):
        data_list = []
        ann = pd.read_table(anno, sep='\t', header=0)
        ann = ann.fillna('-')
        ll = re.compile('\.[0-9]*')
        ann['Pfam_id'] = ann['Pfam_id'].replace(ll, '')
        if len(ann) < 1:
            return

        ###  zouguanqing
        keep_info = {
                "phi":['Gene ID', 'PHI ID','Pathogen Species','Host Species','Gene Function'],
                "vfdb":['Gene ID', 'VFDB ID', 'VFs','Description'],
                "card": ['Gene ID', 'ARO_Accession', 'ARO_description'],## qingchen.zhang 去掉ARO_category @20201109
                "tmhmm" : ['Gene ID', 'Number of predicted TMHs'],
                "tcdb": ['Gene ID', 'TCDB Description', 'TCDB Class'],
                "cazy":['Gene ID','Family','Class','Class_description'],  # 缺 family desc
                "promoter": ['Gene ID','Upstream Pos','Promoter Len'],
                "island":['Gene Id','GI No.'],
                "prephage":['Gene Id', 'Ph No.'], # 缺 'Possible Phage'
                "two_sys": ['Gene', 'Type'] , #
                "antismash":['Gene ID', 'Type', 'Most Similar Cluster','Cluster ID'],
                "signal_p":  ['Gene ID','D-score'] ,  #
                "signal_n":  ['Gene ID','D-score'] ,
            }

        rename_merge_column = {         #
            "phi":{"Gene Function":"PHI Function"},
            "vfdb": {'Description': "VFDB Description"},
            "island":{'Gene Id':'Gene ID'},
            "prephage":{'Gene Id':'Gene ID'},
            "two_sys" : {"Gene":"Gene ID", 'Type': "Two_sys_Type"},
            "antismash" : {"Type": "antismash_type", "Most Similar Cluster":"antismash_Most_similar_cluster","Cluster ID":"antismash_cluster"},
            "signal_p":{"D-score": "signal_p"},
            "signal_n":{"D-score": "signal_n"},

        }
        for k in keep_info.keys():
            if k in rename_merge_column.keys():
                for id,v in enumerate(keep_info[k]):
                    if v in rename_merge_column[k].keys():
                        keep_info[k][id] = rename_merge_column[k][v]

        mongo_key = {       ##顺序和keep_info 的信息对应,第一个元素没有用到
                "phi":['Gene ID', 'phi_id','pathogen','host_species','phi_function'],
                "vfdb":['Gene ID', 'vfdb_id', 'vfs','vfs_des'],
                "card": ['Gene ID', 'aro_accession', 'aro_desc'],## qingchen.zhang 去掉ARO_category @20201109
                "tmhmm" : ['Gene ID', 'tmh'],
                "tcdb": ['Gene ID', 'tcdb_desc','tcdb_class'],
                "cazy":['Gene ID','cazy_family','cazy_class','cazy_class_des'],  # 缺 family desc
                "promoter": ['Gene ID','pro_upstream','promoter_len'],
                "island":['Gene Id','gi_id'],
                "prephage":['Gene Id', 'ph_id'], # 缺 'Possible Phage'
                "two_sys": ['Gene', 'twosys_type'] , #
                "antismash":['Gene ID', 'anti_type', 'anti_mskc','anti_cluster'],
                "signal_p":  ['Gene ID','signal_p'] ,  #
                "signal_n":  ['Gene ID','signal_n'] ,

        }

        ### 有些分析没有结果
        current_key = []
        for k in keep_info.keys():
            if keep_info[k][1] in ann.columns:
                current_key.append(k)
        self.bind_object.logger.info(current_key)
        self.bind_object.logger.info(dict)

        for i in range(len(ann)):
            self.bind_object.logger.info(ann["Location"][i])
            data = {
                "summary_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id":  ann["Gene ID"][i],
                "location": ann["Location"][i],
                "nr_desc": ann["NR Description"][i],
                "gene_name": ann["Gene Name"][i],
                "swiss_des": ann["Swiss-Prot Description"][i],
                "cog": ann["COG ID"][i],
                "go": ann["GO ID"][i],
                "cog_type": ann["COG Type"][i],
                "ko": ann["KO ID"][i],
                "pfam": ann["Pfam_id"][i],
                "go_desc": ann["GO Description"][i],
                "cog_desc":ann["COG Description"][i],
                "domain":ann["Pfam Domain"][i],
                "domain_desc": ann["Domain Description"][i],
                "ko_desc": ann['KO Description'][i],
                "map_id": ann['Pathway'][i]
            }
            if dict:
                data['seq_type'] = dict[ann["Location"][i]]
            ##zouguanqing

            for data_name in sorted(current_key):
                for id, m_key in enumerate(mongo_key[data_name][1:],1):
                    h_key = keep_info[data_name][id]
                    data[m_key] = ann[h_key][i]


            if "promoter_len" in data.keys():
                if data['promoter_len'] != '-':
                    data['promoter'] = 'YES'
                else:
                    data['promoter'] = 'NO'


            if 'gene_ori_name' in ann.columns:
                data['gene_name'] = ann['gene_ori_name'][i]
                data['gene_des'] = ann['gene_ori_desc'][i]
            else:
                data['gene_des'] = data['nr_desc']
            ###

            data_son = SON(data)
            data_list.append(data_son)

        data_stat = {
            "summary_id": ObjectId(inserted_id),
            "specimen_id": specimen_id,
            "nr": ann['NR Description'][ann['NR Description']!='-'].count(),
            "swiss":ann['Swiss-Prot Description'][ann['Swiss-Prot Description']!='-'].count(),
            "pfam" :ann['Pfam_id'][ann['Pfam_id']!='-'].count(),
            "cog" :ann['COG ID'][ann['COG ID']!='-'].count(),
            "go" :ann['GO ID'][ann['GO ID']!='-'].count(),
            "kegg" :ann['KO Description'][ann['KO Description']!='-'].count()
        }

        try:
            collection = self.db["anno_summary_detail"]
            collection_stat = self.db["anno_summary_stat"]
            collection_stat.insert_one(data_stat)
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)

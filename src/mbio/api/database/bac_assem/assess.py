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


class Assess(Base):
    def __init__(self, bind_object):
        super(Assess, self).__init__(bind_object)
        self._project_type = "bac_assem"

    @report_check
    def add_draft_assess(self, db, params, task_id=None, project_sn=None, name=None, samp_list=[]):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': 'draft_assess',
            'created_ts': created_ts,
            'name': name if name else 'draft_assess',
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end',
            'samp_list': ",".join(samp_list)
        }
        collection = self.db['db']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': main_id},{'$set':{'main_id':main_id}})
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
    def add_assess_kmer(self, main_id, samp, table):
        main_id = self.check_id(main_id)
        data = pd.read_table(table, header=0, names=["depth", "frequency"])
        data["samp"] = samp
        data["draft_id"] = main_id
        data_list = data.to_dict(orient="records")
        self.export_db("assess_kmer", data_list, table)

    @report_check
    def add_draft_assess_pca(self, main_id, samp, table):
        main_id = self.check_id(main_id)
        # data = pd.read_table(table, header=0).reindex(columns=["Sample_ID", "PC1", "PC2", "cov"])
        data = pd.read_table(table, header=0, names=["label", "pc1", "pc2", "cov"])
        # data.columns = ["label", "pc1", "pc2"]
        data["samp"] = samp
        data["assess_id"] = main_id
        # data["cov"] = 1  # 暂时没有cov值
        data_list = data.to_dict(orient="records")
        self.export_db("draft_assess_pca", data_list, table)

    @report_check
    def add_assess_gc(self, table_name,main_id, samp, local_path, remote_path):
        main_id = self.check_id(main_id)
        data_list = []
        str_map = {1000: "1k", 3000: "3k", 5000: "5k", 8000: "8k", 10000: "10k"}
        for step in [1000, 3000, 5000, 8000, 10000]:
            if os.path.isdir(os.path.join(local_path, "depth_gc_%s" % step)):
                insert_data = {
                    "assess_id": main_id,
                    "window": str_map[step],
                    "path": os.path.join(remote_path, "depth_gc_%s" % step),
                    "samp": samp
                }
                data_list.append(insert_data)
            else:
                break
        if len(data_list) > 0:
            self.export_db(table_name, data_list, local_path)

    @report_check
    def add_assess_size(self, main_id, samp, table):
        main_id = self.check_id(main_id)
        data = pd.read_table(table, header=0, names=["method", "kmer", "k_indiv", "k_depth", "size", "r_size", "h_rate",
                                                     "repeat", "same_kmer", "accuracy_rate", "u_base", "aver_depth"])
        data["draft_id"] = main_id
        data.loc[data["h_rate"]<0.01, "h_rate"] = 0.01
        data["samp"] = samp
        data_list = data.to_dict(orient="records")
        self.export_db("assess_size", data_list, table)

    @report_check
    def add_draft_assess_16s(self, table_name, main_id, samp, table):
        main_id = self.check_id(main_id)
        data = pd.read_table(table, header=0, names=["org", "ident", "cov", "evalue", "score"])
        data["samp"] = samp
        data["organism_id"] = main_id
        if data.empty:
            return
        data_list = data.to_dict(orient="records")
        self.export_db(table_name, data_list, "16stable")

    @report_check
    def add_draft_assess_hk(self, table_name, main_id, samp, table, listname):
        main_id = self.check_id(main_id)
        data = pd.read_table(table, header=0, names=["house", "org", "ident", "cov", "evalue", "score"])
        data["samp"] = samp
        data["organism_id"] = main_id
        if data.empty:
            return
        data_list = data.to_dict(orient="records")
        self.export_db(table_name, data_list, "hk_table")
        collection = self.db['organism']
        collection.update({'_id': main_id}, {'$set': {'names': listname}})


    def get_org(self, num=10):
        taxon_batch = ["Bacillus_subtilis_BEST7613",
                       "Streptococcus_troglodytae",
                       "Arthrobacter_sp_Hiyo4",
                       "Arthrobacter_sp_Hiyo8",
                       "Streptomyces_laurentii",
                       "Petrimonas_sp_IBARAKI",
                       "Candidatus_Hodgkinia_cicadicola",
                       "Aggregatibacter_actinomycetemcomitans_D11S_1",
                       "Lactobacillus_brevis_BSO_464",
                       "synthetic_Escherichia_coli_C321deltaA",
                       "Clostridium_botulinum_CDC_297",
                       "Candidatus_Hodgkinia_cicadicola"
                       ]
        first_t = taxon_batch[random.randint(0,10)]
        taxons = [first_t] * (num -1)
        taxons.append(taxon_batch[random.randint(0,10)])
        return taxons

    def get_n_random(self, min, max, num=10):
        value = random.sample(range(0, 100), num)
        nums = pd.Series(value) / 100 * (max - min) + min
        nums_l = nums.tolist()
        return nums_l

    @report_check
    def fake_add_draft_assess_16s(self, main_id, samp):
        main_id = self.check_id(main_id)
        diction = {
            "assess_id": [main_id] * 10,
            "samp": [samp] * 10,
            "org": self.get_org(10),
            "ident": self.get_n_random(80,100),
            "cov": self.get_n_random(90,100),
            "evalue": "1e-05",
            "score": sorted(self.get_n_random(300,2000), reverse=True)
        }
        data_list = pd.DataFrame(diction).to_dict(orient="records")
        self.export_db("draft_assess_16s", data_list, "16stable")

    @report_check
    def fake_add_draft_assess_hk(self, main_id, samp):
        main_id = self.check_id(main_id)
        for house in ["dnaG", "frr", "infC", "nusA", "pgk", "pyrG", "rplA", "rplB", "rplC"]:
            diction = {
                "assess_id": [main_id] * 10,
                "samp": [samp] * 10,
                "house": [house] * 10,
                "org": self.get_org(10),
                "ident": self.get_n_random(80,100),
                "cov": self.get_n_random(90, 100),
                "evalue": "1e-05",
                "score": sorted(self.get_n_random(300,2000), reverse=True)
            }
        data_list = pd.DataFrame(diction).to_dict(orient="records")
        self.export_db("draft_assess_hk", data_list, "hk_table")
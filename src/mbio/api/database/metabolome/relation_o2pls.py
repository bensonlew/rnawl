# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo '
#  20210409

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import pandas as pd
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class RelationO2pls(ApiBase):
    def __init__(self, bind_object):
        super(RelationO2pls, self).__init__(bind_object)
        self._project_type = "metabolome"

    def add_o2plsda(self, main_id=None, project_sn='metabolome', task_id='metabolome', params=None):
        if main_id is None:
            name = "relation_o2plsda" + "_"
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Relation_o2plsda',
                params=params,
                status="start")
            main_id = self.create_db_table('relation_o2plsda', [main_info])
        else:
            main_id = ObjectId(main_id)
            try:
                self.update_db_record('relation_o2plsda', main_id)
            except Exception as e:
                self.bind_object.logger.error("导入relation_o2plsda数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入relation_o2plsda数据成功")
        return main_id

    def add_o2plsda_detail(self, main_id, file1, file2, file3, name, metab_desc, trans_desc):
        df_metab_desc = pd.read_table(metab_desc,'\t')
        df_metab_desc = df_metab_desc.fillna("-")
        dict_metab = df_metab_desc.set_index(["metab_id"])["Metabolite"].to_dict()
        if trans_desc:
            df_trans_desc = pd.read_table(trans_desc,'\t')
            df_trans_desc = df_trans_desc.fillna("-")
            dict_trans = df_trans_desc.set_index(["gene_id"])["gene_name"].to_dict()
            self.bind_object.logger.info("dict_trans:{}".format(dict_trans))
        data1 = pd.read_table(file1, sep='\t', header=0, index_col=0)
        data_list1 = []
        for i in data1.index:
            insert_data = {
                    "o2pls_id": main_id,
                    "id": i,
                    "x": data1.loc[i, "pq[1]"],
                    "y": data1.loc[i, "pq[2]"],
                    "omics": data1.loc[i, "Omics"]
                    }
            if data1.loc[i, "Omics"] == "Metabolome":
                insert_data["name"] = dict_metab[i]
            if data1.loc[i, "Omics"] == "Transcript":
                if trans_desc:
                    if i not in dict_trans.keys():
                        insert_data["name"] = "-"
                    else:
                        insert_data["name"] = dict_trans[i]
                else:
                    insert_data["name"] = "-"

            data_list1.append(insert_data)
        try:
            collection = self.db['relation_o2pls_loading']
            collection.insert(data_list1, check_keys=False)
        except Exception as e:
            self.bind_object.logger.error("导入relation_o2plsda_loading数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入relation_o2plsda_loading数据成功")

        data2 = pd.read_table(file2, sep='\t', header=0, index_col=0)
        data_list2 = []
        for j in data2.index:
            insert_data = {
                "o2pls_id": main_id,
                "name": data2.loc[j, "Sample"],
                "x": data2.loc[j, "tu[1]"],
                "y": data2.loc[j, "tu[2]"],
                "omics": data2.loc[j, "Omics"],
                "category": data2.loc[j, 'Group']
            }
            data_list2.append(insert_data)
        try:
            collection1 = self.db['relation_o2pls_scores']
            collection1.insert(data_list2, check_keys=False)
        except Exception as e:
            self.bind_object.logger.error("导入relation_o2plsda_scores数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入relation_o2plsda_scores数据成功")
            
        data3 = open(file3, "r")
        lines = data3.readlines()
        data_list3 = []
        for line in lines[1:]:
            sp = line.strip().split("\t")
            insert_data = {
                    "o2pls_id": main_id,
                    "R2X": sp[0],
                    "R2Y": sp[1],
                    "R2Xjoint": sp[2],
                    "R2Yjoint": sp[3],
                    "R2Xhat": sp[4],
                    "R2Yhat": sp[5],
                    "R2Xpred": sp[6],
                    "R2Ypred": sp[7]
                }
        data_list3.append(insert_data)
        try:
            collection2 = self.db['relation_o2pls_summary']
            collection2.insert(data_list3, check_keys=False)
        except Exception as e:
            self.bind_object.logger.error("导入relation_o2pls_summary数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入relation_o2pls_summary数据成功")        
        
        
        
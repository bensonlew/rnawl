# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
#20191104

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON

from api_base import ApiBase


class PanVenn(ApiBase):
    def __init__(self, bind_object):
        super(PanVenn, self).__init__(bind_object)
        self._project_type = "tool_lab"


    def add_venn_detail(self, inserted_id, venn_path,collection_type='pgap'):
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        data_list = []
        with open(venn_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                if (re.search(r"Group_name", line)) or (re.search(r"Sample_name", line)):
                    pass
                else:
                    line = line.strip().split("\t")
                    index_name = line[0].split(" only")[0].strip()
                    gene_num = line[1]
                    try:
                        gene_list = line[2]
                    except:
                        gene_list = "-"
                    data = {
                        "pang_id": new_inserted_id,
                        "label": index_name,
                        "cluster_list": gene_list,
                        "gene_num": gene_num,
                        }
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            if collection_type == 'pgap':
                collection = self.db["pangenome_pgap_venn"]
            else:
                collection = self.db["pangenome_hom_venn"]
            collection.insert_many(data_list)

        except Exception as e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (venn_path, e))

    def add_venn_detail_tool(self, inserted_id, venn_path,collection_type='pgap'):
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        data_list = []
        with open(venn_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                index_name = line[0]
                gene_list = ",".join(line[1].split(";"))
                data = {
                    "pang_id": new_inserted_id,
                    "name": index_name,
                    "seqs": gene_list,
                    "type": "venn"
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            if collection_type == 'pgap':
                collection = self.db["pangenome_pgap_venn"]
                main_collection = self.db["pangeome_pgap"]
            else:
                collection = self.db["pangenome_hom_venn"]
                main_collection= self.db['pangemome_hom']

            venn_data = json.dumps({"category":"name","names":"seqs","condition" : {"type":"venn"}})
            main_collection.update({"_id":new_inserted_id},{"$set":{"venn_data": venn_data}})
            collection.insert_many(data_list)


        except Exception as e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (venn_path, e))

    def add_venn2_detail_tool(self, inserted_id, venn_path,collection_type='pgap'):
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        data_list = []
        with open(venn_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                if (re.search(r"#group_name", line)) :
                    pass
                else:
                    line = line.strip().split("\t")
                    #index_name = line[0].split(" only")[0].strip()
                    index_name = line[0]

                    try:
                        gene_list = line[1]
                    except:
                        gene_list = ""
                    if gene_list:
                        gene_list = ",".join(gene_list.split(";"))
                    data = {
                        "pang_id": new_inserted_id,
                        "name": index_name,
                        "seqs": gene_list,#
                        #"gene_num": gene_num,
                        "type" : "venn",
                        "class" : "hgene"
                        }
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            if collection_type == 'pgap':
                collection = self.db["pangenome_pgap_venn"]
                main_collection = self.db["pangeome_pgap"]
            else:
                collection = self.db["pangenome_hom_venn"]
                main_collection= self.db['pangemome_hom']

            venn_data = json.dumps({"category":"name","names":"seqs","condition" : {"type":"venn"}})
            main_collection.update({"_id":new_inserted_id},{"$set":{"venn_data": venn_data}})
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (venn_path, e))
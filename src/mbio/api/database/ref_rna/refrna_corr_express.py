# -*- coding: utf-8 -*-
from pymongo import MongoClient
from bson.objectid import ObjectId
import types
from types import StringTypes
import re
import json, time
import pandas as pd
import numpy as np
import datetime, os
from bson.son import SON
from collections import Counter
import glob
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config

class RefrnaCorrExpress(Base):
    def __init__(self, bind_object):
        super(RefrnaCorrExpress, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #db = Config().MONGODB + '_ref_rna'
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        #print db
    @report_check
    def add_pca_table(self, pca_path, group_id=None,group_detail=None,name=None, params=None, express_id=None, detail=True, seq_type=None):
        # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        params = {}
        # params ['group_id'] = "58f01bbca4e1af488e52de3d"
        # group_detail = {"A":["58d8a96e719ad0adae70fa14","58d8a96e719ad0adae70fa12"], "B":["58d8a96e719ad0adae70fa11","58d8a96e719ad0adae70fa13"]}
        # params["group_detail"] = group_detail
        if group_id and group_detail:
            params['group_id']=group_id
            params['group_detail'] = group_detail
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        with open(pca_path + "/pca_importance.xls",'r+') as f1:
            pc_num = []
            pc = {}
            f1.readline()
            for lines in f1:
                line = lines.strip().split("\t")
                pc_num.append(line[0])
                pc[line[0]]=round(float(line[1]),6)
                
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "type": seq_type,
            "name": name if name else "pca_origin_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            "status": "end",
            "desc": "样本间相关性PCA分析",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "express_id": ObjectId(express_id),
            "pc_num":pc_num
            # "pc_num":["PC1","PC2","PC3","PC4"],
            # "PC1":str(0.54047),
            # "PC2":str(0.25556),
            # "PC3":str(0.20397),
            # "PC4":str(0)
        }
        if pc:
            for keys,values in pc.items():
                insert_data[keys]=values
        collection = self.db["sg_express_pca"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if detail:
            pca_file = os.path.join(pca_path, 'pca_importance.xls')
            pca_rotation = os.path.join(pca_path, 'pca_rotation.xls')
            site_file = os.path.join(pca_path, 'pca_sites.xls')
            if os.path.exists(pca_path):
                self.add_pca(pca_file=pca_file, correlation_id=inserted_id)
                #self.add_pca_rotation(input_file=pca_rotation, db_name='sg_express_pca_rotation', correlation_id=inserted_id)
                self.add_pca_rotation(input_file=site_file, db_name='sg_express_pca_rotation', correlation_id=inserted_id)
        print "pca主表导入成功！{}".format(inserted_id)
        return inserted_id

    @report_check
    def add_pca(self, pca_file, correlation_id=None):
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        data_list = []
        correlation_id = ObjectId(correlation_id)
        with open(pca_file, "r") as f:
            f.readline()
            for line in f:
                line = line.strip().split("\t")
                data = {
                    "pca_id": correlation_id,
                    line[0]: line[1]
                }
                data_list.append(data)
        try:
            collection = self.db["sg_express_pca_sites"]
            result = collection.insert_many(data_list)
        except Exception, e:
            print ("导入sg_express_correlation_pca数据出错:%s" % e)
        else:
            print ("导入sg_express_correlation_pca数据成功")

    @report_check
    def add_pca_rotation(self,input_file, db_name, correlation_id=None):
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        data_list = []
        correlation_id = ObjectId(correlation_id)
        if db_name == "sg_denovo_pca_rotation":
            col_name = "gene_id"
        else:
            col_name = "species_name"
        with open(input_file, "r") as f:
            pcas = f.readline().strip().split("\t")[1:]
            for line in f:
                line = line.strip().split("\t")
                data = {
                    "pca_id": correlation_id,
                    col_name: line[0]
                }
                for n, p in enumerate(pcas):
                    data[p] = line[n + 1]
                data_list.append(data)
        try:
            collection = self.db[db_name]
            collection.insert_many(data_list)
        except Exception, e:
            print ("导入%s数据出错:%s" % db_name, e)
        else:
            print ("导入%s数据成功" % db_name)
    
    @report_check
    def add_correlation_table(self, correlation, group_id=None,group_detail=None,name=None, params=None, express_id=None, detail=True, seq_type=None):
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        correlation_tree = glob.glob("{}/*.tre".format(correlation))
        with open(correlation_tree[0], "r") as t:
            correlation_tree = t.readline().strip()
            raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', correlation_tree)
            tree_list = [i[1] for i in raw_samp]
        # if not params:
        #     params = dict()
        #     params['express_id'] = express_id
        params = {}
        # params ['group_id'] = "589c22a27e9e39880600002a"
        # group_detail = {"a1":["589c113da4e1af32c1c0929c","589c113da4e1af32c1c0929d"], "a2":["589c113da4e1af32c1c0929e","589c113da4e1af32c1c0929f"]}
        # params["group_detail"] = group_detail
        if group_id and group_detail:
            params['group_id']=group_id
            params['group_detail'] = group_detail
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn    
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "type": seq_type,
            "name": name if name else "correlation_origin_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            "status": "end",
            "desc": "",
            "correlation_tree": correlation_tree,
            "tree_list": tree_list,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_express_correlation"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if detail:
            self.add_correlation_detail(collection=correlation, correlation_id=inserted_id, updata_tree=True)
        return inserted_id
    
    @report_check
    def add_correlation_detail(self, collection, correlation_id=None, updata_tree=False):
        # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        correlation_id = ObjectId(correlation_id)
        correlation_matrix = collection + "/correlation_matrix.xls"
        data_list = []
        with open(correlation_matrix, "r") as m:
            samples = m.readline().strip().split()
            for line in m:
                data = {
                    "correlation_id": correlation_id
                }
                line = line.strip().split()
                for i, s in enumerate(samples):
                    data["specimen_name"] = line[0]
                    data[s] = line[i+1]
                # bind_object.logger.error data
                data_list.append(data)
        if updata_tree:
            col_tree = collection + "/corr_col.tre"
            row_tree = collection + "/corr_row.tre"
            f = open(col_tree, "r")
            r = open(row_tree, "r")
            tree_col = f.readline().strip()
            tree_row = r.readline().strip()
            raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', tree_col)
            tree_list = [i[1] for i in raw_samp]
            f.close()
            r.close()
            collection_first = self.db['sg_express_correlation']
            collection_first.update({"_id": ObjectId(correlation_id)}, {"$set": {"correlation_tree": tree_row, "row_tree": tree_row, "col_tree": tree_col, "specimen": tree_list}})
        try:
            collection = self.db["sg_express_correlation_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            print ("导入相关系数分析数据出错:%s" % e)
        else:
            print ("导入相关系数分析数据成功")
    
    @report_check
    def add_venn(self, group_id,group_detail, express_id, venn_table=None, venn_graph_path=None, params=None, query_type=None,analysis_name="express"):
        #db = Config().MONGODB + '_ref_rna'\
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        if not isinstance(express_id, ObjectId):
            express_id = ObjectId(express_id)
        params={}
        params['group_id'] = group_id
        params['group_detail'] = group_detail
        # params ['group_id'] = "58f01bbca4e1af488e52de3d"
        # group_detail = {"A":["58d8a96e719ad0adae70fa14","58d8a96e719ad0adae70fa12"], "B":["58d8a96e719ad0adae70fa11","58d8a96e719ad0adae70fa13"]}
        # params["group_detail"] = group_detail
        params["type"] = query_type
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "express_id": str(express_id),
            "name": "venn_table_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
            "status": "end",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "params": (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
        }
        if analysis_name=="express":
            collection = self.db['sg_express_venn']
        if analysis_name =='geneset':
            collection = self.db["sg_geneset_venn"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if venn_table:
            self.add_venn_detail(venn_table, inserted_id,project='ref',analysis_name=analysis_name)
        if venn_graph_path:
            self.add_venn_graph(venn_graph_path=venn_graph_path, venn_id=inserted_id, project='ref',analysis_name=analysis_name)
        print inserted_id
        return inserted_id

    def add_venn_detail(self, venn_table, venn_id, project = None,analysis_name="express"):
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        data_list = []
        if not isinstance(venn_id, ObjectId):
            if isinstance(venn_id, StringTypes):
                venn_id = ObjectId(venn_id)
            else:
                raise Exception("venn_id必须为ObjectId对象或其对应的字符串!")
        with open(venn_table, "r") as f:
            for lines in f:
                line = lines.strip('\n').split("\t")
                insert_data = {
                    'venn_id': venn_id,
                    'label': line[0],
                    'num': line[1],
                    'gene_ids': line[2]
                }
                data_list.append(insert_data)
        try:
            if project == 'denovo':
                collection = self.db["sg_denovo_venn_detail"]
            if project == 'ref':
                if analysis_name =="express":
                    collection = self.db['sg_express_venn_detail']
                if analysis_name == "geneset":
                    collection = self.db['sg_geneset_venn_detail']
            collection.insert_many(data_list)
            print 'haha1'
        except Exception, e:
            print ("导入Venn数据出错:%s" % e)
        else:
            print ("导入Venn数据成功")

    def add_venn_graph(self,venn_graph_path, venn_id, project='meta',analysis_name="express"):
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        data_list = []
        if not isinstance(venn_id, ObjectId):
            if isinstance(venn_id, StringTypes):
                venn_id = ObjectId(venn_id)
            else:
                raise Exception("venn_id必须为ObjectId对象或其对应的字符串!")
        with open(venn_graph_path, "r") as f:
            f.readline()
            for lines in f:
                line = lines.strip().split("\t")
                insert_data = {
                    'venn_id': venn_id,
                    # 'otu_id': ObjectId(otu_id),
                    'category_name': line[0]
                }
                if project == 'meta':
                    insert_data['otu_names'] = line[1]
                if project == 'denovo' or project == 'ref':
                    insert_data['gene_list'] = line[1]
                data_list.append(insert_data)
        try:
            if project == 'meta':
                collection = self.db["sg_otu_venn_graph"]
            if project == 'denovo':
                collection = self.db["sg_denovo_venn_graph"]
            if project == 'ref':
                if analysis_name == "express":
                    collection = self.db['sg_express_venn_graph']
                if analysis_name == 'geneset':
                    collection = self.db['sg_geneset_venn_graph']
            collection.insert_many(data_list)
        except Exception, e:
            print ("导入Venn画图数据出错:%s" % e)
        else:
            print ("导入Venn画图数据成功")
        
if __name__ == "__main__":
    db = MongoClient("192.168.10.189:27017").tsanger_ref_rna
    venn_path = "/mnt/ilustre/users/sanger-dev/workspace/20170401/Single_sample_venn/VennTable/output"
    #add_denovo_venn(express_id="58dc5688a4e1af190ea558d6", venn_table=venn_path+"/venn_table.xls", venn_graph_path=venn_path+"/venn_graph.xls", params=None)
    d = RefrnaCorrExpress()
    corr_path = "/mnt/ilustre/users/sanger-dev/workspace/20170425/ExpressCorr_express_corr_4/Correlation/output"
    d.add_correlation_table(corr_path,express_id="58ef0bcba4e1af740ec4c14c", detail=False, seq_type="gene")
    print 'end'
    
    # correaltion_path = "/mnt/ilustre/users/sanger-dev/workspace/20170413/Single_rsem_stringtie_zebra_9/Express/output/correlation/genes_correlation"
    # add_correlation_table( correlation=correaltion_path, name=None, params=None, express_id=ObjectId("58ef0bcba4e1af740ec4c14c"), detail=True, seq_type="gene")
    # add_correlation_detail( collection=correaltion_path, correlation_id="58df1c26a4e1af0c7809a64e", updata_tree=True)
    # pca_path = "/mnt/ilustre/users/sanger-dev/workspace/20170315/Single_small_rsem_stringtie_5/Express/output/pca/genes_pca"
    # pca_path = "/mnt/ilustre/users/sanger-dev/workspace/20170413/Single_rsem_stringtie_zebra_9/Express/output/pca/genes_pca"
    # add_pca_table( pca_path=pca_path, name=None, params=None, express_id="58ef0bcba4e1af740ec4c14c", detail=True, seq_type="gene")
    # add_pca(pca_file=pca_path + "/pca_importance.xls", correlation_id="58df1f04a4e1af062f91b438")
    # add_pca_rotation( input_file=pca_path + "/pca_rotation.xls", db_name="sg_express_pca_rotation", correlation_id="58df1f04a4e1af062f91b438")
    # add_pca_rotation( input_file=pca_path + "/pca_sites.xls", db_name="sg_express_pca_sites", correlation_id="58df1f04a4e1af062f91b438")

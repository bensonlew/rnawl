# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
#20190902

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
import pandas as pd
import numpy as np

from api_base import ApiBase


class Pan(ApiBase):
    def __init__(self, bind_object):
        super(Pan, self).__init__(bind_object)
        self._project_type = "tool_lab"


    def add_pan_detail(self, inserted_id, cluster_data,collection_type='pgap'):
        """
        导表的cluster表
        :param inserted_id: 主表id
        :param cluster: cluster表
        :return:
        """
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        self.bind_object.logger.info("cluster_data:path %s" % (cluster_data))
        with open(cluster_data, "r") as f:
            data_list = []
            header = f.readline().strip().split("\t")
            sample_list = header[3:]
            sample_num = len(sample_list)

            for line in f:
                line = line.strip().split("\t")
                data = {
                    "cluster_id": line[0],
                    "sample_n": line[1]+'/'+str(sample_num),
                    "gene_n": str(line[2]),
                    "pang_id": new_inserted_id,
                }
                if int(line[1]) == len(sample_list) and int(line[1]) == int(line[2]):
                    one_copy = 1
                else:
                    one_copy = 0

                data['one_copy'] = one_copy
                for i in range(3,len(sample_list)+3):
                    sample = sample_list[i-3].replace('.', '_')
                    if line[i] != '-':
                        gene_list = line[i].split(",")
                        data[sample] = str(len(gene_list))
                    else:
                        data[sample] = '0'

                data_son = SON(data)
                data_list.append(data_son)
        try:
            if collection_type == 'pgap':
                collection = self.db["pangenome_pgap_detail"]
            else:
                collection = self.db["pangenome_hom_detail"]

            collection.insert_many(data_list)

        except Exception as e:
            self.bind_object.logger.info("导入cluster结果表出错:%s" % (e))

    def add_pan_genome(self, inserted_id, type, genomes, genome_data,collection_type='pgap'):
        """
        导表计算pan、core、newgenes的数据和公式
        """
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        data_list = []
        data = pd.read_table(genome_data, sep="\t", header=0)###画图的10条数据
        data1 = data.T
        sample_list = list(data.columns)
        index_num = len(data1.index)
        max_x  = 0
        max_y = 0

        for i in range(index_num):

            index_name = str(sample_list[i]).strip("g")

            new_list = list(data1.iloc[i])

            # num_list = map(lambda x:str(i), new_list)
            #num_list =  [str(x) for x in new_list]
            for n in new_list:
                data2 = {
                    "pang_id": new_inserted_id,
                    "x" : int(index_name),
                    "y" : float(n),
                    "sub_type" : type
                }
                if float(index_name) > max_x:
                    max_x = float(index_name)

                if float(n) > max_y:
                    max_y = float(n)

                data2["category"] = 'all'
                if type in ["newgene"]:
                    data2["type"] = "column"
                    data_son = SON(data2)
                    data_list.append(data_son)
                    break

                else:
                    data2["type"] = "scatter"

                    data_son = SON(data2)
                    data_list.append(data_son)


        coresize = "-"
        pansize = "-"
        if type in ["pan"]:
            pan_number = list(data1.iloc[index_num - 1])
            pan_number = [float(x) for x in pan_number]
            pansize = float(pan_number[0])

        elif type in ["core"]:
            core_number = list(data1.iloc[index_num - 1])
            core_number = [float(x) for x in core_number]
            coresize = float(core_number[0])


        with open(genomes, 'r') as f: ###画图的公式和划图线的一条数据
            data = {}
            lines = f.readlines()
            if type in ["newgene"]:
                if len(lines) >= 4:
                    formula = lines[2].strip()
                    formula = re.sub("\"", "", formula)
                    lin = lines[5].strip().split("\t")
                    # self.bind_object.logger.info("cluster_data:path %s" % (lin))
                    for i in range(len(sample_list)):
                        sample = sample_list[i].replace('g', '')
                        data['x'] = int(sample)
                        data['y'] = float(lin[i])

                        data["pang_id"] = new_inserted_id
                        data['sub_type'] = type
                        data["type"] = "line"
                        data["name"] = ''
                        data["category"] = 'all'
                        data_son = SON(data)
                        data_list.append(data_son)
            elif type in ["pan", "core"]:
                if len(lines) >= 4:
                    if not lines[2].startswith("#"):
                        formula = lines[2].strip()
                        formula = re.sub("\"", "", formula)
                    else:
                        formula = lines[3].strip()
                        formula = re.sub("\"", "", formula)
                    for li in lines[4:]:
                        if li.startswith("#"):
                            continue
                        if li.startswith("#") or li.startswith("~"):
                            pass
                        else:
                            new_lin = li.strip().split("\t")
                            geno_name = int(new_lin[0])
                            gene_num = float(new_lin[1])
                            # data[geno_name] = str(float(gene_num))
                            data['x'] = geno_name
                            data['y'] = gene_num
                            data["pang_id"] = new_inserted_id
                            data['sub_type'] = type
                            data["type"] = "line"
                            data["name"] = ''
                            data["category"] = 'all'

                            data_son = SON(data)
                            data_list.append(data_son)
            try:
                if collection_type == 'pgap':
                    #collection = self.db["pangenome_pgap_graph"]
                    collection_stat = self.db["pangenome_pgap_stat"]
                    main_collection = self.db["pangenome_pgap"]
                else:
                    #collection = self.db["pangenome_hom_graph"]
                    collection_stat = self.db["pangenome_hom_stat"]
                    main_collection = self.db["pangenome_hom"]

                #collection.insert_many(data_list)


                #if len(lines) >= 4:
                if type in ["pan"]:
                    self.bind_object.logger.info("开始更新主表的公式")
                    if collection_type == 'pgap':
                        if len(lines) >= 4:
                            #collection.insert_one({"text": formula, 'pang_id':new_inserted_id, "type":"text","name":"hom","sub_type":"pan","x":float(max_x)/2,"y":float(max_y),"location":"up"})
                            #collection.insert_one({"text": formula, 'pang_id':new_inserted_id, "type":"text","name":"hom","sub_type":"pan","x":float(max_x)/2,"y":0,"location":"down"})
                            main_collection.update({"_id": new_inserted_id}, {'$set': {'main_id':new_inserted_id}})
                        else:
                            main_collection.update({"_id": new_inserted_id}, {'$set': {'main_id':new_inserted_id}})
                    else:
                        main_collection.update({"_id": new_inserted_id}, {'$set': {'main_id':new_inserted_id}})


                    collection_stat.insert_one({
                        "pang_id" : new_inserted_id,
                        "pan_size" : str(pansize),
                        "core_size" : '',
                    })

                elif type in ["core"]:
                    #if collection_type == 'pgap':
                    #    self.bind_object.logger.info("开始更新主表的公式")
                        #if len(lines) >= 4:
                            #collection.insert_one({"text": formula, 'pang_id':new_inserted_id, "type":"text","name":"hom","sub_type":"core","x":float(max_x)/2,"y":float(max_y),"location":"up"})
                            #collection.insert_one({"text": formula, 'pang_id':new_inserted_id, "type":"text","name":"hom","sub_type":"core","x":float(max_x)/2,"y":0,"location":"down"})

                    collection_stat.update({"pang_id": new_inserted_id}, {'$set': {"core_size": coresize}})

                # elif type in ["newgene"]:
                #     if collection_type == 'pgap':
                #         self.bind_object.logger.info("开始更新主表的公式")
                        #if len(lines) >= 4:
                            #collection.insert_one({"text": formula, 'pang_id':new_inserted_id, "type":"text","name":"newgene_formula","sub_type":"newgene"})


            except Exception as e:
                self.bind_object.logger.info("导入pangenomes结果表出错:%s" % (e))
            else:
                self.bind_object.logger.info("导入pangenomes结果表成功!")

    def update_main(self,inserted_id,cluster_data, collection_type='pgap', png_info=None):
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id

        if collection_type == 'pgap':
            main_collection = self.db["pangenome_pgap"]
        else:
            main_collection = self.db["pangenome_hom"]


        stat_table = [
                    {"filter": "false", "field": "pan_size", "type": "string", "sort": "false", "title": "Pangenome Size"},
                    {"filter": "false", "field": "core_size", "type": "string", "sort": "false", "title": "Coregenome Size"}
                ]

        detail_table = [
                    {"filter": "false", "field": "cluster_id", "type": "string", "sort": "false", "title": "Cluster ID"},
                    {"filter": "false", "field": "sample_n", "type": "string", "sort": "false", "title": "Sample No."},
                    {"filter": "false", "field": "gene_n", "type": "string", "sort": "false", "title": "Gene No."}
                ]

        tmp_sample_list = []
        with open(cluster_data, "r") as f:
            header = f.readline().strip().split("\t")
            sample_list = header[3:]
            for s in sample_list:
                sample = s.replace('.', '_')
                tmp_sample_list.append([sample, s])

        for new_s,old_s in tmp_sample_list:
            detail_table.append({"filter": "false", "field": new_s, "type": "string", "sort": "false", "title": old_s})


        #scatter_data =  json.dumps({"name":"name","data":["x","y"],"category":"category","condition":{"type":"scatter"}})
        #line_data = json.dumps({"name":"name","x":"x","y":"y","condition":{"type":"line"}})
        #column_data = json.dumps({"name":"x","data":"y","condition":{"type":"column"}})
        venn_data = json.dumps({"category":"name","names":"seqs","condition" : {"type":"venn"}})
        #text_data = json.dumps({"name":"name","condition" : {"type":"text"}})
        update_data = {"stat_table":json.dumps({"column":stat_table,"condition":{}}),
                        # "line_data" :line_data,
                        # "scatter_data":scatter_data,
                        # "column_data": column_data,
                         "venn_data" :venn_data,
                        # "text_data" :text_data,
                         "detail_data" :json.dumps({"column":detail_table,"condition":{}}),
                         #"pan_formula" :  json.dumps([{"field":"hom", "title":"公式1"}]),
                         #"core_formula" : json.dumps([{"field":"hom", "title":"公式1"}])
                    }
        if png_info:
            #update_data["remote_pic_path"] = png_info
            if 'core' in png_info:
                update_data['core_path'] = png_info['core']
            if 'pan' in png_info:
                update_data['pan_path'] = png_info['pan']
            if 'new' in png_info:
                update_data['new_path'] = png_info['new']

        main_collection.update({"_id":new_inserted_id},{"$set": update_data})





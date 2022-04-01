# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190313

import os
import re
import json
from biocluster.config import Config
from bson.objectid import ObjectId

project_type = "basic_molecular"
client = Config().get_mongo_client(mtype=project_type)
db = client[Config().get_mongo_dbname(project_type)]

def export_decarrier_params(data, option_name, dir_path, bind_obj=None):
    """
    导出去载体的samplelist
    """
    decarrier_id = ObjectId(data)
    sample_txt = os.path.join(dir_path,"sample.txt")
    result = db["sg_decarrier"].find_one({"_id": decarrier_id})
    task_info = result["task_sn"]
    sample_list = json.loads(result["sample_list"])
    for sn_info in sample_list:
        if sn_info["seq_path"] == "" or sn_info["ab1_path"] == "":
            continue
        # self.sn, anlysis_type, ab1_s3_path,seq_s3_path,self.info = fd[0],fd[1],fd[2],fd[3],fd[4]
        with open(sample_txt,"a") as st:
            st.write("{}\t{}\t{}\t{}\t{}\n".format(
                sn_info["sample_name"],sn_info["anlysis_type"],sn_info["seq_path"],sn_info["ab1_path"],task_info))
    return sample_txt

def export_lmplot_params(data, option_name, dir_path, bind_obj=None):
    """
    导出标准曲线的lmplot_dir
    """
    os.system('cd {} && mkdir lmplot_dir'.format(dir_path))
    lmplot_id = ObjectId(data)
    # lm_txt = os.path.join(dir_path,"lmplot.txt")
    result = db["sg_lmplot"].find_one({"_id": lmplot_id})
    json_format = json.loads(result["format_data"])   
    name_txt = os.path.join(dir_path,"name.txt")
    gene_txt = os.path.join(dir_path,"gene.txt")
    lmplot_dir = os.path.join(dir_path+"/lmplot_dir")
    get_lmplot_dir = [] 
    for lmplot_info in json_format:
        with open(name_txt,"w+") as r:
            r.write("{}\n".format(lmplot_info["id"]))
        with open(name_txt, 'r') as n:
            name = n.readlines()
            name_info = name[0].strip().split('\t')[0]
        for gene_info in lmplot_info["format_data"]:
            with open(name_txt,"w+") as r1:
                r1.write("{}\n".format(gene_info["gene"]))
            with open(name_txt, 'r') as n1:
                gene = n1.readlines()
                gene_info = gene[0].strip().split('\t')[0]
        lm_txt = os.path.join(dir_path+"/lmplot_dir", "{}_{}.txt".format(gene_info,name_info))
        get_line_info = []
        print lmplot_info
        for line_info in lmplot_info["format_data"]:
            with open(lm_txt,"a") as st:
                st.write("{}\t{}\n".format(
                    line_info["ct"],line_info["co_number"]))
            get_line_info.append(line_info)
        with open(lm_txt, 'r+') as st1:
            content = st1.read()
            st1.seek(0,0)
            st1.write('CT\tQ\n'+content)
        get_lmplot_dir.append(lmplot_info)
    return (lmplot_dir)

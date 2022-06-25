# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# 2021年2月1日15:09:15

import glob
import os
import re
import shutil

import numpy as np
import pandas as pd
from biocluster.config import Config
from bson.objectid import ObjectId
from mbio.packages.medical_transcriptome.copy_file import CopyFile

class GeneInfoSupple(object):
    def __init__(self,project_type="medical_transcriptome",task_id = None,id2namedict = None ,add_columns =[],level = None,file_path =None,split = ";",annotype = "go"):
        self.level = level
        # self.database = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        if id2namedict:
            self.id2name_dict = id2namedict
        else:
            if not task_id:
                raise ValueError("task_id 不能为空值")
            else:
                self.task_id = task_id
                self.id2name_dict = self.get_id2name_by_task_id(task_id)
        self.add_columns =add_columns
        self.file_path = file_path
        self.split = split


    def get_id2name_by_task_id(self,task_id):
        annot_table = self.database['sg_exp']
        annot_main = annot_table.find_one({"task_id": task_id, "level": "T", 'status': 'end'})
        annot_main_id = annot_main['main_id']
        annot_detail = self.database['sg_exp_detail']
        query_dict = dict(exp_id=annot_main_id, )
        result_dict = dict(_id=0, gene_name=1, gene_id=1, transcript_id=1)
        result = annot_detail.find(query_dict, result_dict)
        gene_annot = pd.DataFrame(list(result))
        if self.level.lower() == 't':
            id2gene_name = dict(zip(gene_annot['transcript_id'], [x if x else '-' for x in gene_annot['gene_name']]))
        else:
            id2gene_name = dict(zip(gene_annot['gene_id'], [x if x else '-' for x in gene_annot['gene_name']]))
        return id2gene_name


    def add_gene_name(self,file_path,split,annot_type,add_columns,outfile_path):
        def convert_id2name(series,col_name,split,annot_type):
            if annot_type.lower() in  ["go","reactome","do","kegg_enrich"]:
                ids = series[col_name].split(split)
                names = [self.id2name_dict[x] for x in ids if x!= ""]
                name = split.join(names)
                return name
            if annot_type.lower() == "kegg":
                ids = series[col_name].split(split)
                gene_ids = [x.split("(")[0] for x in ids if x != ""]
                names = [self.id2name_dict[x] for x in gene_ids if x in self.id2name_dict]
                name = split.join(names)
                return name

        raw_df = pd.read_table(file_path, sep="\t", index_col=False)
        for column in add_columns:
            if not column in raw_df.columns :
                raise ValueError("{} 必须为待改列列名".format(column))
            new_column = column + "_name"
            raw_df[column].fillna("",inplace=True)
            raw_df[new_column] = raw_df.apply(convert_id2name,args=(column,split,annot_type),axis=1)
        raw_df.to_csv(outfile_path,sep ="\t",index = False)
















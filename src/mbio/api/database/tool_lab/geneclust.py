# -*- coding: utf-8 -*-



import os
from bson.objectid import ObjectId
import pandas as pd
import glob
from api_base import ApiBase
import json
import numpy as np
from collections import OrderedDict

class Geneclust(ApiBase):
    def __init__(self, bind_object):
        super(Geneclust, self).__init__(bind_object)

    def add_detail(self,main_id,  result_path):
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        data_list = []

        data = pd.read_table(result_path, sep='\t')

        fix_cols =['Genecluster','Gene id','Gene name','Location','Start','End','Strand','Description']
        others = set(data.columns.tolist()) - set(fix_cols)

        field_map = {
            'Genecluster': 'clust',
            'Gene id' : 'gene_id',
            'Gene name': 'gene_name',
            'Location':'location',
            'Start' : 'start',
            'End' : 'end',
            'Strand' : 'strand',
            'Description' : 'desc'
        }

        label_num = 0
        for other in others:
            label_num +=1
            field_map[other] = 'label' +str(label_num)

        data.fillna("",inplace=True)
        for index in data.index:
            insert_data = {"clust_id":main_id}
            for col in data.columns:
                k = field_map[col]
                v = data[col][index]
                insert_data[k] = v
            data_list.append(insert_data)

        try:
            collection = self.db["geneclust_detail"]
            collection.insert_many(data_list)
            #if 'Gene id' in data.columns:
            arrows_axis_data = OrderedDict([("name", "gene_name"),("id", "gene_id"), ("start", "start"), ("end" ,"end"),("group", "label1"),("y","clust"),("condition",{})])
            #else:
            #    arrows_axis_data = OrderedDict([("name", "gene_name"),("id", "gene_name"), ("start", "start"), ("end" ,"end"),("group", "label1"),("y","clust"),("condition",{})])
            arrows_axis_data = json.dumps(arrows_axis_data)
            groups = json.dumps([{"title": k, "field": field_map[k] } for k in others])
            names = json.dumps([{"title": "Gene name", "field":"gene_name"},{"title": "Gene id", "field":"gene_id"}])
            clusts = json.dumps(list(set(data['Genecluster'].tolist())))
            self.db["geneclust"].update({"_id": main_id}, {"$set": {"arrows_axis_data":arrows_axis_data,"groups":groups,'names':names,'clusts':clusts, 'main_id': main_id}})
        except Exception as e:
            self.bind_object.logger.error("导入geneclust_detail数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入geneclust_detail数据成功")


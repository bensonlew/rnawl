# !/usr/bin/python
# -*- coding: utf-8 -*-
from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check
import json
import pandas as pd
import datetime
import glob
import re


class Procrustes(Base):
    def __init__(self, bind_object):
        super(Procrustes, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_procrustes_detail(self, summaryfile, main_id):
        data_list=[]
        with open(summaryfile, 'rb') as r:
            for line in r:
                if re.match("^#",line):
                    continue
                line = line.strip().split("\t")
                insert_data = {
                    'procrus_id': ObjectId(main_id),
                    'dim': int(line[2]),
                    'num': int(line[4]),
                    'm2': float(line[5]),
                    'p': float(line[3])
                }
                data_list.append(insert_data)
        try:
            collection = self.db["metabset_procrustes_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入procrustes数据出错: %s" % e)
            self.bind_object.set_error("Error: failed to import procrusts_detail data: %s" % e)
        else:
            self.bind_object.logger.info("导入procrustes detail数据成功")


    @report_check
    def add_procrustes_graph(self, ref, query, main_id):
        importance_list = []
        site_list = []
        site_ref = {}
        site_query = {}
        importance = 0
        site = 0
        '''
        读ref结果文件：
        获取解释度数据并写入mongo库；获取ref的样本坐标数据
        '''
        with open(ref, 'rb') as r:
            for line in r:
                arr = line.strip().split("\t")
                if not len(line) or line=="\n" or line=="\r\n":
                    importance = 0
                    site = 0
                    continue
                elif re.match("Proportion", line):
                    importance = 1
                    site = 0
                    continue
                elif importance==1:
                    for i in range(0, len(arr)):
                        name = "PC" + str(i+1)
                        #self.bind_object.logger.info("importance_name: %s" % name)
                        importance_data = {
                            'procrus_id': ObjectId(main_id),
                            'name': name,
                            'proportion_variance': float(arr[i]),
                            'type': 'importance'
                        }
                        importance_list.append(importance_data)
                    importance = 0
                elif re.match("Site", line):
                    site = 1
                    continue
                elif site==1:
                    site_tmp=[]
                    for i in range(1, len(arr)):
                        site_tmp.append(arr[i])
                    site_ref[arr[0]]=site_tmp
        '''
        读query结果文件：
        读取query的样本坐标数据
        '''
        with open(query, 'rb') as r:
            for line in r:
                arr = line.strip().split("\t")
                if not len(line) or line=="\n" or line=="\r\n":
                    site = 0
                    continue
                elif re.match("Site", line):
                    site = 1
                    continue
                elif site==1:
                    site_tmp1=[]
                    for i in range(1, len(arr)):
                        site_tmp1.append(arr[i])
                    site_query[arr[0]]=site_tmp1
        '''
        遍历ref和query的样本坐标数据，写入mongo库
        '''
        for sample in site_ref.keys():
            site_data={
                'procrus_id': ObjectId(main_id),
                'specimen': sample,
                'type': 'specimen',
            }
            if len(site_ref[sample]) != len(site_query[sample]):
                self.bind_object.logger.info("错误：ref和query的site坐标维度必须相同")
                self.bind_object.set_error("ERROR: site length of ref and query must be equal!")
            for j in range(0, len(site_query[sample])):
                site_data['p' + str(j+1)] = site_ref[sample][j]
                site_data['c' + str(j+1)] = site_query[sample][j]

            site_list.append(site_data)

        try:
            collection = self.db["metabset_procrustes_graph"]
            collection.insert_many(importance_list)
            collection.insert_many(site_list)
        except Exception, e:
            self.bind_object.logger.info("导入procrustes作图样本数据出错" % e)
            self.bind_object.set_error("Error: failed to import procrusts_graph data: %s" % e)
        else:
            self.bind_object.logger.info("导入procrustes作图样本数据成功")










if __name__ == '__main__':
    Procrustes = Procrustes(None)
    Procrustes.add_procrustes_detail("/mnt/ilustre/users/sanger-dev/sg-users/liulinmeng/metabolome/api/procrustes/outputwf/procrustes_summary.txt", "5bc43a21a4e1af13c6100574")
    Procrustes.add_procrustes_graph("/mnt/ilustre/users/sanger-dev/sg-users/liulinmeng/metabolome/api/procrustes/outputwf/metab.transformed_reference.txt","/mnt/ilustre/users/sanger-dev/sg-users/liulinmeng/metabolome/api/procrustes/outputwf/asso.transformed_query.txt","5bc43a21a4e1af13c6100574")

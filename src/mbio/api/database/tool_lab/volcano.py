# -*- coding: utf-8 -*-
# __author__ = 'zzg'
from bson.objectid import ObjectId
import datetime
import json
import re
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase
import math


class Volcano(ApiBase):
    def __init__(self, bind_object):
        super(Volcano, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_main(self, file, x_method,y_data,y_method,vip_method,p_method1,p_method2,p_method3,vip_method1,vip_method2,vip_method3,up_method1,up_method2,diff_gene,diff_gene_num,diff_gene_name,vip=None,project_sn='tool_lab', main_id=None, task_id='tool_lab', params=None):
        if main_id is None:
            name = "volcano"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='volcano',
                params= params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('volcano', [main_info])
        else:
            main_id = ObjectId(main_id)

        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        insert_data = []
        scatter_group = []

        diff_gene_name_list = []
        if diff_gene == "true":
            data = []
            with open(file, "r") as v:
                data1 = v.readlines()
                for x in data1[1:]:
                    data.append(abs(float(x.strip().split("\t")[2])))
            data.sort(reverse=False)
            diff_num = data[int(diff_gene_num)-1]
            if diff_gene_name:
                diff_gene_name_list = diff_gene_name.strip().split("\n")
        else:
            diff_num = 0

        with open(file,"r") as f:
            data = f.readlines()
            for i in data[1:]:
                num = 0
                if float(i.strip().split("\t")[2]) <= diff_num or i.strip().split("\t")[0] in diff_gene_name_list:
                    show = "true"
                else:
                    show = "false"
                if x_method == "log2":
                    x_data = math.log(float(i.strip().split("\t")[1]),2)
                elif x_method == "log10":
                    x_data = math.log(float(i.strip().split("\t")[1]), 10)
                else:
                    x_data = float(i.strip().split("\t")[1])
                if y_data == "P_value":
                    y_data1 = float(i.strip().split("\t")[2])
                elif y_data == "FDR":
                    y_data1 = float(i.strip().split("\t")[3])
                elif y_data == "VIP_pred_OPLS-DA":
                    y_data1 = float(i.strip().split("\t")[4])
                elif y_data == "VIP_PLS-DA":
                    y_data1 = float(i.strip().split("\t")[5])
                if y_method == "log10":
                    y_data2 = math.log(float(y_data1),0.1)
                elif y_method == "log2":
                    y_data2 = math.log(float(y_data1),0.5)
                else:
                    y_data2 = float(y_data1)
                if p_method1 == "P":
                    p_data = float(i.strip().split("\t")[2])
                else:
                    p_data = float(i.strip().split("\t")[3])
                if p_method2 == "less":
                    if p_data < float(p_method3):
                        num += 1
                else:
                    if p_data <= float(p_method3):
                        num += 1
                if vip_method1 == "vip_oplsda":
                    vip_data = float(i.strip().split("\t")[4])
                else:
                    vip_data = float(i.strip().split("\t")[5])
                if vip_method2 == "more":
                    if vip_data > float(vip_method3):
                        num += 1
                else:
                    if vip_data >= float(vip_method3):
                        num += 1
                if up_method1 == "more":
                    if num == 2:
                        if float(i.strip().split("\t")[1]) > float(up_method2):
                            category = "up"
                            scatter_group.append("up")
                        if float(i.strip().split("\t")[1]) <  1 / float(up_method2):
                            category = "down"
                            scatter_group.append("down")
                    else:
                        category = "nosig"
                        scatter_group.append("nosig")
                else:
                    if num == 2:
                        if float(i.strip().split("\t")[1]) >= float(up_method2):
                            category = "up"
                            scatter_group.append("up")
                        if float(i.strip().split("\t")[1]) <= 1 / float(up_method2):
                            category = "down"
                            scatter_group.append("down")
                    else:
                        category = "nosig"
                        scatter_group.append("nosig")

                if vip_method == "true":
                    if vip == "VIP_pred_OPLS-DA":
                        vip_data = float(i.strip().split("\t")[4])
                        size = float(i.strip().split("\t")[4])
                    else:
                        vip_data = float(i.strip().split("\t")[5])
                        size = float(i.strip().split("\t")[5])
                else:
                    size = 1.0
                if diff_gene == "true":
                    insert_data.append({
                        "volcano_id": main_id,
                        "metabolite": re.sub(r",|\"|'|/", "", i.split("\t")[0]),
                        "fc": round(x_data, 3),
                        "pvalue": round(y_data2, 3),
                        "category": category,
                        "size": size,
                        "vip_plsda": round(float(vip_data), 3),
                        "type": "scatter",
                        "show_name": show
                    })
                else:
                    insert_data.append({
                        "volcano_id": main_id,
                        "metabolite": re.sub(r",|\"|'|/", "", i.split("\t")[0]),
                        "fc": round(x_data, 3),
                        "pvalue": round(y_data2, 3),
                        "category": category,
                        "size": size,
                        "vip_plsda": round(float(vip_data), 3),
                        "type": "scatter"
                    })
            name1 = x_method +"(" + data[0].split("\t")[1] + ")"
            name2 = "-" + y_method + "(" + y_data + ")"
            if vip_method== "true":
                name3 = vip
            table_dict = {
                "column": [{"field": "metabolite", "filter": "false", "sort": "false", "title": "Metabolite", "type": "int"},
                           {"field": "fc", "filter": "false", "sort": "false", "title": name1,
                            "type": "int"},
                           {"field": "pvalue", "filter": "false", "sort": "false", "title": name2, "type": "int"},
                           ], "condition": {}}
            if vip_method== "true":
                table_dict["column"].append({"field": "vip_plsda", "filter": "false", "sort": "false", "title": name3, "type": "int"})
            table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))
            scatter_dict = {"category": "category", "data":["fc","pvalue",""],"name": "metabolite","size": "size","showlabel":"show_name","condition": {"type":"scatter"}}
            scatter_info = json.dumps(scatter_dict, sort_keys=True, separators=(',', ':'))
            axis_dict = {"data":[name1,name2]}
            axis_info = json.dumps(axis_dict, sort_keys=True, separators=(',', ':'))
            scatter_group_dict = {"data":list(set(scatter_group))}
            scatter_group_info = json.dumps(scatter_group_dict, sort_keys=True, separators=(',', ':'))
        try:
            self.create_db_table('volcano_detail', insert_data)
            self.update_db_record('volcano', main_id, status="end", main_id=main_id, stat_table=table_info,scatter_data=scatter_info,axis_data=axis_info,scatter_group=scatter_group_info)
        except Exception as e:
            self.bind_object.logger.error("导入volcano数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入volcano数据成功")

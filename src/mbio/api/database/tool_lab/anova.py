# -*- coding: utf-8 -*-


import os
import json
from collections import defaultdict
import datetime
from bson.objectid import ObjectId
from types import StringTypes
import pandas as pd
import glob
import re
from api_base import ApiBase

class Anova(ApiBase):
    def __init__(self, bind_object):
        super(Anova, self).__init__(bind_object)

    def add_detail_box(self,main_id, result_path,box_file,methor=None,method=None,correction=None,group_names=None):
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        data_list = []
        box_list = []
        significant_list = []
        pvalue_text_list = []
        result_files = glob.glob(result_path + '/*_result.xls')
        result_file = result_files[0]
        dtype_dic = {'subject_id': str,'subject_number': 'float'}
        rdata = pd.read_table(result_file, sep='\t',index_col=0,dtype = dtype_dic)
        merge_text = []
        merge_significant = []
        group_ishap = []

        #获取group列表
        if not group_names:
            group_names = []
            for col in rdata.columns:
                mat_result = re.match('^(.*)-mean$', col)
                if mat_result:
                    group_names.append(mat_result.group(1))

        #获取丰度排名
        sum_abund = rdata[[g+'-mean' for g in group_names]].apply(lambda x: sum(x),axis=1)
        rdata['rank'] = sum_abund.rank(ascending=False)
        rdata['p_rank'] = rdata['qvalue'].rank(ascending=True,method='first')

        bdata = pd.read_table(box_file, sep='\t')
        b_list = bdata.columns.tolist()
        #r_list = rdata.columns.tolist()
        #rdata.rename(columns={r_list[0]: "name"}, inplace=True)
        bdata.rename(columns={b_list[0]: "name"}, inplace=True)
        final_data = pd.merge(rdata, bdata, on="name")
        #final_data = rdata.join(bdata)
        final_data.fillna("NA", inplace=True)

        final_data.set_index(['name'], inplace=True)
        final_data["name"] = final_data.index
        for index in final_data.index:
            cur_index_value = final_data.loc[index]
            rank = cur_index_value['rank']
            insert_data = {
                "anova_id" : main_id,
                "name": index,
                "p_value" : cur_index_value['pvalue'],
                "q_value" : cur_index_value['qvalue'],
                "rank" : rank,
                "p_rank" : cur_index_value["p_rank"]
            }

            sig_data = {
                "anova_id" : main_id,
                "name": index,
                "value" :  cur_index_value['qvalue'],
                "x" : index+"_"+group_names[0]+"|"+index+"_"+group_names[-1],
                "type": "significant_column",
                "factor": "0.05|0.01|0.001"
            }
            significant_list.append(sig_data)

            text_data = {
                "anova_id": main_id,
                "name": index,
                "text": str(cur_index_value['qvalue']),
                "x": index,
                "type": "text"
            }
            pvalue_text_list.append(text_data)

            for g in group_names:
                g_mean = g+'-mean'
                g_sd = g+'-sd'

                g_id_mean_sd = ','.join([str(cur_index_value[g_mean]),str(cur_index_value[g_sd])])
                insert_data[g] = float(g_id_mean_sd.split(',')[0])
                insert_data[g+'_sd'] = float(g_id_mean_sd.split(',')[1])
                insert_data[g+'_sd_bk'] = insert_data[g+'_sd']

                q1 = final_data[g+'--25%'][index]
                q3 = final_data[g+'--75%'][index]
                mean = final_data[g+'--50%'][index]
                #min = final_data[g+'--min'][index]
                #max = final_data[g+'--max'][index]
                sub_min = final_data[g+'--min'][index]
                sub_max = final_data[g+'--max'][index]
                abnormal = final_data[g+'--abnormal'][index]
                if ':' not in str(abnormal):
                    abnormals = []
                else:
                    abnormals = abnormal.split(';')

                box2_data = {
                    "anova_id": main_id,
                    "name": index,
                    "category": g,
                    "rank": rank,
                    "q1": q1,
                    "q3": q3,
                    "median": mean,
                    "min": sub_min,
                    "max": sub_max,
                    "type": "box",
                    "x": index+"_"+g

                }
                box_list.append(box2_data)

                for ab in abnormals:
                    abv = ab.split(':')
                    scater_data = {
                        "anova_id": main_id,
                        "x": index,
                        "y": float(abv[1]),
                        "category": g,
                        "name": index+'_'+g,
                        "type": "scatter"
                    }
                    box_list.append(scater_data)

            insert_data['relation_sign'] = cur_index_value['name']
            insert_data["merge_scatter"] = {}
            for ci_file in os.listdir(result_path):
                if methor in ci_file:
                    #self.bind_object.logger.info("ci_file:%s" % ci_file)
                    with open(result_path + "/" + ci_file) as f:
                        tmp = f.readlines()
                        for x in tmp[1:]:
                            effectsize = 0
                            lowerci = 0
                            upperci = 0
                            low_name = ""
                            upper_name = ""
                            if x.split("\t")[0] == cur_index_value['name']:
                                for xx in range(len(tmp[0].strip("\n").split("\t")) - 1):
                                    if tmp[0].strip("\n").split("\t")[xx + 1].endswith("effectsize"):
                                        effectsize = float(x.strip().split("\t")[xx + 1])
                                        insert_data[tmp[0].strip("\n").split("\t")[xx + 1]] = effectsize
                                        insert_data["merge_scatter"][
                                            tmp[0].strip("\n").split("\t")[xx + 1].split("_effectsize")[0]] = float(
                                            x.strip().split("\t")[xx + 1])
                                        group_ishap.append(tmp[0].strip("\n").split("\t")[xx + 1].split("_effectsize")[0])
                                    elif tmp[0].strip("\n").split("\t")[xx + 1].endswith("lowerCI"):
                                        low_name = tmp[0].strip("\n").split("\t")[xx + 1]
                                        lowerci = float(x.strip().split("\t")[xx + 1])
                                    elif tmp[0].strip("\n").split("\t")[xx + 1].endswith("upperCI"):
                                        upper_name = tmp[0].strip("\n").split("\t")[xx + 1]
                                        upperci = float(x.strip().split("\t")[xx + 1])
                                    else:
                                        insert_data[tmp[0].strip("\n").split("\t")[xx + 1]] = x.strip().split("\t")[
                                            xx + 1]
                                        if tmp[0].strip("\n").split("\t")[xx + 1].endswith("pvalue"):
                                            sig_data = {
                                                "anova_id": main_id,
                                                "name": cur_index_value['name'],
                                                "value": float( x.strip().split("\t")[xx + 1].split(" ")[-1]) if " " in str(
                                                    x.strip().split("\t")[xx + 1]) else round(
                                                    float(x.strip().split("\t")[xx + 1]), 4),
                                                "x": cur_index_value['name'] + "_" +
                                                     tmp[0].strip("\n").split("\t")[xx + 1].split("_pvalue")[
                                                         0].split("-")[0] + "|" +
                                                     cur_index_value['name'] + "_" +
                                                     tmp[0].strip("\n").split("\t")[xx + 1].split("_pvalue")[
                                                         0].split("-")[1],
                                                "type": "significant",
                                                "factor": "0.05|0.01|0.001"
                                            }
                                            significant_list.append(sig_data)
                                    insert_data[low_name] = float(abs(effectsize - lowerci))
                                    insert_data[upper_name] = float(abs(upperci - effectsize))

                # box1_data = {
                #     "anova_id" : main_id,
                #     "name": index,
                #     "category": g,
                #     "rank" : rank,
                #     "q1" : q1,
                #     "q3" : q3,
                #     "median" : mean,
                #     "min" : min,
                #     "max" : max,
                #     "type" : "box1"
                #
                # }
                # box_list.append(box1_data)
            data_list.append(insert_data)
        mul_result_files = glob.glob(result_path + '/anova_*.xls')
        if len(mul_result_files) > 1:
            for mul_result_file in mul_result_files:
                if "anova_boxfile.xls" in mul_result_file or "anova_result.xls" in mul_result_file:
                    pass
                else:
                    mul_df1 = pd.read_table(mul_result_file, sep='\t', index_col=0)
                    mul_col1 = mul_df1.columns.tolist()
                    mul_p1 = [i for i in mul_col1 if i.find('pvalue') != -1]
                    new_col1 = []
                    for n1 in mul_df1[mul_p1[0]].tolist():
                        try:
                            float(n1.split(" ")[-1])
                        except:
                            new_col1.append(0.1)
                        else:
                            new_col1.append(float(n1.split(" ")[-1]))
                    rdata[mul_p1[0]] = new_col1
                    for sample in mul_df1[mul_p1[0]].index:
                        pvalue = mul_df1[mul_p1[0]].loc[sample]
                        sig_data5 = {
                            "anova_id": main_id,
                            "name": sample,
                            "value": round(float(pvalue.split(" ")[-1]) if " " in str(pvalue) else pvalue,4),
                            "x": "_".join(os.path.basename(mul_result_file).split(".xls")[0].split("_")[2:]),
                            "type": "merge_significant",
                            "factor": "0.05|0.01|0.001"
                        }
                        merge_significant.append(sig_data5)
                        pvalue = mul_df1[mul_p1[0]].loc[sample]
                        sig_data6 = {
                            "anova_id": main_id,
                            "name": sample,
                            "text": pvalue if " " in str(pvalue) else str(round(float(pvalue),4)), #round(float(pvalue.split(" ")[-1]) if " " in str(pvalue) else pvalue,4), #round(float(pvalue.split(" ")[-1]) if " " in str(pvalue) else pvalue,4)
                            "x": "_".join(os.path.basename(mul_result_file).split(".xls")[0].split("_")[2:]),
                            "type": "merge_text",
                        }
                        merge_text.append(sig_data6)



        try:
            collection = self.db["anova_detail"]
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.info("导入anova_detail数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入anova_detail数据成功")

        self.bind_object.logger.info("merge_text:%s" % merge_text)
        self.bind_object.logger.info("merge_significant:%s" % merge_significant)
        try:
            self.db["anova_box"].insert_many(box_list)
            if merge_text:
                self.db["anova_box"].insert_many(merge_text)
            if merge_significant:
                self.db["anova_box"].insert_many(merge_significant)
        except Exception as e:
            self.bind_object.logger.info("导入anova_box数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入anova_box数据成功")

        insert_pvalue = []

        for file in os.listdir(result_path):
            if '_result.xls' in file or '_boxfile.xls' in file or  'box_result.xls' in file:
                continue

            with open(result_path+'/'+file) as fr:
                diff = fr.readline().rstrip().split('\t')[4].rsplit('_pvalue',1)[0]

                for line in fr:
                    spline = line.strip().split('\t')
                    name = spline[0]
                    pvalue = spline[4]
                    eff_size = float(spline[1])
                    low_ci = float(spline[2])
                    up_ci = float(spline[3])

                    tmp = {
                        "diff" : diff,
                        'pvalue' : pvalue,
                        'anova_id' : main_id,
                        'name' : name,
                        'eff_size' :eff_size,
                        'low_ci' :low_ci,
                        'up_ci': up_ci
                    }
                    insert_pvalue.append(tmp)


                    sig_data = {
                        "anova_id" : main_id,
                        "name": name,
                        "value" :  pvalue,
                        "x" : diff.replace('-','|'),  #组名没有-
                        "type": "significant_single",
                        "factor": "0.05|0.01|0.001"
                    }
                    significant_list.append(sig_data)

        try:
            collection = self.db["anova_relation"]
            collection.insert_many(insert_pvalue)
            self.db["anova_box"].insert_many(significant_list)
            self.db["anova_box"].insert_many(pvalue_text_list)
        except Exception as e:
            self.bind_object.logger.info("导入anova_relation数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入anova_relation数据成功")


        #更新主表
        main_collection = self.db["anova"]

        column = [
            {"field": "name", "type": "string", "sort": "false", "title": "Name"},
            {"field": "p_value", "type": "float", "sort": "false", "title": "Pvalue"}
        ]
        if correction != "none":
            column.append({"field": "q_value", "type": "float", "sort": "false", "title": "Corrected Pvalue"})

        for g in group_names:
            if method == "T":
                column.append({"field": g, "type": "string", "sort": "false", "title": g+'-mean(%)'})
                column.append({"field": g + '_sd', "type": "string", "sort": "false", "title": g+'-SD(%)'})
            else:
                column.append({"field": g, "type": "string", "sort": "false", "title": g + '-mean'})
                column.append({"field": g + '_sd', "type": "string", "sort": "false", "title": g + '-SD'})

        data2 = []
        for g in group_names:
            data2.append([g + '_effectsize', g + '_lowerCI', g + '_upperCI'])
        ishape_data2 = {"name": "name", "data": data2, "category": "group", "condition": {}}

        table_data  = {'column':column,'condition': {}}

        data = []
        for g in group_names:
            data.append(g)
        column_data = {"name":"name","data":data}

        #工型数据
        data = []
        for g in group_names:
            data.append([g,g+'_sd',g+'_sd_bk'])
        ishape_data = {"name":'name', "data":data}

        #箱线图
        box_data = {"name":"name", "x":"x","condition":{"type":"box"}}

        #箱线图散点
        scatter_data ={"name":"name","data":["x", "y"],"category":"category", "condition":{"type":"scatter"}}

        #单变量柱形图
        data = []
        for g in group_names:
            data.append(g)
        column_data1 = column_data

        data3 = []
        for g in set(group_ishap):
            data3.append([ g + '_effectsize', g + '_lowerCI', g + '_upperCI'])
        merge_ishape_data = {"name": "name", "data": data3}

        merge_text_data = {"name": "name", "condition": {}}

        merge_significant_data = {"name": "name", "value": "value", "condition": {}}

        ishape_data1 = ishape_data

        box_data1 = box_data
        scatter_data1 = scatter_data

        multi_significant_data = {"name": "name","value":"value", "condition":{"type": "significant_column"}}
        multi_significant_data1 = {"name": "name","value":"value", "condition":{"type": "significant"}}
        single_significant_data = {"name": "name","value":"value", "condition":{"type": "significant"}}
        pvalue_text_data = {"name":"name","condition":{"type":"text"}}

        update_info = {
            "main_id": main_id,
            "table_data"  : json.dumps(table_data),
            #"single_column_data" : json.dumps(column_data),
            #"single_ishape_data" : json.dumps(ishape_data),
            #"single_box_data" : json.dumps(box_data),
            #"single_scatter_data" : json.dumps(scatter_data),
            "multi_column_data" : json.dumps(column_data1),
            "multi_ishape_data" : json.dumps(ishape_data1),
            "multi_box_data" : json.dumps(box_data1),
            "multi_scatter_data" : json.dumps(scatter_data1),
            "box_significant_data": json.dumps(multi_significant_data1),
            "multi_significant_data" : json.dumps(multi_significant_data),
            #"single_significant_data" : json.dumps(single_significant_data),
            "merge_ishape_data": json.dumps(merge_ishape_data),
            "merge_text_data": json.dumps(merge_text_data),
            "merge_significant_data": json.dumps(merge_significant_data),
            "pvalue_text_data" : json.dumps(pvalue_text_data),
            "merge_scatter_data": json.dumps(ishape_data2),
        }
        main_collection.update({"_id": main_id},{"$set": update_info})
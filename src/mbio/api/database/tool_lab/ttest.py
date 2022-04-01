# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20180608
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import glob
import re
import pandas as pd
from api_base import ApiBase
import json


class Ttest(ApiBase):
    def __init__(self, bind_object):
        super(Ttest, self).__init__(bind_object)

    def add_detail_box(self, main_id, result_path, box_file, method="T",mul_test="fdr"):
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        data_list = []
        box_list = []
        significant_list = []
        pvalue_text_list = []
        result_files = glob.glob(result_path + '/*_result.xls')
        cifiles = glob.glob(result_path + '/*_CI.xls')
        result_file = result_files[0]
        ci_file = cifiles[0]
        rdata = pd.read_table(result_file, sep='\t',index_col=0)

        #获取group列表
        group_names = []
        for col in rdata.columns:
            mat_result = re.match('^(.*)-mean$', col)
            if mat_result:
                group_names.append(mat_result.group(1))

        #获取丰度排名
        sum_abund = rdata[[g+'-mean' for g in group_names]].apply(lambda x: sum(x),axis=1)
        rdata['rank'] = sum_abund.rank(ascending=False,method='first')
        rdata['p_rank'] = rdata['pvalue'].rank(ascending=True,method='first')

        bdata = pd.read_table(box_file, sep='\t',index_col=0)
        final_data = rdata.join(bdata)
        ci_data = pd.read_table(ci_file, sep='\t',index_col=0)
        final_data = final_data.join(ci_data)
        final_data.fillna(0, inplace=True)
        for index in final_data.index:
            cur_index_value = final_data.loc[index]
            rank = cur_index_value['rank']
            p_rank = cur_index_value['p_rank']
            g0 = group_names[0]
            g1 = group_names[1]
            g0_mean = cur_index_value[g0+'-mean']
            g1_mean = cur_index_value[g1+'-mean']
            if g0_mean!=0 and g1_mean!=0:
                fc = round(float(g0_mean)/g1_mean,5)
                dfc = round(float(g1_mean)/g0_mean,5)
                mfc = max([fc,dfc])
            else:
                fc = '-'
                dfc = '-'
                mfc = '-'
            insert_data = {
                "ttest_id" : main_id,
                "name": index,
                "p_value" : cur_index_value['pvalue'],
                "q_value" : cur_index_value['corrected_pvalue'],
                "rank" : rank,
                'p_rank': p_rank,
                "effect" : cur_index_value['effectsize'],
                "low_ci" : cur_index_value['lowerCI'],
                "up_ci" : cur_index_value['upperCI'],
                g0  : g0_mean,
                g1  : g1_mean,
                g0+'_sd' : cur_index_value[g0+'-sd'],
                g0+'_sd_bk' : cur_index_value[g0+'-sd'],  #画工型用
                g1+'_sd' : cur_index_value[g1+'-sd'],
                g1+'_sd_bk' : cur_index_value[g1+'-sd'], #画工型用
                "fc" :fc,
                'd_fc' : dfc,
                'mix_fc' : mfc,
                "effectsize":cur_index_value['effectsize'],
                "lowerCI": abs(float(cur_index_value['effectsize']) - float(cur_index_value['lowerCI'])),
                "upperCI":abs(float(cur_index_value['effectsize']) - float(cur_index_value['upperCI']))
            }

            sig_data = {
                "ttest_id" : main_id,
                "name": index,
                "value" :  round(float(cur_index_value['corrected_pvalue']),4),
                "x" : index+"_"+group_names[0]+"|"+index+"_"+group_names[1],
                "type": "significant",
                "factor": "0.05|0.01|0.001"
            }
            significant_list.append(sig_data)

            text_data = {
                "ttest_id" : main_id,
                "name": index,
                "x" : index,
                "text" : '%.4g' % float(cur_index_value['corrected_pvalue']),
                "type" : "text"
            }
            pvalue_text_list.append(text_data)

            sig_data2 = {
                "ttest_id" : main_id,
                "name": index,
                "value" :  round(float(cur_index_value['corrected_pvalue']),4),
                "x" : "{0}|{1}".format(g0,g1),
                "type": "significant_single",
                "factor": "0.05|0.01|0.001"
            }
            significant_list.append(sig_data2)

            for g in group_names:
                # g_mean = g+'-mean'
                # g_sd = g+'-sd'
                #
                # g_id_mean_sd = ','.join([str(cur_index_value[g_mean]),str(cur_index_value[g_sd])])
                # insert_data[g+'_mean'] = float(g_id_mean_sd.split(',')[0])
                # insert_data[g+'_sd'] = float(g_id_mean_sd.split(',')[1])

                q1 = final_data[g+'--25%'][index]
                q3 = final_data[g+'--75%'][index]
                mean = final_data[g+'--50%'][index]
                min_v = final_data[g+'--min'][index]
                max_v = final_data[g+'--max'][index]
                sub_min = final_data[g+'--min'][index]
                sub_max = final_data[g+'--max'][index]
                abnormal = final_data[g+'--abnormal'][index]

                if ':' not in str(abnormal):
                    abnormals = []
                else:
                    abnormals = abnormal.split(';')

                # box1_data = {
                #     "ttest_id" : main_id,
                #     "name": index,
                #     "category": g,
                #     "rank" : rank,
                #     "q1" : q1,
                #     "q3" : q3,
                #     "median" : mean,
                #     "min" : min_v,
                #     "max" : max_v,
                #     "type" : "box1"
                #
                # }
                # box_list.append(box1_data)

                box2_data = {
                    "ttest_id" : main_id,
                    "name": index,
                    "category": g,
                    "rank" : rank,
                    "q1":q1,
                    "q3": q3,
                    "median": mean,
                    "min":sub_min,
                    "max":sub_max,
                    "type" : "box",
                    "x": index+"_"+g
                }
                box_list.append(box2_data)

                for ab in abnormals:
                    abv = ab.split(':')
                    scater_data = {
                        "ttest_id" : main_id,
                        "x": index,
                        "y":  float(abv[1]),
                        "category": g,
                        "name" : index+"_"+g,
                        "type" : "scatter"
                    }
                    box_list.append(scater_data)

            data_list.append(insert_data)

        try:
            collection = self.db["ttest_detail"]
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.info("导入ttest_detail数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入ttest_detail数据成功")

        try:
            self.db["ttest_box"].insert_many(box_list)
        except Exception as e:
            self.bind_object.logger.info("导入ttest_box数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入ttest_box数据成功")

        #插入星号数据
        try:
            self.db["ttest_box"].insert_many(significant_list)
            self.db["ttest_box"].insert_many(pvalue_text_list)
        except Exception as e:
            self.bind_object.logger.info("星号数据导入ttest_box数据出错:%s" % e)
        else:
            self.bind_object.logger.info("星号数据导入ttest_box数据成功")

        #更新主表字段
        self.update_main(group_names,main_id,method,mul_test)


    def update_main(self, group_names,main_id,method="T",mul_test="fdr"):
        if method == "T":
            column = [
                {"field": "name", "type": "string", "sort": "false", "title": "Name"},
                {"field": "fc", "type": "string", "sort": "false",
                 "title": "Fold Change(%s/%s)" % (group_names[0], group_names[1])},
                {"field": group_names[0], "type": "string", "sort": "false", "title": group_names[0] + '-Mean(%)'},
                {"field": group_names[1], "type": "string", "sort": "false", "title": group_names[1] + '-Mean(%)'},
                {"field": group_names[0] + "_sd", "type": "string", "sort": "false", "title": group_names[0] + '-Sd(%)'},
                {"field": group_names[1] + "_sd", "type": "string", "sort": "false", "title": group_names[1] + "-Sd(%)"},
                {"field": "low_ci", "type": "string", "sort": "false", "title": "Lower ci"},
                {"field": "up_ci", "type": "string", "sort": "false", "title": "Upper ci"},
                {"field": "effect", "type": "string", "sort": "false", "title": "Effect size"},
                {"field": "p_value", "type": "float", "sort": "false", "title": "P-value"}
            ]
        else:
            column = [
                {"field": "name", "type": "string", "sort": "false", "title": "Name"},
                {"field": "fc", "type": "string", "sort": "false", "title": "Fold Change(%s/%s)" % (group_names[0], group_names[1])},
                {"field": group_names[0], "type": "string", "sort": "false", "title": group_names[0] + '-Mean'},
                {"field": group_names[1], "type": "string", "sort": "false", "title": group_names[1] + '-Mean'},
                {"field": group_names[0] + "_sd", "type": "string", "sort": "false", "title": group_names[0] + '-Sd'},
                {"field": group_names[1] + "_sd", "type": "string", "sort": "false", "title": group_names[1] + "-Sd"},
                {"field": "low_ci", "type": "string", "sort": "false", "title": "Lower ci"},
                {"field": "up_ci", "type": "string", "sort": "false", "title": "Upper ci"},
                {"field": "effect", "type": "string", "sort": "false", "title": "Effect size"},
                {"field": "p_value", "type": "float", "sort": "false", "title": "P-value"}
            ]
        if mul_test != "none":
            column.append({"field": "q_value", "type": "float", "sort": "false", "title": "corrected pvalue"})
        table_data  = {'column':column,'condition': {}}

        data = []
        for g in group_names:
            data.append(g)
        column_data = {"name":"name","data":data}

        #工型数据
        data = []
        for g in group_names:
            data.append([g, g+'_sd', g+'_sd_bk'])
        ishape_data = {"name": "name", "data": data}

        # 工型数据2
        # data = []
        # for g in group_names:
        #    data.append([g, g+'_sd', g+'_sd_bk'])
        ishape_data2 = {"name": "name", "data": ["effectsize", "lowerCI", "upperCI"]}

        #箱线图
        box_data = {"name":"name", "x":"x", "condition":{"type":"box"}}

        #箱线图散点
        scatter_data ={"name":"name","data":["x", "y"],"category":"category", "condition":{}}

        #多变量星号
        multi_significant_data = {"name": "name","value":"value", "condition":{"type": "significant"}}

        #单变量柱形图
        data = []
        for g in group_names:
            data.append(g)
        column_data1 = column_data

        ishape_data1 = ishape_data

        box_data1 = box_data
        scatter_data1 = scatter_data

        #单变量星号
        single_significant_data = {"name": "name","value":"value", "condition":{"type": "significant"}}

        #pvalue text
        pvalue_text = {"name": "name","condition":{"type":"text"}}

        update_info = {
            "main_id": main_id,
            "table_data"  : json.dumps(table_data),
            "multi_column_data" : json.dumps(column_data),
            "multi_ishape_data" : json.dumps(ishape_data),
            "multi_box_data" : json.dumps(box_data),
            "multi_scatter_data" : json.dumps(scatter_data),
            "multi_significant_data" : json.dumps(multi_significant_data),
            "merge_ishape_data": json.dumps(ishape_data2),
            #"single_column_data" : json.dumps(column_data1),
            #"single_ishape_data" : json.dumps(ishape_data1),
            #"single_box_data" : json.dumps(box_data1),
            #"single_scatter_data" : json.dumps(scatter_data1),
            #"single_significant_data" : json.dumps(single_significant_data),
            "pvalue_text_data" : json.dumps(pvalue_text)
        }
        main_collection = self.db['ttest']
        main_collection.update({"_id": main_id},{"$set": update_info})

# !usr/bin/python
# -*- coding: utf-8 -*-
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
from types import StringTypes


class MannTest(ApiBase):
    def __init__(self, bind_object):
        super(MannTest, self).__init__(bind_object)

    def add_mann_detail(self, main_id, result_path, box_file, group_names,method="T",compare_type="multi",kru_methor="tukeykramer",correction=None):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        data_list = []
        box_list = []
        significant_list = []
        result_files = glob.glob(result_path + '/*_result.xls')
        result_file = result_files[0]
        rdata = pd.read_table(result_file, sep='\t', index_col=0)
        #print (rdata)
        merge_text = []
        merge_significant = []
        group_ishap = []

        # 获取group列表  这里改成根据输入参数获取group列表
        # group_names = []
        #for col in rdata.columns:
            # mat_result = re.match('^(.*)-mean$', col)
            # if mat_result:
            #     group_names.append(mat_result.group(1))

            # 多组比较时两两比较计算Pvalue值
        ## post-hoc图相关数据
        mul_result_files = glob.glob(result_path + '/kru_H_*.xls')
        if len(mul_result_files) > 1:
            for mul_result_file in mul_result_files:
                mul_df1 = pd.read_table(mul_result_file, sep='\t', index_col=0)
                mul_col1 = mul_df1.columns.tolist()
                mul_p1 = [i for i in mul_col1 if i.find('_pvalue') != -1]
                new_col1 = []
                for n1 in mul_df1[mul_p1[0]].tolist():
                    try:
                        float(n1.split(" ")[-1]) if " " in str(n1) else float(n1)
                    except:
                        new_col1.append(0.1)
                    else:
                        new_col1.append(float(n1.split(" ")[-1]) if " " in str(n1) else float(n1))
                rdata[mul_p1[0]] = new_col1
                for sample in mul_df1[mul_p1[0]].index:
                    pvalue = mul_df1[mul_p1[0]].loc[sample]
                    sig_data5 = {
                        "mann_id": main_id,
                        "name": sample,
                        "value": round(float(pvalue.split(" ")[-1]) if " " in str(pvalue) else pvalue,4),
                        "x": "_".join(os.path.basename(mul_result_file).split(".xls")[0].split("_")[3:]),
                        "type": "merge_significant",
                        "factor": "0.05|0.01|0.001"
                    }
                    merge_significant.append(sig_data5)
                    pvalue = mul_df1[mul_p1[0]].loc[sample]
                    sig_data6 = {
                        "mann_id": main_id,
                        "name": sample,
                        "text": pvalue if " " in str(pvalue) else str(round(float(pvalue),4)), # round(float(pvalue.split(" ")[-1]) if " " in str(pvalue) else pvalue,4)
                        "x": "_".join(os.path.basename(mul_result_file).split(".xls")[0].split("_")[3:]),
                        "type": "merge_text",
                    }
                    merge_text.append(sig_data6)
        elif len(mul_result_files) == 1:
            mul_result_file = mul_result_files[0]
            aa = open(mul_result_file)
            data = aa.readlines()
            for x in data[1:]:
                for y in range(len(data[0].strip().split("\t"))):
                    name = data[0].strip().split("\t")[y]
                    if name.endswith("_pvalue"):
                        sig_data5 = {
                            "mann_id": main_id,
                            "name": x.strip().split("\t")[0],
                            "value":  round(float(x.strip().split("\t")[y+1]),4) if x.strip().split("\t")[y+1] != "NaN" else 1.0,
                            "x": name.split("_pvalue")[0],
                            "type": "merge_significant",
                            "factor": "0.05|0.01|0.001"
                        }
                        merge_significant.append(sig_data5)
                        sig_data6 = {
                            "mann_id": main_id,
                            "name": x.strip().split("\t")[0],
                            "text": ("%.4g" % float(x.strip().split("\t")[y+1])) if x.strip().split("\t")[y+1] != "NaN" else "1.0",
                            "x": name.split("_pvalue")[0],
                            "type": "merge_text",
                        }
                        merge_text.append(sig_data6)
        else:
            aa = open(result_file)
            data = aa.readlines()
            for x in data[1:]:
                for y in range(len(data[0].strip().split("\t"))):
                    name = data[0].strip().split("\t")[y]
                    if name.endswith("corrected_pvalue"):
                        """
                        sig_data5 = {
                            "mann_id": main_id,
                            "name": x.strip().split("\t")[1],
                            "value": round(float(x.strip().split("\t")[y+1]) if x.strip().split("\t")[y+1] != "NaN" else 1.0,4),
                            "x": group_names[0] + "-" +group_names[1],
                            "type": "merge_significant",
                            "factor": "0.05|0.01|0.001"
                        }
                        merge_significant.append(sig_data5)
                        """
                        sig_data6 = {
                            "mann_id": main_id,
                            "name": x.strip().split("\t")[1],
                            "text":  '%.4g' % float(x.strip().split("\t")[y+1]) if x.strip().split("\t")[y+1] != "NaN" else 1.0,
                            "x": group_names[0] + "-" +group_names[1],
                            "type": "merge_text",
                        }
                        merge_text.append(sig_data6)

        # 获取丰度排名
        sum_abund = rdata[[g+'-mean' for g in group_names]].apply(lambda x: sum(x), axis=1)
        rdata['rank'] = sum_abund.rank(ascending=False, method='first')
        rdata['p_rank'] = rdata['pvalue'].rank(ascending=True, method='first')
        rdata['cp_rank'] = rdata['corrected_pvalue'].rank(ascending=True, method='first')

        bdata = pd.read_table(box_file, sep='\t')
        b_list = bdata.columns.tolist()
        r_list = rdata.columns.tolist()
        rdata.rename(columns={r_list[0]: "name"}, inplace=True)
        bdata.rename(columns={b_list[0]: "name"}, inplace=True)
        final_data = pd.merge(rdata, bdata, on="name")

        for index in final_data.index:
            cur_index_value = final_data.loc[index]
            rank = cur_index_value['rank']
            p_rank = cur_index_value['p_rank']
            cp_rank = cur_index_value['cp_rank']
            g0 = group_names[0]
            g1 = group_names[1]
            g0_mean = cur_index_value[g0+'-mean']
            g1_mean = cur_index_value[g1+'-mean']
            if g0_mean != 0 and g1_mean != 0:
                fc = float(g0_mean)/g1_mean
            else:
                fc = ''
            insert_data = {
                "mann_id": main_id,
                "name": cur_index_value['name'],
                "p_value": cur_index_value['pvalue'],
                "q_value": cur_index_value['corrected_pvalue'],
                "rank": rank,
                'p_rank': p_rank,
                'cp_rank': cp_rank,
                g0: round(float(g0_mean),4),
                g1: round(float(g1_mean),4),
                g0+'_sd': round(float(cur_index_value[g0+'-sd']),4),
                g0+'_sd_bk': round(float(cur_index_value[g0+'-sd']),4),  # 画工字型用
                g1+'_sd': round(float(cur_index_value[g1+'-sd']),4),
                g1+'_sd_bk': round(float(cur_index_value[g1+'-sd']),4),  # 画工字型用
            }

            if len(group_names) > 2:
                for g in group_names[2:]:
                    g_mean = cur_index_value[g + '-mean']
                    insert_data[g] = round(float(g_mean),4)
                    insert_data[g+'_sd'] = round(float(cur_index_value[g+'-sd']),4)
                    insert_data[g+'_sd_bk'] = round(float(cur_index_value[g+'-sd']),4)
                insert_data['relation_sign'] = cur_index_value['name']
                insert_data["merge_scatter"] = {}

                for ci_file in os.listdir(result_path):
                    if kru_methor in ci_file:
                        with open(result_path+"/"+ci_file) as f:
                            tmp = f.readlines()
                            for x in tmp[1:]:
                                effectsize = 0
                                lowerci = 0
                                upperci = 0
                                low_name = ""
                                upper_name = ""
                                if x.split("\t")[0] == cur_index_value['name']:
                                    for xx in range(len(tmp[0].strip("\n").split("\t"))-1):
                                        if tmp[0].strip("\n").split("\t")[xx+1].endswith("effectsize"):
                                            effectsize = float(x.strip().split("\t")[xx+1])
                                            insert_data[tmp[0].strip("\n").split("\t")[xx + 1]] = float(x.strip().split("\t")[xx + 1])
                                            insert_data["merge_scatter"][tmp[0].strip("\n").split("\t")[xx + 1].split("_effectsize")[0]] = float(x.strip().split("\t")[xx + 1])
                                            group_ishap.append(tmp[0].strip("\n").split("\t")[xx + 1].split("_effectsize")[0])
                                        elif tmp[0].strip("\n").split("\t")[xx+1].endswith("lowerCI"):
                                            low_name = tmp[0].strip("\n").split("\t")[xx + 1]
                                            lowerci = float(x.strip().split("\t")[xx+1])
                                        elif tmp[0].strip("\n").split("\t")[xx+1].endswith("upperCI"):
                                            upper_name = tmp[0].strip("\n").split("\t")[xx + 1]
                                            upperci = float(x.strip().split("\t")[xx+1])
                                        else:
                                            insert_data[tmp[0].strip("\n").split("\t")[xx + 1]] = x.strip().split("\t")[xx + 1]
                                            if tmp[0].strip("\n").split("\t")[xx+1].endswith("_pvalue"):
                                                sig_data = {
                                                    "mann_id": main_id,
                                                    "name": cur_index_value['name'],
                                                    "value": float(
                                                        x.strip().split("\t")[xx + 1].split(" ")[-1]) if " " in str(
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
                                        insert_data[low_name] = float(abs(effectsize-lowerci))
                                        insert_data[upper_name] = float(abs(upperci - effectsize))
                                        #insert_data[low_name] = lowerci
                                        #insert_data[upper_name] = upperci

            else:
                # 两组比较时才有差异倍数  FC
                insert_data["fc"] = fc
                insert_data["effect"] = cur_index_value['Effectsize']
                insert_data["low_ci"] = cur_index_value['Lower ci']
                insert_data["up_ci"] = cur_index_value['Upper ci']
                insert_data["effectsize"] = cur_index_value['Effectsize']
                insert_data["lowerCI"] =  abs(float(cur_index_value['Effectsize']) - float(cur_index_value['Lower ci']))
                insert_data["upperCI"] = abs(float(cur_index_value['Effectsize']) - float(cur_index_value['Upper ci']))
                #with open(result_path+"/compare_CI.xls") as f:
                #    tmp = f.readlines()
                #    for x in tmp[1:]:
                #        if x.split("\t")[0] == cur_index_value['name']:
                #            for xx in range(len(tmp[0].strip("\n").split("\t")) - 1):
                #                insert_data[tmp[0].strip("\n").split("\t")[xx + 1]] = x.strip().split("\t")[xx + 1]

                sig_data = {
                    "mann_id": main_id,
                    "name": cur_index_value['name'],
                    "value": round(float(cur_index_value['corrected_pvalue']), 4),
                    "x": cur_index_value['name'] + "_" + group_names[0] + "|" + cur_index_value['name'] + "_" +
                         group_names[-1],
                    "type": "merge_significant",
                    "factor": "0.05|0.01|0.001"
                }
                significant_list.append(sig_data)

                sig_data = {
                    "mann_id": main_id,
                    "name": cur_index_value['name'],
                    "value": round(float(cur_index_value['corrected_pvalue']), 4),
                    "x": cur_index_value['name'] + "_" + group_names[0] + "|" + cur_index_value['name'] + "_" +
                         group_names[-1],
                    "type": "significant",
                    "factor": "0.05|0.01|0.001"
                }
                significant_list.append(sig_data)

            sig_data = {
                "mann_id": main_id,
                "name": cur_index_value['name'],
                "value": round(float(cur_index_value['corrected_pvalue']), 4),
                "x": cur_index_value['name'] + "_" + group_names[0] + "|" + cur_index_value['name'] + "_" +
                     group_names[-1],
                "type": "significant_column",
                "factor": "0.05|0.01|0.001"
            }
            significant_list.append(sig_data)

            pvaue_data = {
                "mann_id": main_id,
                "name": cur_index_value['name'],
                "text": '%.4g' %float(cur_index_value['corrected_pvalue']),
                "x": cur_index_value['name'],
                "q_text": cur_index_value['corrected_pvalue'],
                "type": "text"
            }
            significant_list.append(pvaue_data)

            if len(group_names) > 2:
                sig_data2 = {
                    "mann_id": main_id,
                    "name": cur_index_value['name'],
                    "value": round(float(cur_index_value['corrected_pvalue']),4),
                    "x": '|'.join(group_names),
                    "type": "significant_single",
                    "factor": "0.05|0.01|0.001"
                }
                significant_list.append(sig_data2)
            else:
                sig_data2 = {
                    "mann_id": main_id,
                    "name": cur_index_value['name'],
                    "value":  round(float(cur_index_value['corrected_pvalue']),4),
                    "x": "{0}|{1}".format(g0, g1),
                    "type": "significant_single",
                    "factor": "0.05|0.01|0.001"
                }
                significant_list.append(sig_data2)
            # 这里还要考虑三组比较的情况

            for g in group_names:
                q1 = final_data[g+'--25%'][index]
                q3 = final_data[g+'--75%'][index]
                mean = final_data[g+'--50%'][index]
                sub_min = final_data[g+'--min'][index]
                sub_max = final_data[g+'--max'][index]
                if method == "T":
                    box_data = {
                        "mann_id": main_id,
                        "name": cur_index_value['name'],
                        "category": g,
                        "rank": rank,
                        'p_rank': p_rank,
                        'cp_rank': cp_rank,
                        "p_value": round(float(cur_index_value['pvalue']),4),
                        "q_value": round(float(cur_index_value['corrected_pvalue']),4),
                        "q1": round(float(q1),4),
                        "q3": round(float(q3),4),
                        "median": round(float(mean),4),
                        "min": round(float(sub_min),4),
                        "max": round(float(sub_max),4),
                        "type": "box",
                        "x": cur_index_value['name'] + "_" + g
                    }
                else:
                    box_data = {
                        "mann_id": main_id,
                        "name": cur_index_value['name'],
                        "category": g,
                        "rank": rank,
                        'p_rank': p_rank,
                        'cp_rank': cp_rank,
                        "p_value": round(float(cur_index_value['pvalue']),4),
                        "q_value": round(float(cur_index_value['corrected_pvalue']),4),
                        "q1": q1,
                        "q3": q3,
                        "median": mean,
                        "min": sub_min,
                        "max": sub_max,
                        "type": "box",
                        "x": cur_index_value['name'] + "_" + g
                    }
                if len(group_names) == 2:
                    box_data["fc"] = fc
                box_list.append(box_data)

            data_list.append(insert_data)

        collection = self.db["mann_test_detail"]
        collection.insert_many(data_list)
        self.db["mann_test_box"].insert_many(box_list)
        self.db["mann_test_box"].insert_many(significant_list)
        if merge_text:
            self.db["mann_test_box"].insert_many(merge_text)
        if merge_significant:
            self.db["mann_test_box"].insert_many(merge_significant)

        # 更新主表字段
        self.update_main(group_names, main_id,method,set(group_ishap),correction)

    def update_main(self, group_names, main_id,method,group_ishap,correction):
        if len(group_names) == 2:
            if method == "T":
                column = [
                    {"field": "name", "type": "string", "sort": "false", "title": "Name"},
                    {"field": "fc", "type": "string", "sort": "false",
                     "title": "Fold Change(%s/%s)" % (group_names[0], group_names[1])},
                    {"field": group_names[0], "type": "float", "sort": "false", "title": group_names[0] + '-Mean(%)'},
                    {"field": group_names[0] + '_sd', "type": "float", "sort": "false", "title": group_names[0] + '-Sd(%)'},
                    {"field": group_names[1], "type": "float", "sort": "false", "title": group_names[1] + '-Mean(%)'},
                    {"field": group_names[1] + '_sd', "type": "float", "sort": "false", "title": group_names[1] + '-Sd(%)'},
                    {"field": "low_ci", "type": "string", "sort": "false", "title": "Lower ci"},
                    {"field": "up_ci", "type": "string", "sort": "false", "title": "Upper ci"},
                    {"field": "effect", "type": "string", "sort": "false", "title": "Effectsize"},
                    {"field": "p_value", "type": "float", "sort": "false", "title": "P-value"},
                ]
            else:
                column = [
                    {"field": "name", "type": "string", "sort": "false", "title": "Name"},
                    {"field": "fc", "type": "string", "sort": "false",
                     "title": "Fold Change(%s/%s)" % (group_names[0], group_names[1])},
                    {"field": group_names[0], "type": "float", "sort": "false", "title": group_names[0] + '-Mean'},
                    {"field": group_names[0] + '_sd', "type": "float", "sort": "false", "title": group_names[0] + '-Sd'},
                    {"field": group_names[1], "type": "float", "sort": "false", "title": group_names[1] + '-Mean'},
                    {"field": group_names[1] + '_sd', "type": "float", "sort": "false", "title": group_names[1] + '-Sd'},
                    {"field": "low_ci", "type": "string", "sort": "false", "title": "Lower ci"},
                    {"field": "up_ci", "type": "string", "sort": "false", "title": "Upper ci"},
                    {"field": "effect", "type": "string", "sort": "false", "title": "Effectsize"},
                    {"field": "p_value", "type": "float", "sort": "false", "title": "P-value"},
                ]
            if correction != "none":
                column.append({"field": "q_value", "type": "float", "sort": "false", "title": "corrected pvalue"})
        else:
            column = [
                {"field": "name", "type": "string", "sort": "false", "title": "Name"},
                {"field": "p_value", "type": "float", "sort": "false", "title": "P-value"}
            ]
            if correction != "none":
                column.append({"field": "q_value", "type": "float", "sort": "false", "title": "corrected pvalue"})
            if method == "T":
                for gn in group_names:
                    column.append({"field": gn, "type": "float", "sort": "false", "title": gn+'-Mean(%)'})
                    column.append({"field": gn+'_sd', "type": "float", "sort": "false", "title": gn+'-Sd(%)'})
            else:
                for gn in group_names:
                    column.append({"field": gn, "type": "float", "sort": "false", "title": gn+'-Mean'})
                    column.append({"field": gn+'_sd', "type": "float", "sort": "false", "title": gn+'-Sd'})
        # 结果表数据
        table_data = {'column': column}
        print column

        data = []
        for g in group_names:
            data.append(g)
        # 差异检验柱形图数据
        column_data = {"name": "name", "data": data}

        # 差异检验柱形图上工字型数据
        data = []
        for g in group_names:
            data.append([g, g+'_sd', g+'_sd_bk'])
        ishape_data = {"name": 'name', "data": data}

        # 差异检验柱形图上星号数据
        column_significant_data = {"name": "name", "value": "value", "condition": {"type": "significant_column"}}

        # 差异检验柱形图上pvalue数据
        column_text_data = {"name": "name", "condition": {"type": "text"}}

        # 箱型图数据
        box_data = {"name": "name", "x":"x","condition": {"type": "box"}}

        # 箱型图上星号数据
        box_significant_data = {"name": "name", "value": "value", "condition": {"type": "significant"}}

        # 箱型图上pvalue数据
        box_text_data = {"name": "name", "condition": {"type": "text"}}

        # 单物种比较柱状图数据
        data = []
        for g in group_names:
            data.append(g)
        column_data1 = column_data

        # 单物种比较柱状图上工字型数据
        ishape_data1 = ishape_data

        # 单物种比较箱型图数据
        box_data1 = box_data

        # 单物种比较箱型图上星号数据
        # single_significant_data = {"name": "name", "value": "value", "condition": {"type": "significant"}}

        # 置信区间柱形图or post-hoc检验图右边的工字数据
        if len(group_names) == 2:
            ishape_data2 = {"name": "name", "data": ["effectsize","lowerCI","upperCI"], "category": "group", "condition": {}}
        else:
            data2 = []
            for g in group_names:
                data2.append([g + '_effectsize', g + '_lowerCI', g + '_upperCI'])
            ishape_data2 = {"name": "name", "data": data2, "category": "group","condition": {}}

        # 工字图所用的一个字段
        if len(group_names) == 2:
            data3 = []
            #for g in group_names:
            #    data3.append([g, g + '_sd', g + '_sd_bk'])
            merge_ishape_data = {"name": "name", "data": ["effectsize","lowerCI","upperCI"]}
        else:
            data3 = []
            for g in group_ishap:
                data3.append([g + '_effectsize',g + '_lowerCI', g + '_upperCI'])
            merge_ishape_data = {"name": "name", "data": data3}

        merge_text_data = {"name": "name", "condition": {}}

        merge_significant_data = {"name": "name", "value": "value", "condition": {}}

        update_info = {
            "main_id": main_id,
            "table_data": json.dumps(table_data),
            "mul_column_data": json.dumps(column_data),
            "mul_ishape_data": json.dumps(ishape_data),
            "column_significant_data": json.dumps(column_significant_data),
            "column_text_data": json.dumps(column_text_data),
            "mul_box_data": json.dumps(box_data),
            "box_significant_data": json.dumps(box_significant_data),
            "box_text_data": json.dumps(box_text_data),
            "single_column_data": json.dumps(column_data1),
            "single_ishape_data": json.dumps(ishape_data1),
            "single_box_data": json.dumps(box_data1),
            "mul_scatter_data": json.dumps(ishape_data2),
            "merge_ishape_data":json.dumps(merge_ishape_data),
            "merge_text_data": json.dumps(merge_text_data),
            "merge_significant_data": json.dumps(merge_significant_data)
        }
        self.update_db_record("mann_test", {"_id": self.check_objectid(main_id)}, update_info)

if __name__ == "__main__":
    a = MannTest(None)
    a.add_mann_detail("5ec3492417b2bf4f5912d5ee",
                      "/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20210331/MannTest_rtmvhaqk6m453e0tij0bri2qv6_0331132843192767_6922/output/mann",
                      "/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20210331/MannTest_rtmvhaqk6m453e0tij0bri2qv6_0331132843192767_6922/output/mann/box_result_group.xls",
                      ['ESpr_A', 'ESpr_J', 'LSpr_A', 'LSpr_J', 'Sum_A', 'Sum_J', 'Aut_A', 'Aut_J', 'Win_A', 'Win_J'])

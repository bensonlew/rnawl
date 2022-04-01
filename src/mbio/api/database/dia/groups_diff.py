# -*- coding: utf-8 -*-


from biocluster.api.database.base import Base, report_check
from mbio.api.database.dia.api_base import ApiBase
import os
import json
from collections import defaultdict
import datetime
from bson.objectid import ObjectId
from types import StringTypes
import pandas as pd
import glob
import re


class GroupsDiff(ApiBase):
    def __init__(self, bind_object):
        super(GroupsDiff, self).__init__(bind_object)
        self._project_type = "dia"

    @report_check
    def add_groups_diff(self, name=None, params=None, metab_table_id=None):
        if not isinstance(metab_table_id, ObjectId):
            metab_table_id = ObjectId(metab_table_id)
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "运行结束",
            "name": name,
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "metab_table_id" : metab_table_id
        }
        # collection = self.db["sg_groups_diff"]
        # inserted_id = collection.insert_one(insert_data).inserted_id
        # collection.update_one({'_id': inserted_id}, {'$set': {'main_id': inserted_id}})
        inserted_id = self.create_db_table('sg_groups_diff', [insert_data])
        return inserted_id

    @report_check
    def add_groups_diff_detail(self, result_path, main_id):
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        data_list = []
        result_files = glob.glob(result_path + '/*_result.xls')
        boxfiles = glob.glob(result_path + '/*_boxfile.xls')
        result_file = result_files[0]
        box_file = boxfiles[0]
        rdata = pd.read_table(result_file, sep='\t',index_col=0)

        #获取group列表
        group_names = []
        for col in rdata.columns:
            mat_result = re.match('^(.*)-mean$', col)
            if mat_result:
                group_names.append(mat_result.group(1))
        ##group 设置编号
        group_map_id = {}
        id_map_group = {}
        for id,g in enumerate(group_names,1):
            group_map_id[g] = 'group' + str(id)
            id_map_group['group' + str(id)] = g
        #获取丰度排名
        sum_abund = rdata[[g+'-mean' for g in group_names]].apply(lambda x: sum(x),axis=1)
        rdata['rank'] = sum_abund.rank(ascending=False)

        bdata = pd.read_table(box_file, sep='\t',index_col=0)
        final_data = rdata.join(bdata)
        box_mat = 'min({}),Q1({}),Median({}),Q3({}),max({})'.split(',')

        #获取当前的可变列
        var_head_map = {"ID":'o_id'}
        cur_var = []
        var_head = []
        for h in var_head_map.keys():
            if h in  rdata.columns:
                cur_var.append(h)
                var_head.append(var_head_map[h])
        # final_data['Metabolite'].fillna(value='',inplace=True)
        for index in final_data.index:
            insert_data = {
                "diff_id" : main_id,
                "protein_id": index,
                # "metab" : final_data.loc[index]['Metabolite'] if final_data.loc[index]['Metabolite'] else index,
                "pvalue" : final_data.loc[index]['pvalue'],
                "fdr" : final_data.loc[index]['qvalue'],
                "abund_rank" : final_data.loc[index]['rank']
            }

            # 添加可变列的数据
            for c_var in cur_var:
                mongo_k = var_head_map[cur_var]
                insert_data[mongo_k] = final_data.loc[index][c_var]

            for g in group_names:
                g_id = group_map_id[g]
                g_mean = g+'-mean'
                g_sd = g+'-sd'
                #insert_data[g_id+'_mean'] = final_data.loc[index][g_mean]
                #insert_data[g_id+'_sd'] = final_data.loc[index][g_sd]
                g_id_mean_sd = ','.join([str(final_data.loc[index][g_mean]),str(final_data.loc[index][g_sd])])
                insert_data[g_id+'_mean'] = float(g_id_mean_sd.split(',')[0])
                insert_data[g_id+'_sd'] = float(g_id_mean_sd.split(',')[1])
                box_value = []
                for b in box_mat:
                    box_value.append(str(final_data.loc[index][b.format(g)]))
                box_value = ','.join(box_value)
                insert_data[g_id+'_box'] = box_value

            data_list.append(insert_data)

        # if table_type:
        #     for data in data_list:
        #         data['table_type'] = table_type

        try:
            # collection = self.db["sg_groups_diff_detail"]
            # collection.insert_many(data_list)
            self.create_db_table('sg_groups_diff_detail', data_list)
            # main_collection = self.db["sg_groups_diff"]
            # update_info = {"group_name": id_map_group, "main_id": main_id}
            if var_head:
                # update_info['var_head'] = ','.join(var_head)
                self.update_db_record('sg_groups_diff', main_id, var_head=var_head)
            # main_collection.update({"_id": main_id},{"$set": update_info})
            self.update_db_record('sg_groups_diff', main_id, group_name=id_map_group)
        except Exception as e:
            self.bind_object.logger.info("导入sg_groups_diff_detail数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入sg_groups_diff_detail数据成功")



        ##导入2组之间的pvlaue

        diff_group_map ={}
        for id in range(len(group_names)):
            g1 = group_names[id]
            for id2 in range(len(group_names)-1):
                g2 = group_names[id2]
                diff_group_map[g1+'-'+g2] = [group_map_id[g1], group_map_id[g2]]
                diff_group_map[g2+'-'+g1] = [group_map_id[g2], group_map_id[g1]]


        insert_pvalue = []

        for file in os.listdir(result_path):
            if '_result.xls' in file or '_boxfile.xls' in file:
                continue

            with open(result_path+'/'+file) as fr:
                diff = fr.readline().rstrip().split('\t')[4].rsplit('_pvalue',1)[0]
                name = '|'.join(sorted(diff_group_map[diff]))
                for line in fr:
                    spline = line.strip().split('\t')
                    protein_id = spline[0]
                    pvalue = spline[4]
                    eff_size = float(spline[1])
                    low_ci = float(spline[2])
                    up_ci = float(spline[3])

                    tmp = {
                        'protein_id' : protein_id,
                        'pvalue' : pvalue,
                        'diff_id' : main_id,
                        'name' : name,
                        'eff_size' :eff_size,
                        'low_ci' :low_ci,
                        'up_ci': up_ci
                    }
                    insert_pvalue.append(tmp)

        # if table_type:
        #     for data in insert_pvalue:
        #         data['table_type'] = table_type


        try:
            # collection = self.db["sg_groups_diff_relation"]
            # collection.insert_many(insert_pvalue)
            self.create_db_table('sg_groups_diff_relation', insert_pvalue)
        except Exception as e:
            self.bind_object.logger.info("导入sg_groups_diff_relation数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入sg_groups_diff_relation数据成功")

        self.update_db_record('sg_groups_diff', main_id, status="end",
                              main_id=main_id)

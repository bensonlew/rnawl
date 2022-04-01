# -*- coding: utf-8 -*-
from api_base import ApiBase
import os
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import glob
import json
import re
import pandas as pd
from api_base import ApiBase
import copy


class PlsdaRopls(ApiBase):
    def __init__(self, bind_object):
        super(PlsdaRopls, self).__init__(bind_object)



    def add_exp_diff_model(self, diff_id, pls_dir):
        diff_id = self.check_id(diff_id, "diff_id")
        if not os.path.exists(pls_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pls_dir))
        data_list = []


        pls_dict = {}
        intercept_files = glob.glob(pls_dir  + '/*intercept.xls')
        for eachfile in intercept_files:
            with open(eachfile, "r") as f:
                head = f.next()
                r2 = f.next().strip().split("\t")[1]
                q2 = f.next().strip().split("\t")[1]
                pls_dict["_r2"] = float(r2)
                pls_dict["_q2"] = float(q2)
        model_files = glob.glob(pls_dir + "/*model.xls")
        self.bind_object.logger.info(model_files)
        for eachfile in model_files:

            data_list = self.diff_model(diff_id, data_list, eachfile, pls_dict)
        try:
            collection = self.db['plsda_mode']
            collection.insert_many(data_list)

        except Exception as e:
            self.bind_object.set_error("导入表格plsda_mode信息出错:%s" , variables=(e))
        else:
            self.bind_object.logger.info("导入表格plsda_mode信息成功!")


    def diff_model(self, diff_id, data_list, eachfile, pls_dict):
        with open(eachfile, 'rb') as f:
            head = f.next()
            for line in f:
                line = line.strip().split('\t')
                a = line[0]
                r2x = float(line[1]) if line[1] != "NA" else "NA"
                r2x_cum = float(line[2])
                r2y = float(line[3]) if line[3] != "NA" else "NA"
                r2y_cum = float(line[4])
                q2 = float(line[5]) if line[5] != "NA" else "NA"
                q2_cum = float(line[6])
                insert_data = {
                    'plsda_id': diff_id,
                    'a': a,
                    'r2x': r2x,
                    'r2x_cum': r2x_cum,
                    'q2': q2,
                    'q2_cum': q2_cum,
                    'r2y_cum': r2y_cum,
                    'r2y': r2y,
                    'r2': pls_dict["_r2"],
                    'pq2': pls_dict["_q2"],
                }
                data_list.append(insert_data)
        return data_list

    def add_ci_ellipse(self,main_id, group_ci_ellipse, all_ci_ellipse):
        main_id = self.check_id(main_id, "main_id")
        insert_data = []
        with open(group_ci_ellipse) as fr, open(all_ci_ellipse) as fr2:
            fr.readline()
            fr2.readline()
            for id, file_o in enumerate([fr,fr2]):
                for line in file_o:
                    spline = line.strip().split('\t')
                    pc_x = 'pc' + re.search('(\d+)\D+(\d+)',spline[0]).group(1)
                    pc_y= 'pc' + re.search('(\d+)\D+(\d+)',spline[0]).group(2)
                    tmp = {
                        'name' : pc_x+'|'+pc_y,
                        'group' : spline[1],
                        #'data' : ','.join([spline[7],spline[6],spline[5],spline[2],spline[4],spline[3]]),
                        #'data': {'cx':float(spline[2]),'c11': float(spline[3]),'c12': float(spline[4]) ,'cy':float(spline[5]) ,'c21': float(spline[6]), 'c22':float(spline[7])},
                        'data' : ','.join([spline[2],spline[3],spline[4],spline[5],spline[6],spline[7]]), #cx,c11,c12,cy,c21,c22
                        'x': pc_x,
                        'y': pc_y,
                        'type' : 'ellipse',
                        'plsda_id' : main_id,
                        'method' : 'ci'
                    }
                    if  id == 0:
                        tmp['method_type'] = 'group'
                    else:
                        tmp['method_type'] = 'sample'
                    insert_data.append(tmp)
        try:
            collection = self.db['plsda_comp']
            collection.insert_many(insert_data)
        except Exception as e:
            self.bind_object.set_error("导入置信椭圆表格plsda_comp信息出错:%s" , variables=(e))
        else:
            self.bind_object.logger.info("导入置信椭圆表格plsda_comp信息成功!")

    def add_exp_diff_comp(self, diff_id, pls_dir, group_file):
        diff_id = self.check_id(diff_id, "diff_id")
        if not os.path.exists(pls_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pls_dir))
        if not os.path.exists(group_file):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(group_file))
        group_dict = {}
        with open(group_file, "r") as f:
            head = f.next()
            for line in f:
                line = line.strip().split('\t')
                sample = line[0]
                group = line[1]
                group_dict[sample] = group
        data_list = []

        comp_files = glob.glob(pls_dir  + '/*sites.xls')
        for eachfile in comp_files:
            with open(eachfile, 'rb') as f:
                head = f.next()
                for line in f:
                    line = line.strip().split('\t')
                    name = line[0]
                    group = group_dict[name]
                    pc1 = float(line[1])
                    pc2 = float(line[2])
                    insert_data = {
                        'plsda_id': diff_id,
                        'group': group,
                        'name': name,
                        'pc1': pc1,
                        'pc2': pc2,
                        'type': "scatter"
                    }
                    data_list.append(insert_data)
        ellipse_files = glob.glob(pls_dir  + '/*ellipse.xls')
        for eachfile in ellipse_files:
            with open(eachfile, 'rb') as f:
                head = f.next()
                for line in f:
                    line = line.strip().split('\t')
                    group = line[0]
                    m1 = line[1]
                    m2 = line[2]
                    c11 = line[3]
                    c12 = line[4]
                    c21 = line[5]
                    c22 = line[6]
                    #ellipse = [m1, m2, c11, c12, c21, c22]
                    #ellipse = ','.join([c22,c21,m2,m1,c12,c11])  #c22,c21,cy,cx,c12,c11
                    ellipse = ','.join( [m1, c11, c12, m2, c21, c22])
                    ## cx,c11,c12,cy,c21,c22
                    insert_data = {
                        'plsda_id': diff_id,
                        'group': group,
                        'name': "p1p2",
                        'x' : 'pc1',
                        'y' : 'pc2',
                        'type': "ellipse",
                        'method': 'trend',
                        'data': ellipse
                    }
                    if group == 'All':
                        insert_data['method_type'] = 'sample'
                    else:
                        insert_data['method_type'] = 'group'
                    data_list.append(insert_data)



        try:
            collection = self.db['plsda_comp']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格plsda_comp信息出错:%s" , variables=(e))
        else:
            self.bind_object.logger.info("导入表格plsda_comp信息成功!")


    def add_exp_diff_scatter(self, diff_id, pls_dir):
        diff_id = self.check_id(diff_id, "diff_id")
        if not os.path.exists(pls_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pls_dir), code="54702019")
        data_list = []

        scatter_files = glob.glob(pls_dir  + '/*permMN.xls')
        for eachfile in scatter_files:
            number_r = 0
            number_q = 0
            with open(eachfile, 'rb') as f:
                head = f.next()
                for line in f:
                    line = line.strip().split('\t')
                    name = line[0]
                    r2y = float(line[0])
                    q2 = float(line[1])
                    sim = float(line[2])
                    insert_data = {
                        'plsda_id': diff_id,
                        'x':sim,
                        'y':r2y,
                        'geom':"r2",
                        'type' : 'scatter'
                    }
                    data_list.append(insert_data)
                    if sim == 1 and number_r == 0:
                        number_r += 1
                        r_line_insert_data = {
                            'plsda_id': diff_id,
                            'x1':sim,
                            'y1':r2y,
                            "name" : 'R2',
                            'category':"R2",
                            'type':'divide_line'
                        }
                        #data_list.append(insert_data)


                    insert_data = {
                        'plsda_id': diff_id,
                        'x':sim,
                        'y':q2,
                        'geom':"q2",
                        'type' : 'scatter'
                    }
                    data_list.append(insert_data)
                    if sim == 1 and number_q == 0:
                        number_q += 1
                        q_line_insert_data = {
                            'plsda_id': diff_id,
                            'x1':sim,
                            'y1':q2,
                            'name' : 'Q2',
                            'category':"Q2",
                            'type':"divide_line"
                        }
                        #data_list.append(insert_data)
        intercept_files = glob.glob(pls_dir  + '/*intercept.xls')
        for eachfile in intercept_files:
            with open(eachfile, "r") as f:
                head = f.next()
                r2 = f.next().strip().split("\t")[1]
                q2 = f.next().strip().split("\t")[1]
                r_line_insert_data.update({
                    'x2': 0,
                    'y2': float(r2)
                })
                data_list.append(r_line_insert_data)
                q_line_insert_data.update({
                    'x2': 0,
                    'y2': float(q2)
                })
                data_list.append(q_line_insert_data)
        try:
            collection = self.db['plsda_scatter']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格plsda_scatter信息出错:%s" , variables=(e), code="54702021")
        else:
            self.bind_object.logger.info("导入表格plsda_scatter信息成功!")

    def add_exp_diff_load(self, diff_id, pls_dir):
        diff_id = self.check_id(diff_id, "diff_id")
        data_list=[]
        load_file = glob.glob(pls_dir  + '/*loading.xls')[0]
        vip_file = glob.glob(pls_dir+'/*.vip.xls')[0]
        metab_name_map_vip = self.get_metab_vip(vip_file)
        with open(load_file, 'rb') as f:
            head = f.next()
            for line in f:
                line = line.strip().split('\t')
                name = line[0]
                insert_data = {
                    "name" : name,
                    "comp1" : float(line[1]),
                    "comp2" : float(line[2]),
                    "plsda_id" :diff_id,
                    "vip" : metab_name_map_vip[name]
                }
                data_list.append(insert_data)
        try:
            collection = self.db['plsda_var']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格plsda_var信息出错:%s" , variables=(e))
        else:
            self.bind_object.logger.info("导入表格plsda_var信息成功!")

        ##update main
        self.update_main_collection(diff_id)  ##20200507

    def get_metab_vip(self,file):
        retd = {}
        with open(file) as f:
            f.readline()
            for line in f:
                line = line.strip()
                spline = line.split('\t')
                retd[spline[0]] = spline[1]
        return retd

    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type))
        return object_id

    def update_main_collection(self, main_id):
        main_id = self.check_id(main_id, 'ObjectId')
        table_data_colum = [
            {"field": "name", "type": "string", "sort": "false", "title": ""},
            {"field": "pc1", "type": "string", "sort": "false", "title": "comp1"},
            {"field": "pc2", "type": "string", "sort": "false", "title": "comp2"}
        ]
        table_data = {'column': table_data_colum,'condition':{"type":"scatter"}}


        table_data2_column = [
            {"field": "name", "type": "string", "sort": "false", "title": ""},
            {"field": "comp1", "type": "string", "sort": "false", "title": "comp1"},
            {"field": "comp2", "type": "string", "sort": "false", "title": "comp2"},
            #{"field": "comp3", "type": "string", "sort": "false", "title": "comp3"}
            ]
        table_data2 = {'column': table_data2_column,'condition':{}}

        table_data3_column = [
            {"field": "a", "type": "string", "sort": "false", "title": "Name"},
            {"field": "r2x", "type": "string", "sort": "false", "title": "R2X"},
            {"field": "r2x_cum", "type": "string", "sort": "false", "title": "R2X(cum)"},
            {"field": "r2y", "type": "string", "sort": "false", "title": "R2Y"},
            {"field": "r2y_cum", "type": "string", "sort": "false", "title": "R2Y(cum)"},
            {"field": "q2", "type": "string", "sort": "false", "title": "Q2"},
            {"field": "q2_cum", "type": "string", "sort": "false", "title": "Q2(cum)"}]
        table_data3 = {'column': table_data3_column,'condition':{}}

        table_data4_column = [
            {"field": "name", "type": "string", "sort": "false", "title": ""},
            {"field": "vip", "type": "string", "sort": "false", "title": "VIP"}]

        table_data4 = {'column': table_data4_column,'condition':{}}

        ellipse_data = {"name":"name","data_option":["ci",'trend'],"data":["data","data"],"condition":{"type":"ellipse"}}
        #scatter_data2 = {"name":"","data":["x","y"], "category":"geom","condition":{"geom":["r2","q2"]}}
        scatter_data2 = {"name":"name","data":["x","y"], "category":"geom","condition":{"type":'scatter'}}
        scatter_data = {"name":"name","data":["pc1","pc2"],"category":"group","condition":{"type":"scatter"}}
        #line_data = {"name":"","condition":{"geom":["lr","lq"]}}
        #lr_line_data = {"name":"lr","condition":{"type": "r_line"}}
        #lq_line_data = {"name":"lq","condition":{"type": "q_line"}}
        divide_line_data = {"name":"name","condition":{"type": "divide_line"}}

        update_info = {
            'main_id' : main_id,
            'sample_table_data' : json.dumps(table_data),
            'varload_table_data' : json.dumps(table_data2),
            'model_table_data' : json.dumps(table_data3),
            'vip_table_data' : json.dumps(table_data4),
            'ellipse_data' : json.dumps(ellipse_data),
            'model_scatter_data' : json.dumps(scatter_data2),
            'pc_scatter_data' : json.dumps(scatter_data),
            'divide_line_data' : json.dumps(divide_line_data)

        }
        self.db['plsda'].update({"_id": main_id}, {"$set": update_info})
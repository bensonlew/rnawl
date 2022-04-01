# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from collections import OrderedDict
from bson import SON
from biocluster.config import Config
import types
import re
import datetime
import os
import json
from api_base import ApiBase

class Dbrda(ApiBase):
    def __init__(self, bind_object):
        super(Dbrda, self).__init__(bind_object)
        self._db_name = "sanger_tool_lab"
        self._tables = []  # 记录存入了哪些表格

    def add_detail(self,dir_path,main_id,sn_group,feature_index):
        # _main_collection = self.db['rda']
        rda_files = os.listdir(dir_path)
        self.feature_index = feature_index
        if 'db_rda_sites.xls' in rda_files:
            method = 'db_rda'
        else:
            self.set_error('DbRDA分析没有生成正确的结果数据')
        feature_insert = [{
            "field" : "scatter",
            "title" :"scatter"
        }]
        update_dict5 = {
            "feature": json.dumps(feature_insert)
        }
        self.update_db_record("dbrda",{'_id': self.check_objectid(main_id)}, update_dict5)
        site_path = dir_path.rstrip('/') + '/{}_sites.xls'.format(method)
        self.insert_table_detail(site_path, 'scatter', coll_name= "dbrda_scatter",update_id=main_id,group = sn_group)
        new_scatter_data = open(site_path,"r").readline().rstrip().split('\t')[1:]
        scatter_data_insert = {
            "name": "name",
            "data": new_scatter_data,
            "category": "group",
            "series":"series",
            "condition": {'type': "scatter"},
        }
        update_dict = {
            "scatter_data": json.dumps(scatter_data_insert)
        }
        self.update_db_record("dbrda",{'_id': self.check_objectid(main_id)}, update_dict)
        filelist = os.listdir(dir_path)
        if '{}_biplot.xls'.format(method) in filelist:
            env_vec_path = dir_path.rstrip('/') + '/{}_biplot.xls'.format(method)
            self.insert_arrow_detail(env_vec_path, update_id=self.check_objectid(main_id),coll_name= "dbrda_arrow",group='env')
            arrow_data_insert = {
                "name": "name",
                # "category": "group",
                "condition": {'type': "arrow"},
            }
            update_dict1 = {
                "arrow_data": json.dumps(arrow_data_insert)
            }
            self.update_db_record("dbrda",{'_id': self.check_objectid(main_id)}, update_dict1)
        if '{}_importance.xls'.format(method) in filelist: 
            importance_path = dir_path.rstrip('/') + '/{}_importance.xls'.format(method)
            # self.insert_table_detail(importance_path, 'importance', update_id=main_id)
            update_importances = []
            with open(importance_path,'r') as ip:
                line = ip.readline()
                fd = line.rstrip().split('\t')
                ind1 = fd.index('CAP1')
                ind2 = fd.index('CAP2')
                while 1:
                    tem_importance = {}
                    tem_importance1 = {}
                    line = ip.readline()
                    if not line:
                        break
                    fd = line.rstrip('\r\n').split('\t')
                    if fd[0] != "Proportion Explained":
                        continue
                    tem_importance['field'] = round(float(fd[ind1])*100,2)
                    tem_importance['title'] = "CAP1"
                    update_importances.append(tem_importance)
                    tem_importance1['field'] = round(float(fd[ind2])*100,2)
                    tem_importance1['title'] = "CAP2"
                    update_importances.append(tem_importance1)
            update_dict3 = {"groups":json.dumps(update_importances)}
            self.update_db_record("dbrda",{'_id': self.check_objectid(main_id)}, update_dict3)
        plot_feature_path = dir_path.rstrip('/') + '/{}_plot_feature_data.xls'.format(method)
        self.insert_table_detail(plot_feature_path, 'scatter', coll_name= "dbrda_scatter",update_id=self.check_objectid(main_id))
        ellipse_path = dir_path.rstrip('/') + '/ellipse.xls'
        self.insert_ellipse_table1(ellipse_path,self.check_objectid(main_id))
        ellipse_data_insert = {
            "name": "name",
            "group": "group",
            "data":"data",
            "condition": {'type': "ellipse"},
            "data_option":["confidence","grouped"],
        }
        update_dict2 = {
            "ellipse_data": json.dumps(ellipse_data_insert)
        }
        self.update_db_record("dbrda",{'_id': self.check_objectid(main_id)}, update_dict2)
        self.bind_object.logger.info('db_RDA分析结果导入数据库完成.')
    
    def insert_table_detail(self, file_path, table_type, update_id,
                            coll_name,
                            main_coll='dbrda', group=None,
                            update_column=True, db=None, fileter_biplot=None, remove_key_blank=False,
                            split_fullname=False, colls=None):
        # self._tables.append(table_type)
        if not db:
            db = self.db
        collection = db[coll_name]
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            columns = all_lines[0].rstrip().split('\t')[1:]
            if remove_key_blank:
                columns = [i.replace(' ', '') for i in columns]
            if colls:
                columns = colls
            data_temp = []
            num = 0
            for line in all_lines[1:]:
                values = line.rstrip().split('\t')
                if fileter_biplot:
                    if not isinstance(fileter_biplot, list):
                        self.bind_object.set_error('需要删除的匹配列必须为列表', code="51000308")
                    flag = 0
                    for i in fileter_biplot:
                        if re.match(r'^{}'.format(i), values[0]):
                            flag = 1
                    if flag:
                        continue
                else:
                    pass
                insert_data = {
                    'dbrda_id': self.check_objectid(update_id),
                    'type': table_type,
                }
                if group:
                    if split_fullname:
                        insert_data['fullname'] = values[0]
                        insert_data['name'] = values[0].split(';')[-1].strip()
                        insert_data['group'] = group[values[0].split(';')[-1].strip()]
                        insert_data['series'] = 'sample'
                        insert_data['num'] = num
                        num += 1
                    else:
                        insert_data['name'] = values[0]
                        insert_data['group'] = group[values[0]]
                        insert_data['series'] = 'sample'
                        insert_data['num'] = num
                        num += 1
                else:
                    if split_fullname:
                        insert_data['fullname'] = values[0]
                        insert_data['name'] = values[0].split(';')[-1].strip()
                        insert_data['group'] = 'feature'
                        insert_data['series'] = 'feature'
                        insert_data['num'] = self.feature_index[values[0]]
                        # insert_data['group'] = group[values[0].split(';')[-1].strip()]
                    else:
                        insert_data['name'] = values[0]
                        insert_data['group'] = 'feature'
                        insert_data['series'] = 'feature'
                        insert_data['num'] = self.feature_index[values[0]]
                        # insert_data['group'] = group[values[0]]
                values_dict = dict(zip(columns, [float(x) for x in values[1:]]))
                data_temp.append(dict(insert_data, **values_dict))
            if data_temp:
                collection.insert_many(data_temp)
            else:
                return None
            
    def insert_arrow_detail(self, file_path, update_id,
                            coll_name,group=None,
                            update_column=True, db=None, fileter_biplot=None, remove_key_blank=False,
                            split_fullname=False, colls=None):
        if not db:
            db = self.db
        collection = db[coll_name]
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            columns = all_lines[0].rstrip().split('\t')[1:]
            data_temp = []
            num = 0
            for line in all_lines[1:]:
                values = line.rstrip().split('\t')[1:]
                if len(columns)%2 == 1:
                    for i in range(len(columns)):
                        if i == len(columns):
                            insert_data = {
                                'dbrda_id': self.check_objectid(update_id),
                                'type': 'arrow',
                            }
                            insert_data['x'] = columns[i]
                            insert_data['y'] = columns[i+1]
                            insert_data['x1'] = 0
                            insert_data['y1'] = 0
                            insert_data['x2'] = float(values[i])
                            insert_data['y2'] = float(values[i+1])
                            insert_data['name'] = line.rstrip().split('\t')[0]
                            insert_data['category'] = group
                            insert_data['num'] = num
                            data_temp.append(insert_data)
                            break
                        if i%2 == 1:
                            continue
                        insert_data = {
                            'dbrda_id': update_id,
                            'type': 'arrow',
                        }
                        insert_data['x'] = columns[i]
                        insert_data['y'] = columns[i+1]
                        insert_data['x1'] = 0
                        insert_data['y1'] = 0
                        insert_data['x2'] = float(values[i])
                        insert_data['y2'] = float(values[i+1])
                        insert_data['name'] = line.rstrip().split('\t')[0]
                        insert_data['category'] = group
                        insert_data['num'] = num
                        data_temp.append(insert_data)
                else:
                    for i in range(len(columns)):
                        if i%2 == 1:
                            continue
                        insert_data = {
                            'dbrda_id': update_id,
                            'type': 'arrow',
                        }
                        insert_data['x'] = columns[i]
                        insert_data['y'] = columns[i+1]
                        insert_data['x1'] = 0
                        insert_data['y1'] = 0
                        insert_data['x2'] = float(values[i])
                        insert_data['y2'] = float(values[i+1])
                        insert_data['name'] = line.rstrip().split('\t')[0]
                        insert_data['category'] = group
                        insert_data['num'] = num
                        data_temp.append(insert_data)
                num += 1
            if data_temp:
                collection.insert_many(data_temp)
            else:
                return None

    # def insert_main_tables(self, tables, update_id, main='dbrda'):
    #     """
    #     """
    #     main_collection = self.db[main]
    #     main_collection.update_one({'_id': update_id},
    #                                {'$set': {'tables': ','.join(tables)}},
    #                                upsert=False) 

    def insert_ellipse_table1(self, infile, main_id):
        insert_data = []
        with open(infile,"r") as f:
            group = f.readline().rstrip("\r\n").split('\t')[1:]
            name = ''
            tmp = {}
            for line in f:
                line = line.strip()
                spline = line.split('\t')
                for value in spline[1:]:

                    name = group[spline[1:].index(value)]
                    tmp = {'name': name,"group": name ,'type': 'ellipse', 'dbrda_id': self.check_objectid(main_id),"method":"grouped"}
                    tmp["data"] = value
                    numbers = re.match(r"(.*\d)(.*\d)",spline[0])
                    tmp['x'] = numbers.group(1)
                    tmp['y'] = numbers.group(2)
                    insert_data.append(tmp)
        try:
            collection = self.db['dbrda_ellipse_detail']
            collection.insert_many(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入dbrda置信椭圆表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入dbrda置信椭圆表格成功")

    def add_ellipse(self, infile, main_id):
        insert_data = []
        with open(infile,"r") as f:
            f.readline()
            name = ''
            tmp = {}
            for line in f:
                line = line.strip()
                spline = line.split('\t')
                name = spline[1]
                tmp = {'name': name,"group": name, 'type': 'ellipse', 'dbrda_id': self.check_objectid(main_id),"method":"confidence"}
                k = spline[1]
                tmp["data"] = ','.join(spline[2:])
                numbers = re.match(r"(.*\d)[_|-](.*\d)",spline[0])
                tmp['x'] = numbers.group(1)
                tmp['y'] = numbers.group(2)
                insert_data.append(tmp)
        try:
            collection = self.db['dbrda_ellipse_detail']
            collection.insert_many(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入dbrda置信椭圆表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入dbrda置信椭圆表格成功")
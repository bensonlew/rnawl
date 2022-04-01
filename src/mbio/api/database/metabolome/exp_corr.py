# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20180522
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
#import json


class ExpCorr(Base):
    def __init__(self, bind_object):
        super(ExpCorr, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_exp_corr(self, specimen=None, name=None, main_id=None, params =None, tree_file=None,list_file=None,table_type=None):
        #metab_table_id = self.check_id(metab_table_id, 'metab_table_id')
        if not main_id:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'ExpCorr_Origin',
                'params': params if params else '',
                'status': 'end',
                #'metab_table_id': metab_table_id,
                'specimen': specimen,
                'tree': '',
                'main_id': ''
            }
            if table_type:
                insert_data['table_type'] = table_type
            try:
                collection = self.db['exp_corr']
                exp_corr_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", exp_corr_id, exp_corr_id)
            except Exception, e:
                self.bind_object.set_error('导入exp_corr主表异常:%s', variables=(e), code="54700501")
        else:
            self.update_table("main_id", main_id, main_id)
            exp_corr_id = main_id
        if tree_file:
            if not os.path.exists(tree_file):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(tree_file), code="54700502")
            with open(tree_file,"r") as f:
                specimen_tree = f.readline().strip()
                specimen =  f.readline().strip().split(";")
            self.update_table("tree", specimen_tree, exp_corr_id)
            self.update_table("specimen", specimen, exp_corr_id)
        if list_file:
            if not os.path.exists(list_file):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(list_file), code="54700503")
            with open(list_file,"r") as f2:
                specimen = f2.readline().strip().split("\t")
                specimen = specimen[0:len(specimen)]
            self.update_table("specimen", specimen, exp_corr_id)
        return exp_corr_id

    @report_check
    def add_exp_corr_detail(self, exp_corr_id, corr_file,p_file):
        exp_corr_id = self.check_id(exp_corr_id,"exp_corr_id")
        if not os.path.exists(corr_file):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(corr_file), code="54700504")
        if not os.path.exists(p_file):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(p_file), code="54700505")
        data_list = []
        result = self.db['exp_corr'].find_one({'_id': exp_corr_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(exp_corr_id), code="54700506")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")
        with open(corr_file, 'rb') as f:
            head = f.next()
            sams = head.strip().split("\t")
            for line in f:
                line = line.strip().split('\t')
                name = line[0]
                insert_data = {
                    'corr_id': exp_corr_id,
                    'name':name,
                    'type': 'corr'
                }
                for i in range(0, len(sams)):
                    sam_corr = float(line[i + 1])
                    '''
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                    '''
                    insert_data[sams[i]] = sam_corr
                data_list.append(insert_data)
            try:
                collection = self.db['exp_corr_detail']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(corr_file, e), code="54700507")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % corr_file)
        data_list = []
        with open(p_file, 'rb') as f:
            head = f.next()
            sams = head.strip().split("\t")
            for line in f:
                line = line.strip().split('\t')
                name = line[0]
                insert_data = {
                    'corr_id': exp_corr_id,
                    'name':name,
                    'type': 'pvalue'
                }
                for i in range(0, len(sams)):
                    sam_corr = float(line[i + 1])
                    '''
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                    '''
                    insert_data[sams[i]] = sam_corr
                data_list.append(insert_data)
            try:
                collection = self.db['exp_corr_detail']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(p_file, e), code="54700508")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % p_file)

    @report_check
    def update_table(self, str, name, main_table_id):
        try:
            self.db['exp_corr'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {str: name}})
        except Exception as e:
            self.bind_object.set_error('更新exp_corr主表%s字段出错:%s', variables=(str,e), code="54700509")

    @report_check
    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type), code="54700510")
        return object_id

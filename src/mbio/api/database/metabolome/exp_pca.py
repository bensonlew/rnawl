# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20180620
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId


#import json


class ExpPca(Base):
    def __init__(self, bind_object):
        super(ExpPca, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_exp_pca(self, name=None, main_id=None, params =None,table_type=None):
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
                'name': name if name else 'ExpPca_Origin',
                'params': params if params else '',
                'status': 'end',
                'main_id': ''
            }
            if table_type:
                insert_data['table_type'] = table_type
            try:
                collection = self.db['exp_pca']
                exp_pca_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", exp_pca_id, exp_pca_id)
            except Exception, e:
                self.bind_object.set_error('导入exp_pca主表异常:%s', variables=(e), code="54701701")
        else:
            self.update_table("main_id", main_id, main_id)
            exp_pca_id = main_id
        return exp_pca_id

    @report_check
    def add_exp_pca_detail(self, exp_pca_id, pca_file, ellipse_file, group_file=None):
        exp_pca_id = self.check_id(exp_pca_id, "exp_pca_id")
        if not os.path.exists(pca_file):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pca_file), code="54701702")
        if not os.path.exists(ellipse_file):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(ellipse_file), code="54701703")
        data_list = []
        result = self.db['exp_pca'].find_one({'_id': exp_pca_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(exp_pca_id), code="54701704")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")
        if group_file:
            group_dict = {}
            with open(group_file, 'rb') as gf:
                head = gf.next()
                for line in gf:
                    line = line.strip().split('\t')
                    s_name = line[0]
                    g_name = line[1]
                    group_dict[s_name] = g_name
        with open(pca_file, 'rb') as f:
            head = f.next()
            for line in f:
                line = line.strip().split('\t')
                name = line[0]
                p1 = line[1]
                p2 = line[2]
                group = group_dict[name]
                insert_data = {
                    'pca_id': exp_pca_id,
                    'name': name,
                    'type': 'specimen',
                    'group': group,
                    'pc1': p1,
                    'pc2': p2,
                }
                data_list.append(insert_data)
        with open(ellipse_file, 'rb') as ef:
            head = ef.next()
            for line in ef:
                line = line.strip().split('\t')
                group = line[0]
                m1 = line[1]
                m2 = line[2]
                c11 = line[3]
                c12 = line[4]
                c21 = line[5]
                c22 = line[6]
                ellipse = [m1, m2, c11, c12, c21, c22]
                insert_data = {
                    'pca_id': exp_pca_id,
                    'name': "p1p2",
                    'type': 'circle',
                    'group': group,
                    'ellipse': ellipse,
                }
                data_list.append(insert_data)
            try:
                collection = self.db['exp_pca_detail']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(pca_file, e), code="54701705")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % pca_file)

    @report_check
    def add_exp_pca_model(self, exp_pca_id, model_file):
        exp_pca_id = self.check_id(exp_pca_id, "exp_pca_id")
        if not os.path.exists(model_file):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(model_file), code="54701706")
        data_list = []
        result = self.db['exp_pca'].find_one({'_id': exp_pca_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(exp_pca_id), code="54701707")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")
        with open(model_file, 'rb') as f:
            head = f.next()
            for line in f:
                line = line.strip().split('\t')
                pc = line[0]
                r2x = line[1]
                r2x_cum = line[2]
                insert_data = {
                    'pca_id': exp_pca_id,
                    'a': pc,
                    'r2x': r2x,
                    'r2x_cum': r2x_cum
                }
                data_list.append(insert_data)
            try:
                collection = self.db['exp_pca_model']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(model_file, e), code="54701708")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % model_file)

    @report_check
    def update_table(self, str, name, main_table_id):
        try:
            self.db['exp_pca'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {str: name}})
        except Exception as e:
            self.bind_object.set_error('更新exp_pca主表%s字段出错:%s', variables=(str, e), code="54701709")

    @report_check
    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type), code="54701710")
        return object_id

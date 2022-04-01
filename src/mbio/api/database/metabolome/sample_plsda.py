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


class SamplePlsda(Base):
    def __init__(self, bind_object):
        super(SamplePlsda, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_sample_plsda(self, metab_table_id=None, name=None, main_id=None, params=None):
        if not main_id:
            metab_table_id = self.check_id(metab_table_id, "metab_table_id")
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'Plsda_Origin',
                'params': params if params else '',
                'status': 'end',
                'main_id': '',
                "version" : "3.0",
                "metab_table" :metab_table_id
            }

            try:
                collection = self.db['sample_plsda']
                diff_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", diff_id, diff_id)
            except Exception as e:
                self.bind_object.set_error('导入主表异常:%s', variables=(e))
        else:
            self.update_table("main_id", main_id, main_id)
            diff_id = main_id
        return diff_id

    @report_check
    def add_exp_diff_model(self, diff_id, pls_dir):
        diff_id = self.check_id(diff_id, "diff_id")
        if not os.path.exists(pls_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pls_dir))
        data_list = []
        result = self.db['sample_plsda'].find_one({'_id': diff_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(diff_id), code="54702008")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")

        pls_dict = {}
        intercept_files = glob.glob(pls_dir  + '/*intercept.xls')
        for eachfile in intercept_files:
            model = 'plsda'
            with open(eachfile, "r") as f:
                head = f.next()
                r2 = f.next().strip().split("\t")[1]
                q2 = f.next().strip().split("\t")[1]
                pls_dict[model + "_r2"] = float(r2)
                pls_dict[model + "_q2"] = float(q2)
        model_files = glob.glob(pls_dir + "/*model.xls")
        self.bind_object.logger.info(model_files)
        for eachfile in model_files:
            model = 'plsda'
            data_list = self.diff_model(diff_id, data_list, eachfile, pls_dict)
        try:
            collection = self.db['sample_plsda_mode']
            if len(data_list) > 0:
                collection.insert_many(data_list)

        except Exception as e:
            self.bind_object.set_error("导入表格sample_plsda_mode信息出错:%s" , variables=(e), code="54702009")
        else:
            self.bind_object.logger.info("导入表格sample_plsda_mode信息成功!")

    @report_check
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
                    'r2': pls_dict["plsda_r2"],
                    'pq2': pls_dict["plsda_q2"],
                }
                data_list.append(insert_data)
        return data_list

    @report_check
    def add_exp_diff_comp(self, diff_id, pls_dir, group_file):
        diff_id = self.check_id(diff_id, "diff_id")
        if not os.path.exists(pls_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pls_dir), code="54702010")
        if not os.path.exists(group_file):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(group_file), code="54702011")
        group_dict = {}
        with open(group_file, "r") as f:
            head = f.next()
            for line in f:
                line = line.strip().split('\t')
                sample = line[0]
                group = line[1]
                group_dict[sample] = group
        data_list = []
        result = self.db['sample_plsda'].find_one({'_id': diff_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(diff_id), code="54702012")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")
        #groups = os.listdir(pls_dir)
        #for eachgroup in groups:
        comp_files = glob.glob(pls_dir  + '/*sites.xls')
        for eachfile in comp_files:
            #model = "plsda"
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
                        'type': "specimen"
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
                    ellipse = [m1, m2, c11, c12, c21, c22]
                    insert_data = {
                        'plsda_id': diff_id,
                        'group': group,
                        'name': "p1p2",
                        'type': "circle",
                        'ellipse': ellipse
                    }
                    data_list.append(insert_data)
        try:
            collection = self.db['sample_plsda_comp']
            if len(data_list) > 0:
                collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格sample_plsda_comp信息出错:%s" , variables=(e), code="54702013")
        else:
            self.bind_object.logger.info("导入表格sample_plsda_comp信息成功!")

    @report_check
    def add_exp_diff_bar(self, diff_id, pls_dir):
        diff_id = self.check_id(diff_id, "diff_id")
        if not os.path.exists(pls_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pls_dir), code="54702014")
        data_list = []
        result = self.db['sample_plsda'].find_one({'_id': diff_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(diff_id), code="54702015")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")
        # groups = os.listdir(pls_dir)
        # for eachgroup in groups:
        comp_files = glob.glob(pls_dir  + '/*PLS-DA.model.xls')
        for eachfile in comp_files:
            number = 0
            with open(eachfile, 'rb') as f:
                head = f.next()
                for line in f:
                    line = line.strip().split('\t')
                    number += 1
                    comp = "comp" + str(number)
                    r2y_cum = float(line[4])
                    q2_cum = float(line[6])
                    if float(line[4]) <0 or float(line[6]) <0:
                        self.bind_object.logger.info("贡献度出现负值，建议重新选择参数运行！")
                    insert_data = {
                        'plsda_id': diff_id,
                        'x': comp,
                        'q2': q2_cum,
                        'r2y': r2y_cum
                    }
                    data_list.append(insert_data)

        try:
            collection = self.db['sample_plsda_bar']
            if len(data_list) > 0:
                collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格sample_plsda_bar信息出错:%s" , variables=(e), code="54702018")
        else:
            self.bind_object.logger.info("导入表格sample_plsda_bar信息成功!")

    @report_check
    def add_exp_diff_scatter(self, diff_id, pls_dir):
        diff_id = self.check_id(diff_id, "diff_id")
        if not os.path.exists(pls_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pls_dir), code="54702019")
        data_list = []
        #groups = os.listdir(pls_dir)
        #mytest = []
        #for eachgroup in groups:
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
                        'geom':"r2"
                    }
                    data_list.append(insert_data)
                    if sim == 1 and number_r == 0:
                        number_r += 1
                        insert_data = {
                            'plsda_id': diff_id,
                            'x':sim,
                            'y':r2y,
                            'geom':"lr"
                        }
                        data_list.append(insert_data)
                    insert_data = {
                        'plsda_id': diff_id,
                        'x':sim,
                        'y':q2,
                        'geom':"q2"
                    }
                    data_list.append(insert_data)
                    if sim == 1 and number_q == 0:
                        number_q += 1
                        insert_data = {
                            'plsda_id': diff_id,
                            'x':sim,
                            'y':q2,
                            'geom':"lq"
                        }
                        data_list.append(insert_data)
        intercept_files = glob.glob(pls_dir  + '/*intercept.xls')
        for eachfile in intercept_files:
            with open(eachfile, "r") as f:
                head = f.next()
                r2 = f.next().strip().split("\t")[1]
                q2 = f.next().strip().split("\t")[1]
                insert_data = {
                    'plsda_id': diff_id,
                    'x': 0,
                    'y': float(r2),
                    'geom':"lr"
                }
                data_list.append(insert_data)
                insert_data = {
                    'plsda_id': diff_id,
                    'x': 0,
                    'y': float(q2),
                    'geom':"lq"
                }
                data_list.append(insert_data)
        try:
            collection = self.db['sample_plsda_scatter']
            if len(data_list) > 0:
                collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格sample_plsda_scatter信息出错:%s" , variables=(e), code="54702021")
        else:
            self.bind_object.logger.info("导入表格sample_plsda_scatter信息成功!")

    @report_check
    def update_table(self, str, name, main_table_id):
        try:
            self.db['sample_plsda'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {str: name}})
        except Exception as e:
            self.bind_object.set_error('更新sample_plsda主表%s字段出错:%s', variables=(str, e), code="54702022")

    @report_check
    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type), code="54702023")
        return object_id

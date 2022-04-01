# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
# last_modify:20180521
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId
from contextlib import nested
from mbio.packages.metabolome.common import check_metab
from mbio.packages.metabolome.common import check_info
import pandas as pd
import numpy as np
import re


class Preprocess(Base):
    # 对应表格链接：http://git.majorbio.com/liu.linmeng/metabolome/wikis/collection/preprocess/metab_table
    def __init__(self, bind_object):
        super(Preprocess, self).__init__(bind_object)
        self._project_type = "metabolome"
        self.sample_list = []

    @report_check
    # def add_metab_table(self, name, type, rere_path, group_path=None, params=None):
    def add_metab_table(self, name, type, rere_path, is_raw=0, params=None):
        """
        工作流对metab_table主表进行导表
        :param name: 工作流需要导两张主表，第一张对应Raw_，为原始数据表；第二张对应Table_，为预处理后的结果
        :param type: 流程的类型, {GC|LC}
        :param params: 工作流选择的运行参数
        :return:
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': '代谢表',
            'name': name,
            'created_ts': created_ts,
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end',
            'type': type,
            'setted': 1,
            'table_path': rere_path,
            'is_raw': is_raw,
            'version' : '3.0'
        }
        # if group_path is None:
        #     insert_data['is_raw'] = 0
        # else:
        #     if not os.path.isfile(group_path):
        #         raise Exception("分组文件路径错误%s" % group_path)
        #     insert_data['is_raw'] = 1
        #     insert_data['group_path'] = group_path
        collection = self.db['exp']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({'_id': main_id}, {'$set': {'main_id': main_id}})
        return main_id

    @report_check
    def add_metab_table_detail(self, main_id, table_path, raw_exp_id=''):
        """
        对阴离子，阳离子和合成表分别进行导表
        :param main_id: 主表id
        :param table_path: 三种类型数据表各自的路径，里面对应存储代谢物明细表和丰度表，逗号分割
        :return:
        """
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54702201")
        # main_col = self.db['metab_table']
        main_col = self.db['exp']
        main_result = main_col.find_one({"_id": main_id})
        table_list = table_path.split(',')
        for one_table in table_list:
            if not os.path.isdir(one_table):
                self.bind_object.set_error("代谢表结果%s,不存在此目录", variables=(one_table), code="54702202")
        if main_result['is_raw'] == 0 and main_result['type'] == "LC":
            if len(table_list) != 3:
                self.bind_object.set_error('预处理LC数据必须提供三种类型的数据结果，现在提供个数：%s', variables=(len(table_list)),
                                           code="54702203")
            type_list = ['pos', 'neg', 'mix']
            self._add_metab_table_detail(main_id, table_list[0], type_list[0])
            self._add_metab_table_detail(main_id, table_list[1], type_list[1])
            self._add_metab_table_detail(main_id, table_list[2], type_list[2])
        elif main_result['is_raw'] == 1 and main_result['type'] == "LC":
            if len(table_list) != 2:
                self.bind_object.set_error('LC数据的原始结果表必须提供两种类型的数据结果，现在提供个数： %s', variables=(len(table_list)),
                                           code="54702204")
            type_list = ['pos', 'neg']
            self._add_metab_table_detail(main_id, table_list[0], type_list[0])
            self._add_metab_table_detail(main_id, table_list[1], type_list[1])
            ##冗余mix表 20210331
            self._add_metab_table_detail(main_id, table_list[0], 'mix')
            self._add_metab_table_detail(main_id, table_list[1], 'mix')

        elif main_result['type'] == 'GC':
            if len(table_list) != 1:
                self.bind_object.set_error("GC数据的结果表只应包含阳离子表一种类型结果，现在提供个数：%s", variables=(len(talbe_list)),
                                           code="54702205")
            self._add_metab_table_detail(main_id, table_list[0], 'pos')
        # main_col.update_one({'_id': main_id}, {'$set': {'table_path': table_path, 'specimen': self.sample_list}})

        update_tmp_dic = {'specimen': self.sample_list}
        ##add exp_stat and exp_cv v3.0 20200306
        tar_dir = os.path.abspath(os.path.join(table_list[0],'..'))
        ##v1 没有
        if raw_exp_id:
            #工作流用，交互分析用
            if str(main_id) != str(raw_exp_id):
                if os.path.exists(tar_dir+'/statistic.xls'):
                    self._add_exp_stat(tar_dir+'/statistic.xls', main_id)

                if os.path.exists(tar_dir + '/cv_percent.xls'):
                    self._add_exp_cv(tar_dir+'/cv_percent.xls',  main_id)
            #工作流用
            else:
                if os.path.exists(tar_dir+'/raw_statistic.xls'):
                    self._add_exp_stat(tar_dir+'/raw_statistic.xls', main_id)

                if os.path.exists(tar_dir + '/raw_cv_percent.xls'):
                    self._add_exp_cv(tar_dir + '/raw_cv_percent.xls',  main_id)

            if not isinstance(raw_exp_id, ObjectId):
                raw_exp_id = ObjectId(raw_exp_id)
            update_tmp_dic['raw_exp_id'] = raw_exp_id
        main_col.update_one({'_id': main_id}, {'$set': update_tmp_dic})

    def interaction_ori_cv(self,table_path, raw_exp_id):
        table_list = table_path.split(',')
        tar_dir = os.path.abspath(os.path.join(table_list[0],'..'))
        if os.path.exists(tar_dir+'/raw_statistic.xls'):
            self._add_exp_stat(tar_dir+'/raw_statistic.xls', raw_exp_id)

        if os.path.exists(tar_dir + '/raw_cv_percent.xls'):
            self._add_exp_cv(tar_dir+'/raw_cv_percent.xls',  raw_exp_id)


    def _add_metab_table_detail(self, main_id, table_path, table_type):
        if table_type in ['pos', 'neg', 'mix']:
            # collection = self.db['metab_table_' + table_type]
            collection = self.db['exp_' + table_type]
        else:
            self.bind_object.set_error("不存在的代谢表类型：%s", variables=(table_type), code="54702206")
        data_list = list()  # 存入表格中的信息，然后用insert_many批量导入
        abund_path = os.path.join(table_path, "metab_abund.txt")
        desc_path = os.path.join(table_path, "metab_desc.txt")
        rsd_path = os.path.join(table_path, "metab_abund.txt_rsd.txt")
        if os.path.exists(rsd_path):
            rsd = pd.read_table(rsd_path,sep='\t',header=None)
            rsd_map = dict(rsd.values.tolist())
        else:
            rsd_map = dict()


        with open(abund_path, 'rb') as f1 : #, open(desc_path, 'rb')) as (f1, f2):
            abund_lines = f1.readlines()
            self.sample_list = abund_lines[0].strip().split('\t')[1:]
            self.qc_list = []
            for s in self.sample_list:
                if re.match('^QC_*\d+$',s):
                    self.qc_list.append(s)
            desc_table = pd.read_table(desc_path,sep='\t',header=0)
            hmdb_head = 'HMDB_ID'
            cas_head = 'CAS number'
            if hmdb_head not in desc_table.columns:
                hmdb_head = 'Library ID'
            if cas_head not in desc_table.columns:
                cas_head = 'CAS ID'


            var_heads = {'ID':'o_id','m/z':'mz','Retention time':'rt','RT (min)':'rt','Mode':'mode','Adducts':'adducts','Formula':'formula',
                         'Fragmentation Score':'frag','Theoretical Fragmentation Score':'theo_frag',
                         'Mass Error (ppm)':'mass_error','HMDB_ID':'hmdb_id','Library ID':'hmdb_id',
                         'CAS number':'cas_id','CAS ID':'cas_id', "Score": "score"}


            current_var_heads = []
            for h in var_heads.keys():
                if h in desc_table.columns:
                    current_var_heads.append(h)

            for i in range(len(desc_table)) :
                abund_line = abund_lines[i+1]
                abund_line = abund_line.strip().split('\t')
                data = [
                    ('exp_id', main_id),
                    ('metab_id', desc_table['metab_id'][i]),
                    ('metab', check_metab(desc_table['Metabolite'][i])),
                    #('mode', desc_table['Mode'][i]),
                    # ('formula', check_info(desc_table['Formula'][i])),
                    # ('element', self.element(desc_table['Formula'][i])),
                    # ('mz', desc_line[4]),
                    # ('mz', ("%.4f" % float(desc_table['m/z'][i]))),
                    # ('rt', desc_line[5]),
                    # ('rt', ("%.4f" % float(desc_table['Retention time'][i]))),
                    ('kegg_id', check_info(desc_table['KEGG Compound ID'][i])),
                    # ('hmdb_id', check_info(desc_table[hmdb_head][i])),
                    # ('cas_id', check_info(desc_table[cas_head][i]))
                ]
                if not re.match('^metab_\d*$',desc_table['Metabolite'][i]):  #用于前端展示有代谢物名称的数据
                    data.append(('has_name',1))

                for c_h in current_var_heads:
                    if c_h == 'ID':
                        data.append((var_heads[c_h],desc_table['ID'][i]))
                    elif c_h == 'm/z':
                        try:
                            mz = "%.4f" % float(desc_table['m/z'][i])
                        except:
                            mz = desc_table['m/z'][i]
                        data.append((var_heads[c_h],mz) )
                    elif c_h =='Retention time' or c_h=='RT (min)':
                        data.append((var_heads[c_h], ("%.4f" % float(desc_table[c_h][i]))))
                    elif c_h == 'Adducts':
                        data.append((var_heads[c_h],desc_table['Adducts'][i]))
                    elif c_h == 'Formula':
                        data.append((var_heads[c_h], check_info(desc_table['Formula'][i])))
                        data.append(('element', self.element(desc_table['Formula'][i])))
                    elif c_h =='Fragmentation Score':
                        try:
                            score = float(desc_table['Fragmentation Score'][i])
                        except:
                            score = desc_table['Fragmentation Score'][i]
                        data.append((var_heads[c_h], score))
                    elif c_h == 'Theoretical Fragmentation Score':
                        try:
                            score = float(desc_table['Theoretical Fragmentation Score'][i])
                        except:
                            score = desc_table['Theoretical Fragmentation Score'][i]
                        data.append((var_heads[c_h], score))
                    elif c_h == 'Mass Error (ppm)':
                        data.append((var_heads[c_h], desc_table['Mass Error (ppm)'][i]))
                    # elif c_h == 'RSD':
                    #     try:
                    #         rsd = float(desc_table['RSD'][i])
                    #     except:
                    #         rsd = desc_table['RSD'][i]
                    #     data.append((var_heads[c_h], rsd))
                    elif c_h == cas_head:
                        data.append((var_heads[c_h], check_info(desc_table[cas_head][i])))
                    elif c_h == hmdb_head:
                        data.append((var_heads[c_h], check_info(desc_table[hmdb_head][i])))
                    elif c_h == 'Mode':
                        data.append((var_heads[c_h], desc_table['Mode'][i]))
                    elif c_h == 'Score':
                        data.append((var_heads[c_h], desc_table['Score'][i]))

                #samples_abund = []
                for index, one_sample in enumerate(self.sample_list):
                    data.append((one_sample, ("%.4f" % float(abund_line[index + 1]))))
                #     if one_sample in self.qc_list:
                #         samples_abund.append(float(abund_line[index + 1]))
                # if len(samples_abund) > 1:
                #     rsd = round(np.std(samples_abund)/np.mean(samples_abund),4)  # 增加rsd。样本数不同计算结果不一样
                #     data.append(('rsd',rsd))
                if desc_table['metab_id'][i] in rsd_map:
                    data.append(('rsd', rsd_map[desc_table['metab_id'][i]]))
                data = SON(data)
                data_list.append(data)
        try:
            collection.insert_many(data_list)
            c_var_head_mongo = []
            for c_h in current_var_heads:
                c_var_head_mongo.append(var_heads[c_h])
            # if 'formula' in c_var_head_mongo:
            #     c_var_head_mongo.append('element')  #v2.0 去掉
            c_var_head_mongo_str = ','.join(c_var_head_mongo)
            if len(self.qc_list) >1:
                c_var_head_mongo_str += ',rsd'
            collection_main = self.db['exp']
            collection_main.update({"_id":main_id},{"$set":{"var_head":c_var_head_mongo_str}})

        except Exception as e:
            self.bind_object.set_error("导入%s信息出错:%s", variables=(table_path, e), code="54702207")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % table_path)

    def element(self, formula):
        formula = check_info(formula)
        if formula == "-":
            return formula
        element_list = []
        for i in formula:
            if i in "+.[]-(); ":
                forward_letter = i
                continue
            if i == 'n':
                if forward_letter in [']', ')']:
                    continue
                elif forward_letter.isupper():
                    element_list[-1] += i
                else:
                    self.bind_object.logger.info("化合物中含有未处理的情况，需注意：%s" % formula)
            try:
                int(i)  # 去掉数字
                forward_letter = i
            except:
                forward_letter = i  # 用于判断n是在]后还是在元素后
                if i.isupper():
                    element_list.append(i)
                elif i.islower():
                    element_list[-1] += i
                else:
                    self.bind_object.set_error('分子式%s中发现特殊字符%s', variables=(formula, i), code="54702208")
        element_set = set(element_list)
        return_str = ':'.join(list(element_set))
        return_str = ':' + return_str + ':'
        return return_str

    def _add_exp_stat(self,table,exp_id):
        self.bind_object.logger.info('exp_stat 开始导入')
        data = pd.read_table(table,sep='\t')
        insert_data = []
        for index in data.index:
            cur = data.loc[index]
            tmp = {
                "ion_type" : cur['type'],
                "all_num" : cur['all'],
                "identified_num" : cur['has_name'],
                "hmdb_num": cur['hmdb'] ,
                "kegg_num" : cur['kegg'],
                "exp_id": exp_id
            }
            # if is_raw:
            #     tmp['type'] = 'raw'
            # else:
            #     tmp['type'] = 'deal'
            insert_data.append(tmp)
        self.db['exp_stat'].insert_many(insert_data)
        self.bind_object.logger.info('exp_stat 导入完成')

    def _add_exp_cv(self, table, exp_id):
        self.bind_object.logger.info("exp_cv 开始导入")
        data = pd.read_table(table,sep='\t',index_col='percent')
        insert_data = []
        for type in data.columns:
            tmp = dict([('per%s'%i, e) for i,e in enumerate(data[type],1)])
            tmp['exp_id'] = exp_id
            tmp['model'] = type
            tmp['group'] = 'QC'
            tmp['per0'] = 0
            # if is_raw:
            #     tmp['type'] = 'raw'
            # else:
            #     tmp['type'] = 'deal'
            insert_data.append(tmp)

        self.db['exp_cv'].insert_many(insert_data)
        self.bind_object.logger.info("exp_cv 完成导表")






    '''
    def check_info(self, info):
        if info in ['_', '-', ' ', '', '\r']:
            return '-'
        else:
            return info

    def check_metab(self, info):
        info = self.check_info(info)
        if info == '-':
            return info
        info = info.replace("α","|alpha|").replace("β", "|beta|").replace("ω","|omega|").replace("->","|right|")
        return info
    '''

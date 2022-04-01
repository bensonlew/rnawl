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


class ExpDiff(Base):
    def __init__(self, bind_object):
        super(ExpDiff, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_exp_diff(self, metab_table_id, table_type, name=None, main_id=None, params =None, diff_dir=None):
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
                #'name': name if name else 'ExpDiff_Origin_' + table_type,
                'name': name if name else 'ExpDiff_Origin',
                'params': params if params else '',
                'status': 'end',
                'main_id': '',
                'diff_dir': diff_dir if diff_dir else '',
                "metab_table_id": metab_table_id,
                "metab_table_main_id": metab_table_id,
                #"table_type": table_type,   # v2 remove
                "version" : "3.0"
            }

            try:
                collection = self.db['exp_diff']
                diff_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", diff_id, diff_id)
            except Exception, e:
                self.bind_object.set_error('导入diff主表异常:%s', variables=(e), code="54702001")
        else:
            self.update_table("main_id", main_id, main_id)
            if not diff_dir:
                self.bind_object.set_error('必须输入diff_dir', code="54702002")
            self.update_table("diff_dir", diff_dir, main_id)
            diff_id = main_id
        return diff_id

    @report_check
    def add_exp_diff_detail(self, diff_id, diff_file_dir, table_type):
        diff_id = self.check_id(diff_id, "diff_id")
        if not os.path.exists(diff_file_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(diff_file_dir), code="54702003")
        result = self.db['exp_diff'].find_one({'_id': diff_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(diff_id), code="54702004")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")
            metab_table_id = eval(result['params'])["metab_table"]
        all_files = os.listdir(diff_file_dir)
        if not all_files:
            self.bind_object.set_error('%s为空目录', variables=(diff_file_dir), code="54702005")
        data_list = []
        mydict = self.get_some_info(metab_table_id, table_type)
        #self.bind_object.logger.info(len(mydict))
        for each in all_files:
            if '.diff.exp' not in each:
                self.bind_object.logger.info('导表，skip: '+ each)
                continue
            group = each.split(".diff.exp")[0]
            each_path = os.path.join(diff_file_dir, each)
            self.bind_object.logger.info('开始导表：' + each_path)
            with open(each_path, 'rb') as f:
                head = f.next()
                for line in f:
                    line = line.strip().split('\t')
                    metab_id = line[0]
                    p_value = float(line[3])
                    fdr = float(line[4])
                    fc = float(line[5])
                    pls_vip = float(line[1])
                    opls_vip = float(line[2])
                    if fc != 0:
                        d_fc = 1/fc
                    else:
                        d_fc = 99999999
                    if mydict.has_key(metab_id):
                        metab = mydict[metab_id]['metab']
                        #metab_id = mydict[metab_id][0]
                        #formula = mydict[metab_id][1]
                        #element = mydict[metab_id][2]
                        insert_data = {
                            'diff_id': diff_id,
                            'metab': metab,
                            'metab_id': metab_id,
                            ##'formula': formula,  #v2.0
                            #'element': element,  #v2.0
                            'pls_vip': pls_vip,
                            'opls_vip': opls_vip,
                            'fc': fc,
                            'd_fc':d_fc,
                            'p_value': p_value,
                            'diff_group': group,
                            'fdr' :fdr,
                            'table_type' : table_type  #v2 add 201907
                        }
                        if 'formula' in mydict[metab_id].keys():
                            insert_data['formula'] = mydict[metab_id]['formula']
                        if 'element' in mydict[metab_id].keys():
                            insert_data['element'] = mydict[metab_id]['element']


                        data_list.append(insert_data)
        try:
            collection = self.db['exp_diff_detail']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(each, e), code="54702006")
        else:
            self.bind_object.logger.info("导入表格%s信息成功!" % each)




    def add_exp_diff_detail_plot(self, diff_id, diff_file_dir, table_type):
        diff_id = self.check_id(diff_id, "diff_id")
        if not os.path.exists(diff_file_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(diff_file_dir), code="54702003")
        result = self.db['exp_diff'].find_one({'_id': diff_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(diff_id), code="54702004")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")
            metab_table_id = eval(result['params'])["metab_table"]
        all_files = os.listdir(diff_file_dir)
        if not all_files:
            self.bind_object.set_error('%s为空目录', variables=(diff_file_dir), code="54702005")
        data_list = []
        for each in all_files:
            if "_plot.xls" not in each:
                continue
            group = each.split("_plot.xls")[0]
            each_path = os.path.join(diff_file_dir, each)
            self.bind_object.logger.info(each_path)
            data = pd.read_table(each_path,sep='\t',header=0)
            data.fillna('',inplace=True)
            g1 = re.sub('-mean$','',data.columns[1])
            g2 = re.sub('-mean$','',data.columns[3])
            head1 = [i+'(%s)'%g1 for i in ['min','Q1','Median','Q3','max']]
            head2 = [i+'(%s)'%g2 for i in ['min','Q1','Median','Q3','max']]
            for i in range(len(data)):
                box1 = ','.join(map(str,data.loc[i,head1].tolist()))
                box2 = ','.join(map(str,data.loc[i,head2].tolist()))
                insert_data = {
                    'diff_id': diff_id,
                    'metab_id': data['Metab'][i],
                    'diff_group': group,
                    'table_type' : table_type,  #v2 add 201907
                    #'type' : 'plot',
                    g1+'_std' : data[g1+'-sd'][i],
                    g2+'_std' : data[g2+'-sd'][i],
                    g1+'_mean' : data[g1+'-mean'][i],
                    g2+'_mean' : data[g2+'-mean'][i],
                    g1+'_box' : box1,
                    g2+'_box' : box2
                }
                data_list.append(insert_data)
        try:
            collection = self.db['exp_diff_plot']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(each, e), code="54702006")
        else:
            self.bind_object.logger.info("导入表格%s信息成功!" % each)




    @report_check
    def add_exp_diff_model(self, diff_id, pls_dir, type='None'):
        diff_id = self.check_id(diff_id, "diff_id")
        if not os.path.exists(pls_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pls_dir), code="54702007")
        data_list = []
        result = self.db['exp_diff'].find_one({'_id': diff_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(diff_id), code="54702008")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")
        groups = os.listdir(pls_dir)
        for eachgroup in groups:
            pls_dict = {}
            intercept_files = glob.glob(pls_dir + "/" + eachgroup + '/*intercept.xls')
            for eachfile in intercept_files:
                model = self.get_model(eachfile)
                with open(eachfile, "r") as f:
                    head = f.next()
                    r2 = f.next().strip().split("\t")[1]
                    q2 = f.next().strip().split("\t")[1]
                    pls_dict[model + "_r2"] = float(r2)
                    pls_dict[model + "_q2"] = float(q2)
            model_files = glob.glob(pls_dir + "/" + eachgroup + '/*model.xls')
            self.bind_object.logger.info(model_files)
            for eachfile in model_files:
                model = self.get_model(eachfile)
                if model == "pca":
                    self.bind_object.logger.info(eachfile)
                    with open(eachfile, 'rb') as f:
                        head = f.next()
                        for line in f:
                            line = line.strip().split('\t')
                            a = line[0]
                            r2x = float(line[1])
                            r2x_cum = float(line[2])
                            insert_data = {
                                'diff_id': diff_id,
                                'diff_group': eachgroup,
                                'model': model,
                                'a': a,
                                'r2x': r2x,
                                'r2x_cum': r2x_cum
                            }
                            data_list.append(insert_data)
                else:
                    data_list = self.diff_model(diff_id, data_list, eachfile, model, eachgroup, pls_dict)
        try:
            collection = self.db['exp_diff_model']
            if type != 'None':
                new_data_list = []
                for i in data_list:
                    i['table_type'] = type
                    new_data_list.append(i)
                collection.insert_many(new_data_list)
            else:
                collection.insert_many(data_list)

        except Exception as e:
            self.bind_object.set_error("导入表格exp_diff_model信息出错:%s" , variables=(e), code="54702009")
        else:
            self.bind_object.logger.info("导入表格exp_diff_model信息成功!")

    @report_check
    def diff_model(self, diff_id, data_list, eachfile, model, diff_group, pls_dict):
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
                    'diff_id': diff_id,
                    'diff_group': diff_group,
                    'model': model,
                    'a': a,
                    'r2x': r2x,
                    'r2x_cum': r2x_cum,
                    'q2': q2,
                    'q2_cum': q2_cum,
                    'r2y_cum': r2y_cum,
                    'r2y': r2y,
                    'r2': pls_dict[model + "_r2"],
                    'pq2': pls_dict[model + "_q2"],
                }
                data_list.append(insert_data)
        return data_list

    @report_check
    def add_exp_diff_comp(self, diff_id, pls_dir, group_file,type='None'):
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
        result = self.db['exp_diff'].find_one({'_id': diff_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(diff_id), code="54702012")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")
        groups = os.listdir(pls_dir)
        for eachgroup in groups:
            comp_files = glob.glob(pls_dir + "/" + eachgroup + '/*sites.xls')
            for eachfile in comp_files:
                model = self.get_model(eachfile)
                with open(eachfile, 'rb') as f:
                    head = f.next()
                    for line in f:
                        line = line.strip().split('\t')
                        name = line[0]
                        group = group_dict[name]
                        pc1 = float(line[1])
                        pc2 = float(line[2])
                        insert_data = {
                            'diff_id': diff_id,
                            'diff_group': eachgroup,
                            'model': model,
                            'group': group,
                            'name': name,
                            'pc1': pc1,
                            'pc2': pc2,
                            'type': "specimen"
                        }
                        data_list.append(insert_data)
            ellipse_files = glob.glob(pls_dir + "/" + eachgroup + '/*ellipse.xls')
            for eachfile in ellipse_files:
                model = self.get_model(eachfile)
                with open(eachfile, 'rb') as f:
                    head = f.next()
                    for line in f:
                        line = line.strip().split('\t')
                        if len(line) < 7:
                            continue
                        group = line[0]
                        m1 = line[1]
                        m2 = line[2]
                        c11 = line[3]
                        c12 = line[4]
                        c21 = line[5]
                        c22 = line[6]
                        ellipse = [m1, m2, c11, c12, c21, c22]
                        insert_data = {
                            'diff_id': diff_id,
                            'diff_group': eachgroup,
                            'model': model,
                            'group': group,
                            'name': "p1p1",
                            'type': "circle",
                            'ellipse': ellipse
                        }
                        data_list.append(insert_data)
        try:
            collection = self.db['exp_diff_comp']
            if type != 'None':
                new_data_list = []
                for i in data_list:
                    i['table_type'] = type
                    new_data_list.append(i)
                if new_data_list:
                    collection.insert_many(new_data_list)
            else:
                if data_list:
                    collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格exp_diff_comp信息出错:%s" , variables=(e), code="54702013")
        else:
            self.bind_object.logger.info("导入表格exp_diff_comp信息成功!")

    @report_check
    def add_exp_diff_bar(self, diff_id, pls_dir, type='None'):
        diff_id = self.check_id(diff_id, "diff_id")
        if not os.path.exists(pls_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pls_dir), code="54702014")
        data_list = []
        result = self.db['exp_diff'].find_one({'_id': diff_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(diff_id), code="54702015")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")
        groups = os.listdir(pls_dir)
        for eachgroup in groups:
            comp_files = glob.glob(pls_dir + "/" + eachgroup + '/*PLS-DA.model.xls')
            for eachfile in comp_files:
                model = self.get_model(eachfile)
                number = 0
                if model == "plsda":
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
                                'diff_id': diff_id,
                                'diff_group': eachgroup,
                                'model': model,
                                'x': comp,
                                'q2': q2_cum,
                                'r2y': r2y_cum
                            }
                            data_list.append(insert_data)
                elif model == "oplsda":
                    with open(eachfile, 'rb') as f:
                        head = f.next()
                        for line in f:
                            line = line.strip().split('\t')
                            comp = "comp(1+" + str(number) + ")"
                            r2y_cum = float(line[4])
                            q2_cum = float(line[6])
                            if number == 0:
                                r2y_cum_p = r2y_cum
                                q2_cum_p = q2_cum
                                r2y_cum_pq = r2y_cum
                                q2_cum_pq = q2_cum
                            else:
                                r2y_cum_pq =  r2y_cum_p + r2y_cum
                                q2_cum_pq = q2_cum_p + q2_cum
                            number += 1
                            if float(line[4]) <0 or float(line[6]) <0:
                                self.bind_object.logger.info("贡献度出现负值，建议重新选择参数运行！")
                            if not line[0] == "sum":
                                insert_data = {
                                    'diff_id': diff_id,
                                    'diff_group': eachgroup,
                                    'model': model,
                                    'x': comp,
                                    'q2': q2_cum_pq,
                                    'r2y': r2y_cum_pq
                                }
                                data_list.append(insert_data)
        try:
            collection = self.db['exp_diff_bar']
            if type != 'None':
                new_data_list = []
                for i in data_list:
                    i['table_type'] = type
                    new_data_list.append(i)
                collection.insert_many(new_data_list)
            else:
                collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格exp_diff_bar信息出错:%s" , variables=(e), code="54702018")
        else:
            self.bind_object.logger.info("导入表格exp_diff_bar信息成功!")

    @report_check
    def add_exp_diff_scatter(self, diff_id, pls_dir,type='None'):
        diff_id = self.check_id(diff_id, "diff_id")
        if not os.path.exists(pls_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pls_dir), code="54702019")
        data_list = []
        result = self.db['exp_diff'].find_one({'_id': diff_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(diff_id), code="54702020")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")
        groups = os.listdir(pls_dir)
        mytest = []
        for eachgroup in groups:
            scatter_files = glob.glob(pls_dir + "/" + eachgroup + '/*permMN.xls')
            for eachfile in scatter_files:
                number_r = 0
                number_q = 0
                model = self.get_model(eachfile)
                with open(eachfile, 'rb') as f:
                    head = f.next()
                    for line in f:
                        line = line.strip().split('\t')
                        name = line[0]
                        r2y = float(line[0])
                        q2 = float(line[1])
                        sim = float(line[2])
                        insert_data = {
                            'diff_id': diff_id,
                            'diff_group': eachgroup,
                            'model': model,
                            'x':sim,
                            'y':r2y,
                            'geom':"r2"
                        }
                        data_list.append(insert_data)
                        if sim == 1 and number_r == 0:
                            number_r += 1
                            insert_data = {
                                'diff_id': diff_id,
                                'diff_group': eachgroup,
                                'model': model,
                                'x':sim,
                                'y':r2y,
                                'geom':"lr"
                            }
                            data_list.append(insert_data)
                        insert_data = {
                            'diff_id': diff_id,
                            'diff_group': eachgroup,
                            'model': model,
                            'x':sim,
                            'y':q2,
                            'geom':"q2"
                        }
                        data_list.append(insert_data)
                        if sim == 1 and number_q == 0:
                            number_q += 1
                            insert_data = {
                                'diff_id': diff_id,
                                'diff_group': eachgroup,
                                'model': model,
                                'x':sim,
                                'y':q2,
                                'geom':"lq"
                            }
                            data_list.append(insert_data)
            intercept_files = glob.glob(pls_dir + "/" + eachgroup + '/*intercept.xls')
            for eachfile in intercept_files:
                model = self.get_model(eachfile)
                with open(eachfile, "r") as f:
                    head = f.next()
                    r2 = f.next().strip().split("\t")[1]
                    q2 = f.next().strip().split("\t")[1]
                    insert_data = {
                        'diff_id': diff_id,
                        'diff_group': eachgroup,
                        'model': model,
                        'x': 0,
                        'y': float(r2),
                        'geom':"lr"
                    }
                    data_list.append(insert_data)
                    insert_data = {
                        'diff_id': diff_id,
                        'diff_group': eachgroup,
                        'model': model,
                        'x': 0,
                        'y': float(q2),
                        'geom':"lq"
                    }
                    data_list.append(insert_data)
        try:
            collection = self.db['exp_diff_scatter']
            if type != 'None':
                new_data_list = []
                for i in data_list:
                    i['table_type'] = type
                    new_data_list.append(i)
                collection.insert_many(new_data_list)
            else:
                collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格exp_diff_scatter信息出错:%s" , variables=(e), code="54702021")
        else:
            self.bind_object.logger.info("导入表格exp_diff_scatter信息成功!")

    ##20200316
    @report_check
    def add_exp_diff_load(self, diff_id, pls_dir, type='None'):
        diff_id = self.check_id(diff_id, "diff_id")
        data_list=[]
        groups = os.listdir(pls_dir)

        for eachgroup in groups:
            load_files = glob.glob(pls_dir + "/" + eachgroup + '/*loadings.xls')
            for eachfile in load_files:
                model = self.get_model(eachfile)
                with open(eachfile, 'rb') as f:
                    head = f.next()
                    for line in f:
                        line = line.strip().split('\t')
                        name = line[0]
                        insert_data = {
                            "metab" : name,
                            "diff_group" : eachgroup,
                            "model" :model,
                            "x" : float(line[1]),
                            "y" : float(line[2]),
                            "diff_id" :diff_id
                        }
                        data_list.append(insert_data)
        try:
            collection = self.db['exp_diff_loadplot']
            if type != 'None':
                for i in data_list:
                    i['table_type'] = type
                collection.insert_many(data_list)
            else:
                collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格exp_diff_loadplot信息出错:%s" , variables=(e))
        else:
            self.bind_object.logger.info("导入表格exp_diff_loadplot信息成功!")

    #20200316
    @report_check
    def add_exp_diff_splot(self, diff_id, pls_dir, type='None'):
        diff_id = self.check_id(diff_id, "diff_id")
        data_list=[]
        groups = os.listdir(pls_dir)

        for eachgroup in groups:
            files = glob.glob(pls_dir + "/" + eachgroup + '/*.splots.xls')
            for eachfile in files:
                model = self.get_model(eachfile)
                if model != 'oplsda':
                    continue
                vip_file = pls_dir + "/" + eachgroup + '/OPLS-DA.vips.xls'
                metab_name_map_vip = self.get_metab_vip(vip_file)
                with open(eachfile, 'rb') as f:
                    head = f.next()
                    for line in f:
                        line = line.strip().split('\t')
                        name = line[0]
                        insert_data = {
                            "metab" : name,
                            "diff_group" : eachgroup,
                            "model" :model,
                            "x" : float(line[1]),
                            "y" : float(line[2]),
                            "diff_id" :diff_id,
                            'vip' : float(metab_name_map_vip[name])
                        }
                        data_list.append(insert_data)
        try:
            collection = self.db['exp_diff_splot']
            if type != 'None':
                for i in data_list:
                    i['table_type'] = type
                collection.insert_many(data_list)
            else:
                collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入表格exp_diff_splot信息出错:%s" , variables=(e))
        else:
            self.bind_object.logger.info("导入表格exp_diff_splot信息成功!")

    @report_check
    def get_metab_vip(self,file):
        retd = {}
        with open(file) as f:
            for line in f:
                line = line.strip()
                spline = line.split('\t')
                retd[spline[0]] = spline[1]
        return retd


    @report_check
    def update_table(self, str, name, main_table_id):
        try:
            self.db['exp_diff'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {str: name}})
        except Exception as e:
            self.bind_object.set_error('更新exp_diff主表%s字段出错:%s', variables=(str, e), code="54702022")

    @report_check
    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type), code="54702023")
        return object_id

    @report_check
    def get_some_info(self, metab_table_id, table_type):
        self.mydict = {}
        coll = "exp_" + table_type
        metab_table_id =  self.check_id(metab_table_id, "metab_table_id")
        exp_find = self.db["exp"].find_one({"_id":metab_table_id})
        main_id = exp_find["main_id"]
        main_id = self.check_id(main_id, "main_id")
        #self.bind_object.logger.info("find metab,formula,element")
        self.bind_object.logger.info(main_id)
        if exp_find['is_raw'] == 1 and table_type=='mix':
            result1 = self.db['exp_pos'].find({'exp_id': main_id})
            result2 = self.db['exp_neg'].find({'exp_id': main_id})
            results = [result1, result2]
        else:
            results = [self.db[coll].find({'exp_id': main_id})]

        for result in results:
            for each in result:
                #each_list = []
                each_dic = {}
                metab_id = each['metab_id']
                #formula = each['formula']
                #element = each['element']
                #metab = each['metab']
                #each_list = [metab] # formula, element]
                each_dic['metab'] = each['metab']
                if 'formula' in each.keys():
                    each_dic['formula'] = each['formula']
                if 'element' in each.keys():
                    each_dic['element'] = each['element']

                #self.mydict[metab_id] = each_list
                self.mydict[metab_id] = each_dic
        return self.mydict

    @report_check
    def get_model(self, eachfile):
        if "PCA." in eachfile:
            model = "pca"
        elif "OPLS-DA." in eachfile:
            model = "oplsda"
        else:
            model = "plsda"
        return model



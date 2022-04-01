# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import os
import json
import datetime
import re
import copy
import numpy as np
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from mainapp.libs.param_pack import group_detail_sort


class BetaDiversityEnv(Base):
    """
    导入环境因子分析的MongoDB表  dbrda,rda_cca模块
    """
    def __init__(self, bind_object):
        super(BetaDiversityEnv, self).__init__(bind_object)
        self._project_type = 'metaasv'
        self._tables = []  # 记录存入了哪些表格

    @staticmethod
    def get_main_table_name(self, analysis_type):
        if analysis_type == 'dbrda':
            return 'db-RDA'
        elif analysis_type == 'rda_cca':
            return 'RDA/CCA'
        else:
            self.bind_object.set_error('错误的分析类型')

    @report_check
    def add_beta_multi_analysis_result(self, dir_path, analysis, main_id=None, main=False, group_id=None,
                                     otu_id=None, name=None, params=None, level=9, remove=None,
                                       spname_spid=None, group_file=None):
        self._tables = []  # 记录存入了哪些表格
        if level and level not in range(1, 10):
            self.bind_object.set_error("level水平错误")
        task_id = self.bind_object.sheet.id
        self.task_id = "_".join(task_id.split('_')[0:2])
        if not isinstance(group_id, ObjectId) and group_id is not None:
            group_id = ObjectId(group_id)
        else:
            if 'group_id' in self.bind_object.sheet.data['options']:
                group_id = self.bind_object.option('group_id')
                if group_id not in ['all', 'All', 'ALL']:
                    group_id = ObjectId(group_id)  # 仅仅即时计算直接绑定workflow对象
            else:
                group_id = 'all'
        if isinstance(otu_id, ObjectId):
            pass
        # elif otu_id is not None:
        #     otu_id = ObjectId(otu_id)
        # else:
        #     otu_id = ObjectId(self.bind_object.option('otu_id'))  # 仅仅即时计算直接绑定workflow对象
        _main_collection = self.db[analysis]
        if main:
            if not isinstance(params, dict):
                params_dict = json.loads(params)
            else:
                params_dict = params
            params_dict['asv_id'] = str(otu_id)  # otu_id在再metabase中不可用
            if spname_spid:
                params_dict['group_id'] = group_id
                group_detail = {'all': [str(i) for i in spname_spid.values()]}
                params_dict['group_detail'] = group_detail_sort(group_detail)
            insert_mongo_json = {
                'project_sn': self.bind_object.sheet.project_sn,
                'task_id': self.task_id,
                'asv_id': otu_id,
                'level_id': int(level),
                'name': BetaDiversityEnv.get_main_table_name(analysis) + '_Origin',
                'group_id': group_id,
                'params': (json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
                           if isinstance(params, dict) else params),
                'status': 'end',
                'desc': '',
                'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            }
            multi_analysis_id = _main_collection.insert_one(insert_mongo_json).inserted_id
            main_id = multi_analysis_id
        else:
            if not main_id:
                self.bind_object.set_error('不写入主表时，需要提供主表ID')
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        result = _main_collection.find_one({'_id': main_id})
        if result:
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
            else:
                main_id = main_id
            if analysis == 'rda_cca':
                if "cca_sites.xls" in os.listdir(dir_path.rstrip('/') + '/'):
                    rda_cca = "cca"
                else:
                    rda_cca = "rda"
                site_path = dir_path.rstrip('/') + '/' + rda_cca + '_sites.xls'
                importance_path = dir_path.rstrip('/') + '/' + rda_cca + '_importance.xls'
                dca_path = dir_path.rstrip('/') + '/' + 'dca.xls'
                plot_species_path = dir_path.rstrip('/') + '/' + rda_cca + '_plot_species_data.xls'
                envfit_path = dir_path.rstrip('/') + '/' + rda_cca + '_envfit.xls'
                if os.path.exists(envfit_path):
                    if  len(open(envfit_path).readlines()) < 2:
                        os.remove(envfit_path)
                    else:##correlation_key="pcoa_id",coll_name="pcoa_detail",main_coll="pcoa"
                        self.insert_table(envfit_path, 'envfit', update_id=main_id, correlation_key="rda_id",coll_name="rda_cca_table",main_coll="rda_cca")
                        self.insert_table_detail(envfit_path, 'envfit', update_id=main_id, correlation_key="rda_id",coll_name="rda_cca_detail",main_coll="rda_cca",type="arrow")
                self.insert_table(plot_species_path, 'plot_species', update_id=main_id, correlation_key="rda_id",coll_name="rda_cca_table",main_coll="rda_cca")
                self.insert_table_detail(site_path, 'specimen', update_id=main_id, correlation_key="rda_id",coll_name="rda_cca_detail",main_coll="rda_cca", type="scatter")
                self.insert_table(importance_path, 'importance', update_id=main_id, correlation_key="rda_id",coll_name="rda_cca_table",main_coll="rda_cca")
                self.insert_table(dca_path, 'dca', update_id=main_id, correlation_key="rda_id",coll_name="rda_cca_table",main_coll="rda_cca")
                filelist = os.listdir(dir_path.rstrip('/'))
                if (rda_cca + '_centroids.xls') in filelist:
                    env_fac_path = dir_path.rstrip('/') + '/' + rda_cca + '_centroids.xls'
                    self.insert_table_detail(env_fac_path, 'factor', update_id=main_id)
                if (rda_cca + '_biplot.xls') in filelist:
                    env_vec_path = dir_path.rstrip('/') + '/' + rda_cca + '_biplot.xls'
                    self.insert_table_detail(env_vec_path, 'amount', update_id=main_id, correlation_key="rda_id",coll_name="rda_cca_detail",main_coll="rda_cca", type="arrow")
                #分组椭圆数据
                if group_id not in ['all', 'All', 'ALL']:
                    circle_path = dir_path.rstrip('/') + '/ellipse.xls'
                    if os.path.exists(circle_path):
                        self.insert_scatter_plugin(circle_path, 'circle', update_id=ObjectId(main_id),correlation_key="rda_id",coll_name="rda_cca_plugin",main_coll="rda_cca")
                        os.remove(circle_path)
                self.update_main_collection(main_id, "rda_cca", "rda_type", rda_cca)
                self.bind_object.logger.info('beta_diversity:RDA/CCA分析结果导入数据库完成.')
            elif analysis == 'dbrda':
                site_path = dir_path.rstrip('/') + '/db_rda_sites.xls'
                self.insert_table_detail(site_path, 'specimen', update_id=main_id, correlation_key="dbrda_id",coll_name="dbrda_detail",main_coll="dbrda", type="scatter")
                filelist = os.listdir(dir_path.rstrip('/') + '/')
                # if 'db_rda_centroids.xls' in filelist:
                #     env_fac_path = dir_path.rstrip('/') + '/db_rda_centroids.xls'
                #     self.insert_table_detail(env_fac_path, 'factor', update_id=main_id)
                if 'db_rda_biplot.xls' in filelist:
                    env_vec_path = dir_path.rstrip('/') + '/db_rda_biplot.xls'
                    self.insert_table_detail(env_vec_path, 'amount', update_id=main_id, correlation_key="dbrda_id",coll_name="dbrda_detail",main_coll="dbrda", type="arrow")
                if 'db_rda_envfit.xls' in filelist:
                    envfit_path = dir_path.rstrip('/') + '/db_rda_envfit.xls'
                    if len(open(envfit_path).readlines()) < 2:
                        os.remove(envfit_path)
                    else:
                        self.insert_table(envfit_path, 'envfit', update_id=main_id, correlation_key="dbrda_id",coll_name="dbrda_table",main_coll="dbrda")
                        self.insert_table_detail(envfit_path, 'envfit', update_id=main_id, correlation_key="dbrda_id",coll_name="dbrda_detail",main_coll="dbrda",type="arrow")
                if 'db_rda_importance.xls' in filelist:          #guanqing.zou 20180417
                    importance_path = dir_path.rstrip('/') + '/db_rda_importance.xls'
                    self.insert_table(importance_path, 'importance', update_id=main_id, correlation_key="dbrda_id",coll_name="dbrda_table",main_coll="dbrda")
                plot_species_path = dir_path.rstrip('/') + '/db_rda_plot_species_data.xls'
                self.insert_table(plot_species_path, 'plot_species', update_id=main_id, correlation_key="dbrda_id",coll_name="dbrda_table",main_coll="dbrda")
                # species_path = dir_path.rstrip('/') + '/db_rda_species.xls'
                # self.insert_table_detail(species_path, 'species', update_id=main_id, split_fullname=True)
                if group_id not in ['all', 'All', 'ALL']:
                    circle_path = dir_path.rstrip('/') + '/ellipse.xls'
                    if os.path.exists(circle_path):
                        self.insert_scatter_plugin(circle_path, 'circle', update_id=ObjectId(main_id),correlation_key="dbrda_id",coll_name="dbrda_plugin",main_coll="dbrda")
                        os.remove(circle_path)
                self.bind_object.logger.info('beta_diversity:db_RDA分析结果导入数据库完成.')
        else:
            self.bind_object.logger.error('提供的_id：%s在分析记录中无法找到表, taskid: %s' % (str(main_id), self.task_id))
            self.bind_object.set_error("找不到表")
        return main_id

    def insert_table_detail(self, file_path, table_type, update_id,correlation_key=None,coll_name='pca_table',
                            main_coll='pca', type="scatter"):
        ##导入scatter表格
        self._tables.append(table_type)
        db = self.db
        collection = db[coll_name]
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            data_temp = []
            for line in all_lines[1:]:
                values = line.strip().split('\t')
                if main_coll in ["rda_cca"]:
                    columns = all_lines[0].strip().split('\t')[1:]
                    columns[-1] = "p_value"
                    columns[-2] = "r2"
                    insert_data = {
                        correlation_key: update_id,
                        "group": table_type,
                        "type":type,
                        "name": values[0].split(';')[-1].strip() ###样本名称
                        }
                    values_dict = dict(zip(columns, values[1:]))
                    data_temp.append(dict(insert_data, **values_dict))
                else:
                    columns = all_lines[0].strip().split('\t')
                    insert_data = {
                        correlation_key: update_id,
                        "group": table_type,
                        "type":type,
                        "name": values[0].split(';')[-1].strip() ###样本名称
                        }
                    values_dict = dict(zip(columns, values[1:]))
                    data_temp.append(dict(insert_data, **values_dict))
            try:
                collection.insert_many(data_temp)
                self.bind_object.logger.info('导入%s成功'%coll_name)
            except Exception, e:
                self.bind_object.logger.error('导入%s成功'%coll_name)
                self.bind_object.set_error("更新主表失败: %s" %e)
            try:
                main_collection = self.db[main_coll]
                scatter_data = {
                    "scatter_data": {"name": "name",
                    "data" : columns,
                    "category": ""}
                }
                scatter_data_json = json.dumps(scatter_data, sort_keys=True, separators=(',', ':'))
                arrow_data = {
                    "arrow_data": {"name": "name",
                    "data" : columns,
                    "category": ""}
                }
                arrow_data_json = json.dumps(arrow_data, sort_keys=True, separators=(',', ':'))
                main_collection.update_one({"_id": update_id}, {"$set": {"scatter_data": scatter_data_json,
                                                                         "arrow_data": arrow_data_json}})
                self.bind_object.logger.info('更新主表%s成功'%coll_name)
            except Exception, e:
                self.bind_object.logger.error('导入%s成功'%coll_name)
                self.bind_object.set_error("更新主表失败: %s" %e)

    def insert_table(self, file_path, table_type, update_id, correlation_key=None, coll_name='pca_table', main_coll='pca'):
        ##导入table的结果表
        self._tables.append(table_type)
        db = self.db
        collection = db[coll_name]
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            header = all_lines[0]
            data_temp = []
            for line in all_lines[1:]:
                if main_coll in ["rda_cca"] and table_type in ["envfit"]:
                    head = re.split("\s+", header.strip())
                    values = re.split("\s+", line.strip())
                    insert_data = {
                        correlation_key: update_id,
                        "type": table_type,
                        "p_value": values[-1],
                        "r2": values[-2],
                    }
                    try:
                        insert_data["name"] =  eval(values[0]) ###样本名称
                    except:
                        insert_data["name"] =  str(values[0]) ###样本名称
                    try:
                        column = [eval(x) for x in head[0:2]]
                    except:
                        column = [str(x) for x in head[0:2]]
                    value_dict = dict(zip(column,values[1:3]))
                    data_temp.append(dict(insert_data, **value_dict))
                elif table_type in ["importance"]:
                    values = line.strip().split('\t')
                    if main_coll in ["rda_cca"]:
                        insert_data = {
                            correlation_key: update_id,
                            "name": values[0], ###样本名称
                            "type": table_type,
                            "proportion_variance": values[1]
                            }
                        data_temp.append(insert_data)
                    else:
                        head = header.strip().split("\t")
                        insert_data = {
                            correlation_key: update_id,
                            "name": values[0], ###样本名称
                            "type": table_type,
                            }
                        value_dict = dict(zip(head,values[1:]))
                        data_temp.append(dict(insert_data, **value_dict))
                elif table_type in ["plot_species"]:
                    if main_coll in ["rda_cca"]:
                        head = header.strip().split()
                        column = head[1:]
                        values = line.strip().split('\t')
                        # self.bind_object.logger.info('head:%s'%column)
                        insert_data = {
                            correlation_key: update_id,
                            "type":table_type,
                            "name": values[0].split(';')[-1].strip() ###物种名称
                            }
                        values_dict = dict(zip(head[1:], values[1:]))
                        data_temp.append(dict(insert_data, **values_dict))
                    else:
                        head = header.strip().split()
                        column = head
                        values = line.strip().split('\t')
                        self.bind_object.logger.info('head:%s'%column)
                        insert_data = {
                            correlation_key: update_id,
                            "type":table_type,
                            "name": values[0].split(';')[-1].strip() ###物种名称
                            }
                        values_dict = dict(zip(head, values[1:]))
                        data_temp.append(dict(insert_data, **values_dict))
                else:
                    column = header.strip().split()
                    values = line.strip().split('\t')
                    insert_data = {
                        correlation_key: update_id,
                        "type":table_type,
                        "name": values[0] ###样本名称
                        }
                    values_dict = dict(zip(column, values[1:]))
                    data_temp.append(dict(insert_data, **values_dict))
            try:
                collection.insert_many(data_temp)
                self.bind_object.logger.info('导入%s成功'%coll_name)
            except Exception, e:
                self.bind_object.set_error("导入%s失败: %s" %(coll_name,e))
            try:
                main_collection = self.db[main_coll]
                settled_params = {"software" : "R-3.3.1 (vegan)"}
                settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
                importance_table_data = {
                        "table_data": ["name", "proportion_variance"],
                        "condition": {"type": "importance"}
                    }
                importance_table_data_json = json.dumps(importance_table_data, sort_keys=True, separators=(',', ':'))
                if table_type in ["dca"]:
                    dca_table_data = {
                            "table_data": ["name"] + column,
                            "condition": {"type": "dca"}
                        }
                    dca_table_data_json = json.dumps(dca_table_data, sort_keys=True, separators=(',', ':'))
                    main_collection.update_one({"_id": update_id}, {"$set": {"dca_table_data": dca_table_data_json}})
                if table_type in ["envfit"]:
                    envif_table_data = {
                            "table_data": ["name"] + column + ["r2", "p_value"],
                            "condition": {"type": "arrow"}
                        }
                    envif_table_data_json = json.dumps(envif_table_data, sort_keys=True, separators=(',', ':'))
                    main_collection.update_one({"_id": update_id}, {"$set": {"envif_table_data": envif_table_data_json}})
                if table_type in ["plot_species"]:
                    species_scatter_data = {
                    "scatter_data": {"name": "name",
                    "data" : column,
                    "category": ""}
                    }
                    species_scatter_data_json = json.dumps(species_scatter_data, sort_keys=True, separators=(',', ':'))
                    main_collection.update_one({"_id": update_id}, {"$set": {"species_scatter_data": species_scatter_data_json}})
                main_collection.update_one({"_id": update_id}, {"$set": {"importance_table_data": importance_table_data_json,
                                                                        "settled_params": settled_params_json,
                                                                         "main_id": update_id}})
                self.bind_object.logger.info('更新主表%s成功'%coll_name)
            except Exception, e:
                self.bind_object.set_error("更新主表失败: %s" %e)

    def insert_scatter_plugin(self, file_path, table_type, update_id,correlation_key=None, coll_name='pca_table', main_coll='pca'):
        self._tables.append(table_type)
        db = self.db
        data_temp = []
        collection = db[coll_name]
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            columns = all_lines[0].rstrip().split('\t')[1:]
            for line in all_lines[1:]:
                values = line.rstrip().split('\t')
                group_name_list =  re.split(r"\d", values[0])
                group_name_number_list = re.findall(r"\w[1-9]",values[0])
                first = group_name_list[0] + group_name_number_list[0][1]
                last = group_name_list[1] + group_name_number_list[1][1]
                insert_data = {
                    correlation_key: update_id,
                    'type': "ellipse",
                    "method": table_type,
                    "name": first + "_"+ last,
                    "x": first,
                    "y": last
                }
                values_dict = dict(zip(columns, values[1:]))
                data_temp.append(dict(insert_data, **values_dict))
        try:
            collection.insert_many(data_temp)
            self.bind_object.logger.info('导入%s成功'%coll_name)
        except Exception, e:
            self.bind_object.logger.error('导入%s成功'%coll_name)
            self.bind_object.set_error("更新主表失败: %s" %e)
        try:
            main_collection = self.db[main_coll]
            ellipse_data = {
                    "ellipse_data": {"data": columns,
                    "data_option" : columns,
                    "condition": {"type": ["circle", "ci_circle"]}}
                }
            ellipse_data_json = json.dumps(ellipse_data, sort_keys=True, separators=(',', ':'))
            main_collection.update_one({"_id": update_id}, {"$set": {"ellipse_data": ellipse_data_json}})
            self.bind_object.logger.info('更新主表%s成功'%main_coll)
        except Exception, e:
            self.bind_object.logger.error('更新主表%s失败'%main_coll)
            self.bind_object.set_error("更新主表失败: %s" %e)

    def insert_text_detail(self, file_path, data_type, main_id,
                           coll_name='sg_beta_multi_analysis_json_detail', db=None):
        if not db:
            db = self.db
        collection = db[coll_name]
        with open(file_path, 'rb') as f:
            data = f.read()
            insert_data = {
                'multi_analysis_id': main_id,
                'type': data_type,
                'json_value': data
            }
            collection.insert_one(insert_data)

    def insert_tree_table(self, file_path, table_type, update_id,correlation_key=None,coll_name='pca_table',
                            main_coll='pca'):
        if update_id is None:
            self.bind_object.set_error("需提供dist_id!")
        else:
            if not isinstance(update_id, ObjectId):
                dist_id = ObjectId(update_id)
        data_list = []
        with open(file_path, 'r') as f:
            if table_type in ["tree"]:
                line = f.readline().strip()
                insert_data = {
                    "name": "",
                    "data": line,
                    "type": table_type,
                    correlation_key:dist_id
                }
                data_list.append(insert_data)
            elif table_type in ["column"]:
                lines = f.readlines()
                columns = lines[0].rstrip().split('\t')[1:]
                for line in lines[1:]:
                    species_name = line[0]
                    insert_data = {
                        "species_name": species_name,
                        "type": table_type,
                        correlation_key:dist_id
                    }
                    values_dict = dict(zip(columns, line[1:]))
                    data_list.append(dict(insert_data, **values_dict))
            try:
                collection = self.db[coll_name]
                collection.insert_many(data_list)
                if type in ["column"]:
                    column_data = {
                        "column_data": {"name": "species_name",
                        "data": columns,
                        "condition": {"type": "column"}}
                    }
                    column_data_json = json.dumps(column_data, sort_keys=True, separators=(',', ':'))
                    main_collection = self.db[main_coll]
                    main_collection.update_one({"_id": dist_id}, {"$set": {"specimen_list": columns,
                                                                           "column_data": column_data_json}})
                elif type in ["tree"]:
                    main_collection = self.db[main_coll]
                    tree_data = {
                        "tree_data": {"name": "name",
                        "condition": {"type": "tree"}}
                    }
                    tree_data_json = json.dumps(tree_data, sort_keys=True, separators=(',', ':'))
                    main_collection.update_one({"_id": dist_id}, {"$set": {"tree_data": tree_data_json}})
            except Exception, e:
                self.bind_object.logger.error("更新%s信息出错:%s" % (file_path, e))
            else:
                self.bind_object.logger.info("更新%s信息成功!" % file_path)

    def insert_box_detail(self, file_path, table_type, update_id,correlation_key=None,coll_name='pca_table',
                            main_coll='pca', group_detail=None):
        ##计算box表，因为不需要呈现结果，只是用于画box图，需要根据heatmap的数据进行重新计算
        self._tables.append(table_type)
        db = self.db
        collection = db[coll_name]
        if group_detail: ## 读取group和样本的对应关系
            with open(self.option('group_file').prop['path'], 'r') as f:
                lines = f.readlines()
                sample_dict = {}
                group_list = []
                for line in lines[1:]:## 获取分组列表和样本与分组对应关系列表
                    line = line.strip().split("\t")
                    if line[1] in group_list:
                        group_list.append(line[1])
                    sample_dict[line[0]] = line[1]
                group_dict = {}
                for group in group_list:## 转为分组与样本的对应关系列表
                    sample_list = []
                    for sample in sample_dict.keys():
                        if sample_dict[sample] in [group]:
                            if sample not in sample_list:
                                sample_list.append(sample)
                    group_dict[group] = sample_list

        with open(file_path, 'rb') as f:##读取样本和样本的对应关系的value
            all_lines = f.readlines()
            columns = all_lines[0].strip().split('\t')[1:]
            all_specimen_dict = {}
            for line in all_lines[1:]:
                values = line.strip().split('\t')
                values_dict = dict(zip(columns, values[1:]))
                all_specimen_dict[line] = values_dict
        data_list = []
        for group in group_list:##获取组内数据
            value_list = []
            group_samples_list = group_dict[group]
            insert_data = {
                correlation_key:update_id,
                "specimen": group
            }
            for sample in group_samples_list:
                for sample_r in group_samples_list:
                    value = all_specimen_dict[sample][sample_r]
                    value_list.append(value)
            unqiue_list = list(set(value_list))
            unqiue_list.pop(0.0)
            box_dict = self.get_box(unqiue_list)
            data_list.append(dict(insert_data, **box_dict))
        if len(group_list) == 1:##获取组间数据
            for group in group_list:
                value_list = []
                group_samples_list = group_dict[group]
                insert_data = {correlation_key:update_id,
                                "specimen": "Between"}
                for sample in group_samples_list:
                    for sample_r in group_samples_list:
                        value = all_specimen_dict[sample][sample_r]
                        value_list.append(value)
                unqiue_list = list(set(value_list))
                unqiue_list.pop(0.0)
                box_dict = self.get_box(unqiue_list)
                data_list.append(dict(insert_data, **box_dict))
        else:##获取组间的数据
            insert_data = {correlation_key:update_id,
                                "specimen": "Between"}
            value_list = []
            for group in group_list:
                new_group_list = copy.deepcopy(group_list)
                new_group_list.remove(group)
                group_samples_list = group_dict[group]
                left_group_sample_list = []
                for new_g in new_group_list:
                    left_group_sample_list += group_dict[new_g]
                for sp in group_samples_list:
                    for sp_name in left_group_sample_list:
                        value = all_specimen_dict[sp][sp_name]
                        value_list.append(value)
            unqiue_list = list(set(value_list))
            unqiue_list.pop(0.0)
            box_dict = self.get_box(unqiue_list)
            data_list.append(dict(insert_data, **box_dict))
        try:
            collection.insert_many(data_list)
            self.bind_object.logger.info('导入%s成功'%coll_name)
        except Exception, e:
            self.bind_object.logger.error('导入%s成功'%coll_name)
            self.bind_object.set_error("更新主表失败: %s" %e)
        try:
            main_collection = self.db[main_coll]
            box_data = {
                        "box_data": {"name": "specimen",
                        "data": ["min", "q1", "median", "q3", "max"]}
                    }
            box_data_json = json.dumps(box_data, sort_keys=True, separators=(',', ':'))
            main_collection.update_one({"_id": update_id}, {"$set": {"box_data": box_data_json}})
            self.bind_object.logger.info('更新%s成功'%coll_name)
        except Exception, e:
            self.bind_object.logger.error('导入%s成功'%coll_name)
            self.bind_object.set_error("更新主表失败: %s" %e)

    def get_box(self, list):
        """
        功能：计算list中的box数据
        :param list: 输入的list
        :return:
        """
        box={}
        filter_list = []
        new_list = [float(x) for x in list]
        sort_list = sorted(new_list, reverse=False)
        # number = len(sort_list)
        # q1 = sort_list[int((number + 1)/4)]
        # median = sort_list[int((number + 1)/2)]
        # q3 = sort_list[int(3*(number + 1)/4)]
        percentiles = np.array([25, 50, 75])
        percent_value = np.percentile(sort_list, percentiles)
        # q1 = sort_list[int((number + 1)/4)]
        # median = sort_list[int((number + 1)/2)]
        # q3 = sort_list[int(3*(number + 1)/4)]
        q1 = percent_value[0]
        median = percent_value[1]
        q3 = percent_value[2]
        q_mean = float(q3 - median)
        max_value = q3 + 1.5*q_mean
        min_value = q1 - 1.5*q_mean
        for i in sort_list:
            if i < min_value and i > max_value:
                filter_list.append(i)
        box["min"] = min_value
        box["q1"] = q1
        box["median"] = median
        box["q3"] = q3
        box["max"] = max_value
        box["filter_data"] = filter_list
        return box

    def update_main_collection(self, update_id, coll_name, key_collection,value, db=None):
        if not db:
            db = self.db
        collection = db[coll_name]
        collection.update_one({"_id": ObjectId(update_id)}, {"$set": {key_collection: value}})
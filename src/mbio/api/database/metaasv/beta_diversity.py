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
from mainapp.libs.param_pack import group_detail_sort,param_pack
from mbio.packages.metaasv.common_function import find_group_name



class BetaDiversity(Base):
    """
    导入无环境因子分析的MongoDB表  pca模块
    """
    def __init__(self, bind_object):
        super(BetaDiversity, self).__init__(bind_object)
        self._project_type = 'metaasv'
        self._tables = []  # 记录存入了哪些表格
        self.main_name = {
            'pca' : 'PCA',
            'pcoa' : 'PCoA',
            'nmds' : 'NMDS',
            'hcluster' : 'Hcluster',
        }

    @report_check
    def add_beta_multi_analysis_result(self, dir_path, analysis, main_id=None, main=False, group_id=None,
                                     otu_id=None, name=None, params=None, level=9, remove=None,
                                       spname_spid=None, group_file=None):
        self._tables = []  # 记录存入了哪些表格
        if level and level not in range(1, 10):
            self.bind_object.set_error("level水平错误")
        task_id = self.bind_object.sheet.id
        self.task_id = "_".join(task_id.split('_')[0:2])
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
                if group_id not in [None, "All", "all", "ALL"]:
                ## 调用common模块，功能将导入的分组方案返回group_detail
                    group_detail = find_group_name(self.task_id)
                else:
                    group_detail = {'All': [str(i) for i in spname_spid.values()]}
                params_dict['group_detail'] = group_detail_sort(group_detail)
            params = param_pack(params_dict)
            insert_mongo_json = {
                'project_sn': self.bind_object.sheet.project_sn,
                'task_id': self.task_id,
                'asv_id': ObjectId(otu_id),
                'level_id': int(level),
                'name': self.main_name[analysis] + '_Origin',
                'group_id': group_id,
                'params': params,
                'status': 'end',
                'desc': 'Job has been finished',
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
            if analysis == 'pca':
                site_path = dir_path.rstrip('/') + '/pca_sites.xls'
                importance_path = dir_path.rstrip('/') + '/pca_importance.xls'
                self.insert_table_detail(site_path, 'specimen', update_id=ObjectId(main_id), correlation_key="pca_id",coll_name="pca_detail",main_coll="pca")
                self.insert_table(importance_path, 'importance', update_id=ObjectId(main_id), correlation_key="pca_id",coll_name="pca_table",main_coll="pca")
                if group_id not in ['all', 'All', 'ALL']:
                    circle_path = dir_path.rstrip('/') + '/ellipse.xls'
                    if os.path.exists(circle_path) :
                        self.insert_scatter_plugin(circle_path, 'circle', update_id=ObjectId(main_id), correlation_key="pca_id",coll_name="pca_scatter_plugin",main_coll="pca")
                        os.remove(circle_path)
                self.bind_object.logger.info('PCA分析结果导入数据库完成!')
            elif analysis == 'pcoa':
                site_path = dir_path.rstrip('/') + '/pcoa_sites.xls'
                self.insert_table_detail(site_path, 'specimen', update_id=ObjectId(main_id), correlation_key="pcoa_id",coll_name="pcoa_detail",main_coll="pcoa")
                importance_path = dir_path.rstrip('/') + '/pcoa_eigenvaluespre.xls'
                self.insert_table(importance_path, 'importance', update_id=ObjectId(main_id), correlation_key="pcoa_id",coll_name="pcoa_table",main_coll="pcoa")
                if group_id not in ['all', 'All', 'ALL']:
                    circle_path = dir_path.rstrip('/') + '/ellipse.xls'
                    if os.path.exists(circle_path) :
                        self.insert_scatter_plugin(circle_path, 'circle', update_id=ObjectId(main_id),correlation_key="pcoa_id",coll_name="pcoa_scatter_plugin",main_coll="pcoa")
                        os.remove(circle_path)
                self.bind_object.logger.info('PCoA分析结果导入数据库完成!')
            elif analysis == 'nmds':
                site_path = dir_path.rstrip('/') + '/nmds_sites.xls'
                self.insert_table_detail(site_path, 'specimen', update_id=ObjectId(main_id), correlation_key="nmds_id",coll_name="nmds_detail",main_coll="nmds")
                if group_id not in ['all', 'All', 'ALL']:
                    circle_path = dir_path.rstrip('/') + '/ellipse.xls'
                    if os.path.exists(circle_path) :
                        self.insert_scatter_plugin(circle_path, 'circle', update_id=ObjectId(main_id), correlation_key="nmds_id",coll_name="nmds_scatter_plugin",main_coll="nmds")
                        os.remove(circle_path)
                nmds_stress = float(open(dir_path.rstrip('/') + '/nmds_stress.xls').readlines()[1])
                _main_collection.update_one({'_id': main_id}, {'$set': {'nmds_stress': nmds_stress}})
                self.bind_object.logger.info('NMDS分析结果导入数据库完成!')
            elif analysis == 'hcluster':
                listdirs = os.listdir(dir_path)
                for file in listdirs:
                    file_path = os.path.join(dir_path, file)
                    if file in ["hcluster.tre"]:
                        self.insert_tree_table(file_path, 'tree', update_id=ObjectId(main_id), correlation_key="hcluster_id",coll_name="hcluster_tree_bar",main_coll="hcluster")
                    elif file in ["barplot_table.xls"]:
                        self.insert_tree_table(file_path, 'column', update_id=ObjectId(main_id), correlation_key="hcluster_id",coll_name="hcluster_tree_bar",main_coll="hcluster")
                    else:
                        self.insert_table_detail(file_path, 'specimen', update_id=ObjectId(main_id), correlation_key="hcluster_id",coll_name="hcluster_heatmap",main_coll="hcluster")
                        self.insert_box_detail(file_path, 'specimen', update_id=ObjectId(main_id), correlation_key="hcluster_id",coll_name="hcluster_box",main_coll="hcluster", group_detail=group_file)

            else:
                self.bind_object.set_error('提供的analysis：%s不存在')
            # self.insert_main_tables(self._tables, update_id=main_id)
        else:
            self.bind_object.logger.error('提供的_id：%s在分析记录中无法找到表, taskid: %s' % (str(main_id), self.task_id))
            self.bind_object.set_error("找不到表")
        return main_id

    def insert_table_detail(self, file_path, table_type, update_id,correlation_key=None,coll_name='pca_table',
                            main_coll='pca'):
        ##导入详情的结果表
        self._tables.append(table_type)
        db = self.db
        collection = db[coll_name]
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            if main_coll in ["nmds"]:
                columns = ["NMDS1", "NMDS2"]
            elif main_coll in ["hcluster", "pcoa"]:
                columns = all_lines[0].strip().split('\t')
            else:
                columns = all_lines[0].strip().split('\t')[1:]
            data_temp = []
            for line in all_lines[1:]:
                values = line.strip().split('\t')
                insert_data = {
                    correlation_key: update_id,
                    table_type: values[0].split(';')[-1].strip() ###样本名称
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
                settled_params = {"software" : "R-3.3.1 (vegan)"}
                settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
                if main_coll in ["hcluster"]:
                    if table_type in ["specimen"]:
                        heatmap_data = {
                            "heatmap_data": {"name": "specimen",
                            "data" : columns,
                            "category": ""}
                            }
                        heatmap_data_json = json.dumps(heatmap_data, sort_keys=True, separators=(',', ':'))
                        main_collection.update_one({"_id": update_id}, {"$set": {"main_id": update_id,
                                                                             "settled_params": settled_params_json,
                                                                             "heatmap_data": heatmap_data_json}})
                    elif table_type in ["column"]:
                        column_data = {
                            "column_data": {"name": "specimen",
                            "data" : columns,
                            "category": "",
                            "condition": {"type": "column"}}
                            }
                        specimen_sort = columns
                        column_data_json = json.dumps(column_data, sort_keys=True, separators=(',', ':'))
                        main_collection.update_one({"_id": update_id}, {"$set": {"main_id": update_id,
                                                                            "specimen_sort":specimen_sort,
                                                                             "settled_params": settled_params_json,
                                                                             "column_data": column_data_json}})
                else:
                    scatter_data = {
                        "scatter_data": {"name": "specimen",
                        "data" : columns,
                        "category": ""}
                    }
                    scatter_data_json = json.dumps(scatter_data, sort_keys=True, separators=(',', ':'))
                    main_collection.update_one({"_id": update_id}, {"$set": {"main_id": update_id,
                                                                             "settled_params": settled_params_json,
                                                                    "scatter_data": scatter_data_json}})
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
            data_temp = []
            for line in all_lines[1:]:
                values = line.strip().split('\t')
                insert_data = {
                    correlation_key: update_id,
                    "name": values[0].split(';')[-1].strip(), ###样本名称
                    "type": table_type,
                    "proportion_variance": float(values[1])
                    }
                data_temp.append(insert_data)
            try:
                collection.insert_many(data_temp)
                self.bind_object.logger.info('导入%s成功'%coll_name)
            except Exception, e:
                self.bind_object.logger.error('导入%s成功'%coll_name)
                self.bind_object.set_error("更新主表失败: %s" %e)
            try:
                main_collection = self.db[main_coll]
                importance_table_data = {
                        "table_data": ["name", "proportion_variance"],
                        "condition": {"type": "importance"}
                    }
                importance_table_data_json = json.dumps(importance_table_data, sort_keys=True, separators=(',', ':'))
                main_collection.update_one({"_id": update_id}, {"$set": {"importance_table_data": importance_table_data_json}})
                self.bind_object.logger.info('更新主表%s成功'%coll_name)
            except Exception, e:
                self.bind_object.logger.error('导入%s成功'%coll_name)
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
                if main_coll in ["pca", "pcoa"]:
                    name = str(values[0])
                    first =  "PC" + str(name[0])
                    last = "PC" + str(name[1])
                elif main_coll in ["nmds"]:
                    name = str(values[0])
                    first = "NMDS" + str(name[0])
                    last = "NMDS" + str(name[1])
                insert_data = {
                    correlation_key: update_id,
                    'type': "ellipse",
                    "method": table_type,
                    "name": first + "_"+last,
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
            else:
                dist_id = update_id
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
                    line = line.strip().split("\t")
                    species_name = line[0].strip().split("; ")[-1]
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
                settled_params = {"software" : "R-3.3.1"}
                settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
                if table_type in ["column"]:
                    column_data = {
                        "column_data": {"name": "species_name",
                        "data": columns,
                        "condition": {"type": "column"}}
                    }
                    column_data_json = json.dumps(column_data, sort_keys=True, separators=(',', ':'))
                    main_collection = self.db[main_coll]
                    main_collection.update_one({"_id": dist_id}, {"$set": {"specimen_list": columns,
                                                                           "main_id": dist_id,
                                                                           "column_data": column_data_json}})
                elif table_type in ["tree"]:
                    main_collection = self.db[main_coll]
                    tree_data = {
                        "tree_data": {"name": "name",
                        "condition": {"type": "tree"}}
                    }
                    tree_data_json = json.dumps(tree_data, sort_keys=True, separators=(',', ':'))
                    main_collection.update_one({"_id": dist_id}, {"$set": {"tree_data": tree_data_json,
                                                                           "main_id": dist_id,
                                                                           "settled_params": settled_params_json}})
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
        with open(group_detail, 'r') as f:## 读取group和样本的对应关系
            lines = f.readlines()
            sample_dict = {}
            group_list = []
            header = lines[0].strip().split("\t")
            group_name = header[1]
            if group_name in ["##empty_group##"]:
                for line in lines[1:]:## 获取分组列表和样本与分组对应关系列表
                    line = line.strip().split("\t")
                    if "all" not in group_list:
                        group_list.append("all")
                    sample_dict[line[0]] = "all"
            else:
                for line in lines[1:]:## 获取分组列表和样本与分组对应关系列表
                    line = line.strip().split("\t")
                    if line[1] not in group_list:
                        group_list.append(line[1])
                    sample_dict[line[0]] = line[1]
                self.bind_object.logger.info("group_list: {}".format(group_list))
            group_dict = {}
            for group in group_list:## 转为分组与样本的对应关系列表
                sample_list = []
                for sample in sample_dict.keys():
                    if sample_dict[sample] in [group]:
                        if sample not in sample_list:
                            sample_list.append(sample)
                group_dict[group] = sample_list
            self.bind_object.logger.info("group_dict: {}".format(group_dict))
        with open(file_path, 'rb') as f:##读取样本和样本的对应关系的value
            all_lines = f.readlines()
            columns = all_lines[0].strip().split('\t')
            all_specimen_dict = {}
            for line in all_lines[1:]:
                values = line.strip().split('\t')
                values_dict = dict(zip(columns,[float(x) for x in values[1:]]))
                all_specimen_dict[values[0]] = values_dict
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
            for value in unqiue_list:
                if value == 0.0:
                    unqiue_list.remove(value)
            if len(group_samples_list) == 1:
                unqiue_list = [1.0]
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
                for value1 in unqiue_list:
                    if value1 == 0.0:
                        unqiue_list.remove(value1)
                box_dict = self.get_box(unqiue_list)
                data_list.append(dict(insert_data, **box_dict))
        else:##获取组间的数据
            insert_data = {correlation_key:update_id,"specimen": "Between"}
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
            for value2 in unqiue_list:
                if value2 == 0.0:
                    unqiue_list.remove(value2)
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
        包含离异点的计算
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
                sort_list.remove(i)
                filter_list.append(i)
        ## 20200701 qinghcen.zhang fix this 原因是：原来将确定为异常值的界限作为最大值和最小值，这是不对的
        new_min_value = min(sort_list)
        new_max_value = max(sort_list)
        box["min"] = new_min_value
        box["q1"] = q1
        box["median"] = median
        box["q3"] = q3
        box["max"] = new_max_value
        box["filter_data"] = filter_list
        return box

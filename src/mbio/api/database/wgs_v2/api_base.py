# !/usr/bin/python
# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base, report_check
# from gevent import monkey; monkey.patch_socket()
from bson.objectid import ObjectId
from types import StringTypes
import datetime
import gevent
import json
import os
import time


class ApiBase(Base):
    def __init__(self, bind_object):
        """
        api基类，可以用于将列表导表入库，以及完成momgo表的查找，更新等
        __author__ = HONGDONG
        __last_modify__ = 20180206
        :param bind_object:
        """
        super(ApiBase, self).__init__(bind_object)
        self._project_type = "dna_wgs_v2"

    def col_find_one(self, collection, query_dic):
        """
        用于mongo查询一条记录，查询字段是一个字典,支持一个字段的查询，也支持多个字段的查询
        example: {"batch_id": ObjectID("5a782577a4e1af477e1ac081"), "type": "pt"}
        :param collection:
        :param query_dic:
        :return:
        """
        collection = self.db[collection]
        result = collection.find_one(query_dic)
        # if not result:
        #     print ("没有找到{}表中的{}结果，请检查".format(collection, query_dic.keys()))
        return result

    def col_find(self, collection, query_dic):
        """
        用于mongo查询多条记录，=查询字段是一个字典,支持一个字段的查询，也支持多个字段的查询
        example: {"batch_id": ObjectID("5a782577a4e1af477e1ac081"), "type": "pt"}
        :param collection:
        :param query_dic:
        :return:
        """
        collection = self.db[collection]
        return collection.find(query_dic)

    def col_insert_data(self, collection, data_list, is_show_log="true"):
        """
        插入多条记录，data_list样例格式[{"1":"2"},{"2":"3"},{"3": "4"}]
        这里只是抽象出，进行大规模的mongo插表,这里兼容了插入一条记录与多条记录
        :param collection:
        :param data_list:
        :param is_show_log: 用于设定是否要显示出插入成功的日志信息, 为false的时候就不显示成功的日志信息
        :return:
        """
        start = time.time()
        if not data_list:
            raise Exception("列表为空，不能进行后面的插表操作")
        table_id = None
        con = self.db[collection]
        record_num = len(data_list)
        try:
            if record_num > 5000000:
                for i in range(0, record_num, 4000000):
                    temp = data_list[i: i + 4000000]
                    con.insert_many(temp)
            else:
                if record_num >= 2:
                    con.insert_many(data_list)
                else:
                    table_id = con.insert_one(data_list[0]).inserted_id
        except Exception as e:
            if record_num >= 2:
                raise Exception("往{}插入多条记录失败{}".format(collection, e))
            else:
                raise Exception("往{}插入一条记录失败{}".format(collection, e))
        else:
            if is_show_log == 'true':
                if record_num >= 2:
                    print "往{}插入多条记录成功".format(collection)
                    # self.bind_object.logger.info("往{}插入多条记录成功".format(collection))
                else:
                    print "往{}插入一条记录成功".format(collection)
                    # self.bind_object.logger.info("往{}插入一条记录成功".format(collection))
                end = time.time()
                # self.bind_object.logger.info("文件导入mongo中花费的时间:{}".format(end - start))
                print ("文件导入mongo中花费的时间:{}".format(end - start))
            return table_id

    def gevet_insert_data(self, collection, data, step=100000):
        """
        实现并发导表，当一个文件比较大的时候，将文件存入data列表中，然后按照step步长进行分割导表
        :param collection:
        :param data:
        :param step:
        :return:
        """
        start = time.time()
        gevent_data = []
        data_length = len(data)
        for i in range(0, data_length, step):
            temp = data[i: i + step]
            gevent_data.append(gevent.spawn(self.col_insert_data, collection, temp, "true"))
        print "gevent data length: {}".format(len(gevent_data))
        gevent.joinall(gevent_data)
        end = time.time()
        print "gevent 导表花费时间：{}".format(end - start)

    def update_db_record(self, collection, query_dict, update_dict, is_show_log="true", upsert=True, multi=True):
        """
        用于更新表格的字段，暂时只是简单的更新一个表的对应字段，后面有需求再进行完善
        :param collection: 集合的名称
        :param query_dict: 查询的字段，可以有多个字段进行联合查询, 是一个字典
        :param update_dict: 表格中要更新进去的字段
        :param is_show_log: 用于设定是否要显示出更新成功的日志信息, 为false的时候就不显示成功的日志信息
        :param upsert: 表中没有查找到对应条件的记录，确认是否要更新 默认为True
        :param multi: 表中按照对应条件查找到多条记录，是否要全部更新 默认为True
        :return:
        """
        try:
            self.db[collection].update(query_dict, {'$set': update_dict}, upsert=upsert, multi=multi)
        except Exception as e:
            raise Exception("更新表格{}失败！{}".format(collection, e))
        else:
            if is_show_log == "true":
                print "更新{}到表格{}成功！".format(update_dict, collection)
                # self.bind_object.logger.info("更新表格{}成功！".format(collection))

    def remove_db_record(self, collection, query_dict):
        pass

    def add_main_table(self, collection, task_id, project_sn, params, name, desc, is_update="true"):
        """
        用于在基础workflow中，插入不同分析对应的主表，该模块格式固定, is_update 用于决定是否要将主表的状态更新为end，
        如果为false的时候就不更新为end
        :return:
        """
        data_list = []
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "start",
            "desc": desc,
            "task_id": task_id,
            "project_sn": project_sn,
            "name": name,
            "params": params
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data(collection, data_list)
        update_dict = {
            "main_id": main_id
        }
        if is_update == "true":
            update_dict.update({"status": "end"})
        self.update_db_record(collection, {"_id": main_id}, update_dict)
        return main_id

    def sg_curve(self, task_id, origin_id, name, categories, types, location=None, other_attr=None, ext_name=None):
        """
        用到导曲线图的数据
        :param task_id:
        :param origin_id: 上一级主表id，如sg_specimen_qc_id,sg_mapping_detail_id等
        :param name: 图的中文名称
        :param categories: x轴
        :param types: 类型，普通曲线
        :param location: 位置，对应的图英文名称，用于定位当前关联ID的哪个位置，通常用于解决一个结果表画相同的类型的图
        :param other_attr: 用于记录一些标记
        :param ext_name: 用于标题中显示名称的，当type为sample，改名时标题也需要改名：
        example：{"type":"sample","title":"样本原名"，"ext_title":"_1"}
        :return:
        """
        if not isinstance(origin_id, ObjectId):
            if isinstance(origin_id, StringTypes):
                origin_id = ObjectId(origin_id)
            else:
                raise Exception("origin_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        insert_data = {
            "task_id": task_id,
            "origin_id": origin_id,
            "name": name,
            "categories": categories,
            "type": types,
            "location": location if location else "",
            "other_attr": other_attr if other_attr else "",
            "ext_name": ext_name if ext_name else "",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_curve", data_list)
        return main_id

    def sg_curve_detail(self, curve_id, name, value, is_show_log="true"):
        """
        添加曲线图细节表
        :param curve_id: 曲线图的主表id
        :param name:  样本名称
        :param value:  值不能是字符串，值的顺序按照主表的categories里面对应的顺序
        :param is_show_log:  true  or false 用于决定是否显示插入成功的log日志，true为显示
        :return:
        """
        data_list = []
        insert_data = {
            "curve_id": curve_id,
            "name": name,
            "value": value
        }
        data_list.append(insert_data)
        self.col_insert_data("sg_curve_detail", data_list, is_show_log)

    def sg_bar(self, task_id, origin_id, name, categories, types, location=None, other_attr=None, ext_name=None):
        """
        用于导柱形图的数据
        :param task_id:
        :param origin_id: 上一级主表id，如sg_specimen_qc_id,sg_mapping_detail_id等
        :param name: 图的中文名称
        :param categories: x轴
        :param types: 类型，普通曲线
        :param location: 位置，对应的图英文名称，用于定位当前关联ID的哪个位置，通常用于解决一个结果表画相同的类型的图
        :param other_attr: 用于记录一些标记
        :param ext_name: 用于标题中显示名称的，当type为sample，改名时标题也需要改名：
        example：{"type":"sample","title":"样本原名"，"ext_title":"_1"}
        :return:
        """
        if not isinstance(origin_id, ObjectId):
            if isinstance(origin_id, StringTypes):
                origin_id = ObjectId(origin_id)
            else:
                raise Exception("origin_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        insert_data = {
            "task_id": task_id,
            "origin_id": origin_id,
            "name": name,
            "categories": categories,
            "type": types,
            "location": location if location else "",
            "other_attr": other_attr if other_attr else "",
            "ext_name": ext_name if ext_name else "",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_bar", data_list)
        return main_id

    def sg_bar_detail(self, bar_id, name, value, is_show_log="true", types='1', categories=None, tooltip=None):
        """
        添加柱形图细节表
        :param bar_id: 柱形图的主表id
        :param name:  图例
        :param value:  值不能是字符串，值的顺序按照主表的categories里面对应的顺序
        :param is_show_log:  true  or false 用于决定是否显示插入成功的log日志，true为显示
        :param types:  1, 2 ,3 对应三种柱状图，1，情况就是普通常见的柱状图，一个图例对应多个柱子，不截去前多少个记录，
        2对应的是有截图条件的，选择前20个记录，3对应的是一个图例对应一个柱子
        :param categories
        :param tooltip
        :return:
        """
        data_list = []
        insert_data = {
            "bar_id": bar_id,
            "name": name,
            "value": value
        }
        if types == "2":
            insert_data.update({"tooltip": tooltip})
        elif types == '3':
            insert_data.update({"categories": categories})
        data_list.append(insert_data)
        self.col_insert_data("sg_bar_detail", data_list, is_show_log)

    def sg_area(self, project_sn, task_id, origin_id, name):
        """
        添加基因组覆盖度分布主表
        :param origin_id: 主表sg_mapping表id
        :param name: 样本名称
        :param project_sn: 项目名
        :param task_id: 任务id
        """
        origin_id = self.check_objectid(origin_id)
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "origin_id": origin_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "name": name,
            "location": "mapping_multi_area"
        }
        main_id = self.db["sg_area"].insert_one(insert_data).inserted_id
        return main_id

    def sg_area_detail(self, area_id, legend, data):
        """
        添加基因组覆盖度分布细节表
        modified by hd 20180313 添加x轴统一化
        :param area_id: 主表sg_area表id
        :param data: 画图数据，dict,如{"chr1": [5, 4, 3], "chr2": [2, 3, 4]}
        :param legend:
        """
        area_id = self.check_objectid(area_id)
        data_list = []
        x_len = 0
        for keys in data.keys():   # 查找到x轴最长的值
            if len(data[keys]) > x_len:
                x_len = len(data[keys])
        # for k in data.keys():
        for k in legend:
            # if len(data[k]) < x_len:
            #     for m in range(len(data[k]) + 1, x_len + 1):
            #         data[k].append(None)
            insert_data = {
                "area_id": area_id,
                "name": k,
                "data": data[k]
            }
            data_list.append(insert_data)
        self.col_insert_data("sg_area_detail", data_list)
        self.db["sg_area"].update({"_id": area_id},
                                  {'$set': {"legend": legend, "x_max": x_len}}, upsert=True, multi=True)

    def sg_distribution(self, task_id, origin_id, name, types='1', start=None, end=None, location=None,
                        categories=None, data_type=None):
        """
        所有的染色体分布图
        :param task_id:
        :param origin_id:
        :param name:
        :param types:
        :param categories:
        :param location:
        :param start:
        :param end:
        :param data_type:
        :return:
        """
        origin_id = self.check_objectid(origin_id)
        insert_data = {
            "task_id": task_id,
            "origin_id": origin_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "name": name,
            "type": types,
            "location": location if location else "",
            "categories": categories if categories else [],
            "start": start,
            "end": end
        }
        if data_type:
            insert_data["data_type"] = data_type
        main_id = self.db["sg_distribution"].insert_one(insert_data).inserted_id
        return main_id

    def sg_distribution_detail(self, distribution_id, name, value, snp_value=None):
        """
        染色体分布图细节表
        :param distribution_id:
        :param name:
        :param value:
        :param snp_value:
        :return:
        """
        distribution_id = self.check_objectid(distribution_id)
        data_list = []
        insert_data = {
            "distribution_id": distribution_id,
            "name": name,
            "value": value,
            "snp_value": snp_value if snp_value else []
        }
        data_list.append(insert_data)
        self.col_insert_data("sg_distribution_detail", data_list, "false")

    def sg_scatter(self, task_id, origin_id, name, types='1', location=None, categories=None):
        """
        散点图主表
        :param task_id:
        :param origin_id:
        :param name:
        :param types:
        :param categories:
        :param location:
        :return:
        """
        origin_id = self.check_objectid(origin_id)
        insert_data = {
            "task_id": task_id,
            "origin_id": origin_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "name": name,
            "type": types,
            "location": location if location else "",
            "categories": categories if categories else []
        }
        main_id = self.db["sg_scatter"].insert_one(insert_data).inserted_id
        return main_id

    def sg_scatter_detail(self, scatter_id, name, value, x_max=None, y_max=None):
        """
        散点图细节表
        :param scatter_id:
        :param name:
        :param value:
        :param x_max:
        :param y_max:
        :return:
        """
        scatter_id = self.check_objectid(scatter_id)
        data_list = []
        insert_data = {
            "scatter_id": scatter_id,
            "name": name,
            "value": value,
            "x_max": x_max,
            "y_max": y_max
        }
        data_list.append(insert_data)
        self.col_insert_data("sg_scatter_detail", data_list, "false")

    def sg_collinearity(self, task_id, origin_id, name, left_categorie, right_categorie, data_type, location=None):
        """
        图谱质量评估-共线性分析图2主表
        :param task_id:
        :param origin_id:
        :param name:
        :param right_categorie:
        :param left_categorie:
        :param location:
        :param data_type:
        :return:
        """
        origin_id = self.check_objectid(origin_id)
        insert_data = {
            "task_id": task_id,
            "origin_id": origin_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "name": name,
            "data_type": data_type,
            "type": 1,
            "location": location if location else "",
            "left_categorie": left_categorie,
            "right_categorie": right_categorie
        }
        main_id = self.db["sg_collinearity"].insert_one(insert_data).inserted_id
        return main_id

    def sg_collinearity_detail(self, collinearity_id, name, left_value, right_value):
        """
        图谱质量评估-共线性分析图2细节表
        :param collinearity_id:
        :param name:
        :param left_value:
        :param right_value:
        :return:
        """
        collinearity_id = self.check_objectid(collinearity_id)
        data_list = []
        insert_data = {
            "collinearity_id": collinearity_id,
            "name": name,
            "right_value": right_value,
            "left_value": left_value
        }
        data_list.append(insert_data)
        self.col_insert_data("sg_collinearity_detail", data_list, "false")

    def add_sg_heatmap(self, task_id, origin_id, name, types=1, location=None, categories=None):
        """
        增加heatmap主表
        :param task_id:
        :param origin_id:
        :param name:
        :param types:
        :param location:
        :param categories:
        :return:
        """
        origin_id = self.check_objectid(origin_id)
        insert_data = {
            "task_id": task_id,
            "origin_id": origin_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "name": name,
            "type": types,
            "location": location if location else "",
            "categories": categories if categories else []
        }
        main_id = self.db["sg_heatmap"].insert_one(insert_data).inserted_id
        return main_id

    def add_sg_heatmap_detail(self, heatmap_id, name, value, is_show_log="false", x_max=None):
        """
        添加heatmap的细节表
        :param heatmap_id:
        :param name:
        :param value:
        :param is_show_log:
        :param x_max:  存储染色体1中所有的markers个数，只是针对遗传图谱中的数据
        :return:
        """
        heatmap_id = self.check_objectid(heatmap_id)
        data_list = []
        insert_data = {
            "heatmap_id": heatmap_id,
            "name": name,
            "value": value,
            "x_max": x_max,
        }
        data_list.append(insert_data)
        self.col_insert_data("sg_heatmap_detail", data_list, is_show_log)

    def data_split(self, list_data, step=400000):
        """
        当一条记录中的，某个字段的值大于16M的时候，我们要进行切割后重组后导表
        :param list_data:
        :param step: 按照400000一切割
        :return:
        """
        new_data = []
        if len(list_data) < step:
            new_data.append(list_data)
        else:
            for i in range(0, len(list_data), step):
                temp = list_data[i: i + step]
                new_data.append(temp)
        return new_data

    def check_objectid(self, id_):
        """
        用于检查并转成成ObjectID
        :param id_:
        :return:
        """
        if not isinstance(id_, ObjectId):
            if isinstance(id_, StringTypes):
                id_ = ObjectId(id_)
            else:
                raise Exception("id必须为ObjectId对象或其对应的字符串!")
        return id_

    def check_exists(self, file_path):
        """
        用于检查文件及文件夹是否存在
        :param file_path:
        :return:
        """
        if not os.path.exists(file_path):
            raise Exception("文件或文件夹{}不存在！".format(file_path))

    def add_sg_tree(self, task_id, origin_id, name, types=1, location=None):
        """
        树的主表
        :param task_id:
        :param origin_id:
        :param name:
        :param types:
        :param location:
        :return:
        """
        origin_id = self.check_objectid(origin_id)
        insert_data = {
            "task_id": task_id,
            "origin_id": origin_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "name": name,
            "type": types,
            "location": location if location else ""
        }
        main_id = self.db["sg_tree"].insert_one(insert_data).inserted_id
        return main_id

    def add_sg_tree_detail(self, tree_id, name, value, is_show_log="false"):
        """
        添加树的细节表
        :param tree_id:
        :param name:
        :param value:
        :param is_show_log:
        :return:
        """
        tree_id = self.check_objectid(tree_id)
        data_list = []
        insert_data = {
            "tree_id": tree_id,
            "name": name,
            "value": value
        }
        data_list.append(insert_data)
        self.col_insert_data("sg_tree_detail", data_list, is_show_log)

    def get_file_len(self, file_path):
        """
        获取文件行数
        :return:
        """
        return len(open(file_path, 'rU').readlines())

# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' @20191230
import os
import re
import datetime
import json
from bson.son import SON
from bson.objectid import ObjectId
import types
import gridfs
from biocluster.api.database.base import Base, report_check
# from biocluster.config import Config


class Picrust2Predict(Base):
    def __init__(self, bind_object):
        super(Picrust2Predict, self).__init__(bind_object)
        self._project_type = 'meta'

    @report_check
    def add_function_prediction(self, name=None, params=None, otu_id=0):
        """
        导入MongoDB的主表
        :param name: 主表name的名称
        :param params: params是页面参数的名称的dict
        :param otu_id: otu表的id
        :return:
        """
        if otu_id != 0 and not isinstance(otu_id, ObjectId):
            if isinstance(otu_id, types.StringTypes):
                otu_id = ObjectId(otu_id)
            else:
                self.bind_object.set_error("otu_id必须为ObjectId对象或其对应的字符串!", code="51008901")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": otu_id})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(otu_id)))
            self.bind_object.set_error("sg_otu表中找不到相应记录", code="51008902")
        project_sn = result['project_sn']
        task_id = result['task_id']
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "name": name if name else "PICRUSt2_Origin",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "status": "end",
            "otu_id": otu_id,
            "desc": "PICRUSt2功能预测主表",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db["sg_picrust2"]
        prediction_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({"_id": prediction_id}, {'$set': {'main_id': prediction_id}})
        return prediction_id

    @report_check
    def update_specimen(self, sample_path, prediction_id):
        """
        用于更新主表id
        :param sample_path: enzyme样本名称表
        :param prediction_id: 主表id
        :return:
        """
        if not isinstance(prediction_id, ObjectId):
            if isinstance(prediction_id, types.StringTypes):
                prediction_id = ObjectId(prediction_id)
            else:
                self.bind_object.set_error("prediction_id必须为ObjectId对象或其对应的字符串！", code="51008903")
        collection = self.db["sg_picrust2"]
        result = collection.find_one({"_id": prediction_id})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_picrust2表里找到相应的记录".format(str(prediction_id)))
            self.bind_object.set_error("在sg_picrust2找不到相应的记录", code="51008904")
        specimen = []
        with open(sample_path, "rb") as t:
            line1 = t.readline().strip().split('\t')
            for i in range(2, len(line1)):
                specimen.append(line1[i])
        collection.update_one({'_id': ObjectId(prediction_id)}, {'$set': {'specimen': specimen}})

    @report_check
    def add_cog_annotation(self, prediction_id, annotation_path, type):
        """
        导cog表详情表
        :param prediction_id: 主表id
        :param function_path: 注释需要导入MongoDB的表（cog和function）
        :param type:（cog和function）
        :return:
        """
        if not isinstance(prediction_id, ObjectId):
            if isinstance(prediction_id, types.StringTypes):
                prediction_id = ObjectId(prediction_id)
            else:
                self.bind_object.set_error("prediction_id必须为ObjectId对象或其对应的字符串！", code="51008905")
        if not os.path.exists(annotation_path):
            self.bind_object.set_error("%s所指定的路径不存在，请检查！", code="51008906")
        if type not in ["cog", "function"]:
            self.bind_object.set_error("所指定的type类型不对，请检查！", code="51008907")
        data_list = []
        with open(annotation_path, "rb") as t:
            lines = t.readlines()
            head_list = lines[0].strip().split("\t")
            for line in lines[1:]:
                line = line.strip().split("\t")
                if type in ['function']:
                    data = [
                        ("prediction_id", prediction_id),
                        ("category", line[0]),
                        ("description", line[1]),
                        ("type", type),
                        ]
                    mylist = []
                    summary = 0
                    for i in range(2, len(line)):
                        data += [
                            (head_list[i], float(line[i]))
                        ]
                        summary += float(line[i])
                        mylist.append(float(line[i]))
                    box = self.get_box_new(mylist)
                    data.append(("box", box))
                    data.append(("calculation", summary))
                    data = SON(data)
                    data_list.append(data)
                elif type in ['cog']:
                    data = [
                        ("prediction_id", prediction_id),
                        ("cog_id", line[0]),
                        ("description", line[1]),
                        ("category", line[2]),
                        ("type", type),
                        ]
                    mylist = []
                    calculation, summary = 0, 0
                    for i in range(3, len(line)):
                        data += [
                            (head_list[i], float(line[i]))
                        ]
                        summary += float(line[i])
                        mylist.append(float(line[i]))
                    box = self.get_box_new(mylist)
                    data.append(("box", box))
                    data.append(("calculation", summary))
                    data = SON(data)
                    data_list.append(data)
        try:
            collection = self.db["sg_picrust2_cog_abu"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_picrust2"]
            main_collection.update_one({"_id": prediction_id}, {'$set': {'main_id': prediction_id}})
        except:
            self.bind_object.logger.error("导入cog功能预测信息出错")
        else:
            self.bind_object.logger.info("导入cog功能预测信息成功!")

    @report_check
    def add_kegg_abu(self, prediction_id, abu_path, type):
        """
        将KO、module、enzyme的丰度表导入MongoDB
        :param prediction_id: 主表的id或者main_id
        :param abu_path: 丰度表的path
        :param type: 类型ko、module、enzyme
        :return:
        """
        if not isinstance(prediction_id, ObjectId):
            if isinstance(prediction_id, types.StringTypes):
                prediction_id = ObjectId(prediction_id)
            else:
                self.bind_object.set_error("prediction_id必须为ObjectId对象或其对应的字符串！", code="51008908")

        if not os.path.exists(abu_path):
            self.bind_object.set_error("所指定的路径不存在，请检查！", code="51008909")
        self.bind_object.logger.info(abu_path)

        data_list = []

        rank_list = []
        with open(abu_path, "rb") as m:
            lins = m.readlines()
            head = lins[0].strip().split("\t")
            for lin in lins[1:]:
                lin = lin.strip().split("\t")
                summary = 0
                for i in range(2, len(head)):
                    summary += float(lin[i])
                rank_list.append((lin[0],summary))
        rank_list.sort(key=self.takeSecond, reverse=True)

        with open(abu_path, "rb") as en:
            lines = en.readlines()
            head_list = lines[0].strip().split("\t")
            for line in lines[1:]:
                line = line.strip().split("\t")
                if type in ["ko", "module", "enzyme"]:
                    data = [
                        ("prediction_id", prediction_id),
                        ("name", line[0]),
                        ("definition", line[1]),
                        ("type", type),
                        ]
                    sum_line = 0
                    for i in range(2, len(head_list)):
                        data += [
                            (head_list[i], float(line[i]))
                        ]
                        sum_line += float(line[i])
                    index = int(rank_list.index((line[0],sum_line))) + 1
                    data.append(("rank", index))
                    data = SON(data)
                    data_list.append(data)
                else:
                    self.bind_object.set_error("传入不正确的type类型", code="51008910")
        try:
            collection = self.db["sg_picrust2_kegg_abu"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_picrust2"]
            main_collection.update_one({"_id": prediction_id}, {'$set': {'main_id': prediction_id}})
        except Exception as e:
            self.bind_object.logger.error("导入kegg注释功能预测信息出错！" )
            self.bind_object.set_error("导入kegg注释功能预测信息出错！", code="51008911")
        else:
            self.bind_object.logger.info("导入kegg注释功能预测信息成功!")

    @report_check
    def add_kegg_level(self, prediction_id, kegg_path):
        if not isinstance(prediction_id, ObjectId):
            if isinstance(prediction_id, types.StringTypes):
                prediction_id = ObjectId(prediction_id)
            else:
                self.bind_object.set_error("prediction_id必须为ObjectId对象或其对应的字符串！", code="51008912")

        if not os.path.exists(kegg_path):
            self.bind_object.set_error("所指定的路径不存在，请检查！", code="51008913")
        else:
            if os.path.isdir(kegg_path):
                self.bind_object.logger.info(kegg_path)
                self.bind_object.logger.info("正确的pathway文件夹")
            else:
                self.bind_object.set_error("pathway文件夹不存在，请检查！", code="51008914")

        data_list = []

        rank_data = []
        for j in [1, 2, 3]:
            rank_list = []
            level_path = kegg_path + "/prediction_pathway.L" + str(j) + ".xls"
            with open(level_path, "rb") as m:
                lins = m.readlines()
                head = lins[0].strip().split("\t")
                for lin in lins[1:]:
                    lin = lin.strip().split("\t")
                    summary = 0
                    if j == 1:
                        add_num = 0
                    elif j == 2:
                        add_num = 1
                    else:
                        add_num = 3
                    for i in range(1+add_num, len(head)):
                        summary += float(lin[i])
                    rank_list.append((lin[0],summary))
            rank_list.sort(key=self.takeSecond, reverse=True)
            rank_data.append(rank_list)
        for j in [1, 2, 3]:
            level_path = kegg_path + "/prediction_pathway.L" + str(j) + ".xls"
            if os.path.exists(level_path):
                with open(level_path, "rb") as le:
                    lines = le.readlines()
                    head_list = lines[0].strip().split('\t')

                    for line in lines[1:]:
                        line = line.strip().split("\t")
                        data = [
                            ("prediction_id", prediction_id),
                            ("level", j),
                            ("name", line[0]),
                        ]
                        if j==2:
                            data.append(("type", "level2"))
                            data.append(('level1', line[1]))
                            add_num = 1
                            sum_line = 0
                            for i in range(1+add_num, len(head_list)):
                                data += [
                                    (head_list[i], float(line[i]))
                                ]
                                sum_line += float(line[i])
                            list_rank = rank_data[1]
                            index = int(list_rank.index((line[0],sum_line))) + 1
                            data.append(("rank", index))
                        elif j==3:
                            data.append(("type", "level3"))
                            data.append(('level2', line[2]))
                            data.append(('level1',line[3]))
                            data.append(('pathway', line[1]))
                            add_num = 3
                            sum_line = 0
                            for i in range(1+add_num, len(head_list)):
                                data += [
                                    (head_list[i], float(line[i]))
                                ]
                                sum_line += float(line[i])
                            list_rank = rank_data[2]
                            index = int(list_rank.index((line[0],sum_line))) + 1
                            data.append(("rank", index))
                        else:
                            data.append(("type", "level1"))
                            add_num = 0
                            sum_line = 0
                            for i in range(1+add_num, len(head_list)):
                                data += [
                                    (head_list[i], float(line[i]))
                                ]
                                sum_line += float(line[i])
                            list_rank = rank_data[0]
                            index = int(list_rank.index((line[0],sum_line))) + 1
                            data.append(("rank", index))

                        data = SON(data)
                        data_list.append(data)
                    self.update_main_type(prediction_id,"level%s"%j)
            else:
                self.bind_object.set_error("kegg功能预测的结果文件中level丰度表不存在，请检查丰度表！", code="51008915")
        try:
            collection = self.db["sg_picrust2_kegg_pathway"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_picrust2"]
            main_collection.update_one({"_id": prediction_id}, {'$set': {'main_id': prediction_id}})
        except:
            self.bind_object.logger.error("导入kegg的pathway预测信息出错")
            self.bind_object.set_error("导入kegg的pathway功能预测出错", code="51008916")
        else:
            self.bind_object.logger.info("导入kegg的pathway功能预测信息成功!")

    @report_check
    def add_metacyc_abu(self, prediction_id, abu_path):
        """
        导入metacyc详情表结果
        :param prediction_id: 主表id
        :param abu_path: 导入表格的路径
        :return:
        """
        if not isinstance(prediction_id, ObjectId):
            if isinstance(prediction_id, types.StringTypes):
                prediction_id = ObjectId(prediction_id)
            else:
                self.bind_object.set_error("prediction_id必须为ObjectId对象或其对应的字符串！", code="51008917")

        if not os.path.exists(abu_path):
            self.bind_object.set_error("所指定的路径不存在，请检查！", code="51008918")

        data_list = []

        rank_list = []
        with open(abu_path, "rb") as m:
            lins = m.readlines()
            head = lins[0].strip().split("\t")
            for lin in lins[1:]:
                lin = lin.strip().split("\t")
                summary = 0
                for i in range(2, len(head)):
                    summary += float(lin[i])
                rank_list.append((lin[0],summary))
        rank_list.sort(key=self.takeSecond, reverse=True)
        with open(abu_path, "rb") as en:
            lines = en.readlines()
            head_list = lines[0].strip().split("\t")
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = [
                    ("prediction_id", prediction_id),
                    ("name", line[0]),
                    ("desc", line[1]),
                ]
                sum_line = 0
                for i in range(2, len(head_list)):
                    data += [
                        (head_list[i], float(line[i]))
                    ]
                    sum_line += float(line[i])
                index = int(rank_list.index((line[0],sum_line))) + 1
                data.append(("rank", index))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_picrust2_metacyc_abu"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_picrust2"]
            main_collection.update_one({"_id": prediction_id}, {'$set': {'main_id': prediction_id}})
        except Exception as e:
            self.bind_object.logger.error("导入metacyc功能预测信息出错！" )
            self.bind_object.set_error("导入metacyc功能预测信息出错！", code="51008919")
        else:
            self.bind_object.logger.info("导入metacyc功能预测信息成功!")

    def get_box_new(self, mylist):
        """
        目的是为了获得箱型图的画图的点和捕获异常点
        :param mylist:
        :return:
        """
        mylist = sorted(mylist)
        length = len(mylist)
        filter_list = []
        out = []
        if length > 3:
            q1_index = 0.25 * (length + 1)
            q2_index = 0.5 * (length + 1)
            q3_index = 0.75 * (length + 1)
            q1_index_int = int(q1_index)
            q2_index_int = int(q2_index)
            q3_index_int = int(q3_index)
            q1 = mylist[q1_index_int - 1] + (mylist[q1_index_int] - mylist[q1_index_int - 1]) * (q1_index - q1_index_int)
            q3 = mylist[q3_index_int - 1] + (mylist[q3_index_int] - mylist[q3_index_int - 1]) * (q3_index - q3_index_int)
            q2 = mylist[q2_index_int - 1] + (mylist[q2_index_int] - mylist[q2_index_int - 1]) * (q2_index - q2_index_int)
            qd = q3 - q1
            max_limit = 1.5 * qd + q3
            min_limit = q1 - 1.5 * qd
            new_list = []
            for i in mylist:
                if i >= min_limit and i <= max_limit:
                    new_list.append(i)
                else:
                    filter_list.append(i)
            max_box = new_list[-1]
            min_box = new_list[0]
        elif(length == 3):
            max_box = mylist[2]
            min_box = mylist[0]
            q3 = float(mylist[1] + mylist[2]) / 2
            q2 = mylist[1]
            q1 = float(mylist[1] + mylist[0]) / 2
        elif(length == 2):
            max_box = mylist[1]
            min_box = mylist[0]
            q2 = float(mylist[1] + mylist[0]) / 2
            q3 = (mylist[1] + q2) / 2
            q1 = (q2 + mylist[0]) / 2
        elif(length == 1):
            max_box = mylist[0]
            min_box = mylist[0]
            q2 = mylist[0]
            q3 = mylist[0]
            q1 = mylist[0]
        else:
            self.bind_object.set_error("cog数据太少，无法画箱线图", code="51008920")
        out.append(min_box)
        out.append(q1)
        out.append(q2)
        out.append(q3)
        out.append(max_box)
        out.append(filter_list)
        return out

    def takeSecond(self, elem):
        return elem[1]

    @report_check
    def update_main_type(self, main_id, type):
        """
        用于更新主表的type字段，供前端进行筛选
        :param main_id: 主表id
        :param type: type类型
        :return:
        """
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                main_id = main_id
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串！", code="51008921")
        collection = self.db["sg_picrust2"]
        try:
            table = collection.find_one({"_id": main_id})
            find_table_type = collection.find({"_id": main_id, "table_type": {"$exists": True}}).count()
            self.bind_object.logger.info(find_table_type)
            if int(find_table_type) == 1:
                table_type = table["table_type"]
                table_type_list = table_type.split("|")
                table_kegg = table_type_list[0]
                table_kegg_list = table_kegg.split(",")
                if type in ["ko", "module", "enzyme", "level1", "level2", "level3"]:
                    new_table_type_list = []
                    table_kegg_list.append(type)
                    new_table_kegg = ",".join(table_kegg_list)
                    new_table_type_list.append(new_table_kegg)
                    if len(table_type_list) == 2:
                        new_table_type_list.append(table_type_list[1])
                    new_table_type = "|".join(new_table_type_list)
                elif type in ['metacyc']:
                    new_table_type_list = []
                    new_table_type_list.append(table_type_list[0])
                    new_table_type_list.append(type)
                    new_table_type = "|".join(new_table_type_list)
                else:
                    self.bind_object.set_error("不存在的type类型", code="51008922")
                self.bind_object.logger.info(new_table_type)
                collection.update_one({'_id': ObjectId(main_id)}, {'$set': {'table_type': new_table_type}})
            else:
                table_type_list = []
                table_type_list.append(type)
                new_table_type = "|".join(table_type_list)
                self.bind_object.logger.info(new_table_type)
                collection.update_one({'_id': ObjectId(main_id)},{'$set': {'table_type': new_table_type}})
        except Exception as e:
            self.bind_object.logger.error("更新主表table_type出错！" )
            self.bind_object.set_error("更新主表table_type出错！", code="51008923")

    def get_database(self, task_id):
        """
        根据task_id获取sg_task表中的database
        :param task_id:
        :return:
        """
        collection = self.db["sg_task"]
        task_info = collection.find_one({"task_id": task_id})
        if task_info:
            return task_info
        else:
            self.bind_object.logger.error("任务不存在:{}".format(str(task_id)))

    def add_path(self, main_id, path):
        """
        添加结果路径到主表
        :param task_id:
        :return:
        """
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                prediction_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("prediction_id必须为ObjectId对象或其对应的字符串！", code="51008903")
        collection = self.db["sg_picrust2"]
        collection.update_one({'_id': ObjectId(main_id)}, {'$set': {'path': path}})
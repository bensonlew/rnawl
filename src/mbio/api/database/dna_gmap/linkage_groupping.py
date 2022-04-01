# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# modified 2018.06.26

from biocluster.api.database.base import Base, report_check
from api_base import ApiBase
import copy
from bson.objectid import ObjectId
from types import StringTypes
import datetime
import json
import os
import re


class LinkageGroupping(ApiBase):
    def __init__(self, bind_object):
        """
        用于连锁分群中的所有的导表
        """
        super(LinkageGroupping, self).__init__(bind_object)
        self._project_type = "dna_gmap"

    def add_sg_lg(self, task_id, project_sn, params):
        """
        sg_lg
        add by hongdong@20180705
        :return:
        """
        main_id = self.add_main_table("sg_lg", task_id, project_sn, params, "origin_linkage_grouping",
                                      "开始进行图谱评估")
        self.update_db_record("sg_lg", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_evalutaion(self, sg_lg_id):
        """
        图谱评估的主表
        add by hongdong@20180705
        :return:
        """
        lg_id = self.check_objectid(sg_lg_id)
        result = self.col_find_one("sg_lg", {"_id": lg_id})
        try:
            name = result['name']
            task_id = result['task_id']
            project_sn = result['project_sn']
        except:
            raise Exception("sg_lg中没有找到name or task_id or project_id")
        main_id = self.add_main_table("sg_evalutaion", task_id, project_sn, "origin_params", name, "开始进行图谱评估")
        self.update_db_record("sg_evalutaion", {"_id": main_id}, {"origin_id": lg_id, "main_id": main_id})
        return main_id

    def add_sg_evalutaion_stat(self, file_path, evalutaion_id, data_type="total"):
        """
        图谱评估-基本信息统计表
        total.mapstat, male.mapstat, female.mapstat, sexAver.mapstat
        add by hongdong@20180705
        :return:
        """
        evalutaion_id = self.check_objectid(evalutaion_id)
        self.check_exists(file_path)
        data_list = []
        with open(file_path, 'r') as r:
            data = r.readlines()[1:]
            for m in data:
                temp = m.strip().split('\t')
                insert_data = {
                    "lg": temp[0],
                    "marker_num": int(temp[1]),
                    "unique_num": int(temp[2]),
                    "map_distance": float("%.3f" % float(temp[3])),
                    "average_distance": float(temp[4]),
                    "gap_less_5cm": float("%.3f" % float(temp[5])),
                    "max_gap": float("%.3f" % float(temp[6])),
                    "evalutaion_id": evalutaion_id,
                    "data_type": data_type
                }
                data_list.append(insert_data)
        if len(data_list) != 0:
            self.col_insert_data("sg_evalutaion_stat", data_list)
        else:
            self.bind_object.logger.info("文件{}为空，不进行导表！".format(file_path))

    def add_yichuantupu(self, task_id, evalutaion_id, file_path, data_type, updata_chr_list=False):
        """
        谱图评估中的遗传图谱图数据, 当是f1群体的时候，要导入雌性，中性，雄性，在workflow中循环到3次，name,分别为female，male，
        sexaver，如果非f1群体，就是导入total,name为total
        nocp:
        total.map
        cp:
        total.female.map， total.male.map ， total.sexAver.map
        add by hongdong@20180705
        :param task_id:
        :param file_path:
        :param evalutaion_id: 连锁分群的主表 sg_lg
        :param data_type:female，male,sexaver, total
        :param updata_chr_list:True or False
        :return:
        """
        evalutaion_id = self.check_objectid(evalutaion_id)
        self.check_exists(file_path)
        lds = []
        chrs_ = []
        scas_ = []
        all_value = []
        dis_min = 10000000000000000
        dis_max = 0
        with open(file_path, 'r') as r:
            for line in r:
                temp = line.strip().split('\t')
                if temp[0] == 'group':
                    lds.append(int(temp[1]))
                    all_value.append({temp[0]: int(temp[1])})
                else:
                    temp_ = temp[0].split("_")[0]
                    if re.match(r"^chr.*", temp_.lower()):
                        if temp_ not in chrs_:
                            try:
                                chrs_.append(int(temp_[3:]))
                            except:
                                chrs_.append(temp_[3:])
                    elif re.match(r'^sca.*', temp_.lower()):
                        if temp_ not in scas_:
                            try:
                                scas_.append(int(temp_[3:]))
                            except:
                                scas_.append(temp_[3:])
                    if float("%.4f" % float(temp[1])) > dis_max:
                        dis_max = float("%.4f" % float(temp[1]))
                    if float("%.4f" % float(temp[1])) < dis_min:
                        dis_min = float("%.4f" % float(temp[1]))
                    all_value.append({temp[0]: temp[1]})
        # print "lds:", lds
        cat_lds = list(set(lds))
        cat_lds.sort()
        # print "cat_lds:", cat_lds
        chrs_.sort()
        scas_.sort()
        chrs = list(set(chrs_))
        chr_sca = self.set_chr_sca(chrs, [])
        main_id = self.sg_distribution(task_id, evalutaion_id, data_type, "2", dis_min, dis_max,
                                       "evalutaion_distribution", ["LG"+str(i) for i in cat_lds])
        value_dict = {}
        for m in range(0, len(lds)):
            value = []
            tem_value = []
            if m == len(lds) - 1:
                if {"group": lds[m]} in all_value:
                    tem_value = all_value[all_value.index({"group": lds[m]}) + 1:]
            else:
                if {"group": lds[m]} in all_value and {"group": lds[m + 1]} in all_value:
                    tem_value = all_value[all_value.index({"group": lds[m]}) + 1: all_value.index(
                        {"group": lds[m + 1]})]
            if tem_value:
                for n in tem_value:
                    value.append([n.keys()[0].split("_")[0], float("%.4f" % float(n.values()[0]))])
                # print "value:", value
            value_dict[lds[m]] = value
            m += 1
        # print "----------------------"
        # print "value_dict:", value_dict
        # print "----------------------"
        for n in cat_lds:
            if n in value_dict.keys():
                self.sg_distribution_detail(main_id, 'LG{}'.format(n), value_dict[n])
        print "yichuantupu ok"
        if updata_chr_list:
            self.update_db_record("sg_evalutaion", {"_id": evalutaion_id}, {"chr_list": chr_sca,
                                                                            "lgs": ["LG"+str(i) for i in cat_lds]})

    def add_yichuantupu_old(self, task_id, evalutaion_id, file_path, data_type, updata_chr_list=False):
        """
        谱图评估中的遗传图谱图数据, 当是f1群体的时候，要导入雌性，中性，雄性，在workflow中循环到3次，name,分别为female，male，
        sexaver，如果非f1群体，就是导入total,name为total
        nocp:
        total.map
        cp:
        total.female.map， total.male.map ， total.sexAver.map
        add by hongdong@20180705
        :param task_id:
        :param file_path:
        :param evalutaion_id: 连锁分群的主表 sg_lg
        :param data_type:female，male,sexaver, total
        :param updata_chr_list:True or False
        :return:
        """
        evalutaion_id = self.check_objectid(evalutaion_id)
        self.check_exists(file_path)
        chrs_ = []
        scas_ = []
        all_value = []
        dis_min = 10000000000000000
        dis_max = 0
        with open(file_path, 'r') as r:
            for line in r:
                if re.match(r'^group.*', line):
                    pass
                else:
                    temp = line.strip().split('\t')
                    temp_ = temp[0].split("_")[0]
                    if re.match(r"^chr.*", temp_):
                        if temp_ not in chrs_:
                            chrs_.append(int(temp_[3:]))
                        all_value.append({temp_: temp[1]})
                    elif re.match(r'^sca.*', temp_):
                        if temp_ not in scas_:
                            scas_.append(int(temp_[3:]))
                        all_value.append({temp_: temp[1]})
                    if float("%.4f" % float(temp[1])) > dis_max:
                        dis_max = float("%.4f" % float(temp[1]))
                    if float("%.4f" % float(temp[1])) < dis_min:
                        dis_min = float("%.4f" % float(temp[1]))
        # print "dis_max:{}, dis_min:{}".format(dis_max, dis_min)
        chrs_.sort()
        scas_.sort()
        chrs = list(set(chrs_))
        scas = list(set(scas_))
        chr_sca = self.set_chr_sca(chrs, scas)
        main_id = self.sg_distribution(task_id, evalutaion_id, data_type, "2", dis_min, dis_max,
                                       "evalutaion_distribution", chr_sca)
        for m in chr_sca:
            value = []
            for n in all_value:
                if n.keys()[0] == m:
                    value.append(float("%.4f" % float(n.values()[0])))
            self.sg_distribution_detail(main_id, m, value)
        print "yichuantupu ok"
        if updata_chr_list:
            self.update_db_record("sg_evalutaion", {"_id": evalutaion_id}, {"chr_list": chr_sca})

    def set_chr_sca(self, chrs_, scas_):
        """
        合并chr与sca，完成排序
        add by hongdong@20180705
        :return:
        """
        chr_sca = []
        for m in chrs_:
            chr_sca.append("chr{}".format(m))
        if scas_:
            for n in scas_:
                chr_sca.append("sca{}".format(n))
        return chr_sca

    def collinearity_old(self, task_id, file_path, evalutaion_id, spearman_path, data_type):
        """
        遗传图谱中所有共线性分析图的--导表的逻辑有变化，不在按照染色体去进行展示，按照group后面的值连锁群去展示
        所有染色体的共线性图
        单个染色体的共线性图
        共线性分析表
        add by hongdong@20180705
        :param task_id:
        :param file_path: total.female.map， total.male.map ， total.sexAver.map
        :param evalutaion_id:
        :param spearman_path:
        :param data_type: female，male,sexaver, total
        :return:
        """
        evalutaion_id = self.check_objectid(evalutaion_id)
        lds = []
        all_value = []
        chrs_ = []
        scas_ = []
        pos_max = 0
        value_max = 0
        chr_value = []
        sca_value = []
        with open(file_path, 'r') as r:
            for line in r:
                temp = line.strip().split('\t')
                if temp[0] == "group":
                    lds.append(temp[1])
                    all_value.append({temp[0]: temp[1]})
                    try:
                        pos = temp[0].split("_")[1].split("-")[0]
                    except:
                        if len(temp[0].split("_")) == 3 and temp[0].split("_")[1] == 'random':
                            pos = temp[0].split("_")[2].split("-")[0]
                        else:
                            continue
                    chr_ = temp[0].split("_")[0]
                    if int(pos) > pos_max:
                        pos_max = int(pos)
                    if float(temp[1]) > value_max:
                        value_max = float(temp[1])
                    if re.match(r"^chr.*", chr_.lower()):
                        if chr_ not in chrs_:
                            try:
                                chrs_.append(int(chr_[3:]))
                            except:
                                chrs_.append(chr_[3:])
                        chr_value.append({temp[0]: temp[1]})
                    elif re.match(r"sca.*", chr_.lower()):
                        if chr_ not in scas_:
                            try:
                                scas_.append(int(chr_[3:]))
                            except:
                                scas_.append(chr_[3:])
                        sca_value.append({temp[0]: temp[1]})
        if chrs_:
            all_value = chr_value
            scas_ = []
        else:
            all_value = sca_value
        chrs_.sort()
        scas_.sort()
        chrs = list(set(chrs_))
        scas = list(set(scas_))
        chr_sca = self.set_chr_sca(chrs, scas)
        scatter_id = self.sg_scatter(task_id, evalutaion_id, data_type, types='2', location="evalutaion_scatter",
                                     categories=chr_sca)
        for m in chr_sca:
            chr_value = []
            left_categorie = []   # pos
            right_categorie = []   # 相关系数
            left_value = []
            right_value = []
            chr_value_x = []
            chr_value_y = []
            for n in all_value:
                if n.keys()[0].split("_")[0] == m:
                    chr_value_x.append(float(n.values()[0]))
                    chr_value_y.append(int(n.keys()[0].split("_")[-1].split("-")[0]))
                    chr_value.append([round(float(n.values()[0]), 2), int(n.keys()[0].split("_")[-1].split("-")[0])])
                    left_categorie.append(n.keys()[0].split("_")[-1].split("-")[0])
                    left_value.append(float(n.keys()[0].split("_")[-1].split("-")[0]))
                    right_categorie.append(n.values()[0])
                    right_value.append(float(n.values()[0]))
            x_max_ = self.get_max(chr_value_x)
            y_max_ = self.get_max(chr_value_y)
            self.sg_scatter_detail(scatter_id, m, chr_value, x_max_, y_max_)
            right_value_max = self.get_max(right_value)
            left_value_max = self.get_max(left_value)
            collinearity_id = self.sg_collinearity(task_id, evalutaion_id, m, left_categorie, right_categorie,
                                                   data_type, "collinearity_single_chr")
            self.sg_collinearity_detail(collinearity_id, m, self.divide(left_value, left_value_max),
                                        self.divide(right_value, right_value_max))
        self.add_sg_collinearity_stat(evalutaion_id, spearman_path, data_type)
        print "相关性分析结果{} ok".format(data_type)

    def collinearity(self, task_id, file_path, evalutaion_id, spearman_path, data_type):
        """
        遗传图谱中所有共线性分析图的 --暂时弃用
        所有染色体的共线性图
        单个染色体的共线性图
        共线性分析表
        add by hongdong@20180705
        :param task_id:
        :param file_path: total.female.map， total.male.map ， total.sexAver.map
        :param evalutaion_id:
        :param spearman_path:
        :param data_type: female，male,sexaver, total
        :return:
        """
        evalutaion_id = self.check_objectid(evalutaion_id)
        chrs_ = []
        scas_ = []
        pos_max = 0
        value_max = 0
        chr_value = []
        sca_value = []
        with open(file_path, 'r') as r:
            for line in r:
                if re.match(r'^group.*', line):
                    pass
                else:
                    temp = line.strip().split('\t')
                    pos = temp[0].split("_")[-1].split("-")[0]
                    chr_ = temp[0].split("_")[0]
                    if int(pos) > pos_max:
                        pos_max = int(pos)
                    if float(temp[1]) > value_max:
                        value_max = float(temp[1])
                    if re.match(r"^chr.*", chr_.lower()):
                        if chr_ not in chrs_:
                            try:
                                chrs_.append(int(chr_[3:]))
                            except:
                                chrs_.append(chr_[3:])
                        chr_value.append({temp[0]: temp[1]})
                    elif re.match(r"sca.*", chr_.lower()):
                        if chr_ not in scas_:
                            try:
                                scas_.append(int(chr_[3:]))
                            except:
                                scas_.append(chr_[3:])
                        sca_value.append({temp[0]: temp[1]})
        if chrs_:
            all_value = chr_value
            scas_ = []
        else:
            all_value = sca_value
        chrs_.sort()
        scas_.sort()
        chrs = list(set(chrs_))
        scas = list(set(scas_))
        chr_sca = self.set_chr_sca(chrs, scas)
        scatter_id = self.sg_scatter(task_id, evalutaion_id, data_type, types='2', location="evalutaion_scatter",
                                     categories=chr_sca)
        for m in chr_sca:
            chr_value = []
            left_categorie = []  # pos
            right_categorie = []  # 相关系数
            left_value = []
            right_value = []
            chr_value_x = []
            chr_value_y = []
            for n in all_value:
                try:  # 这里兼容chr01导致的匹配失败的bug
                    new_chr_name = n.keys()[0].split("_")[0][0:3] + str(int(n.keys()[0].split("_")[0][3:]))
                except:
                    new_chr_name = n.keys()[0].split("_")[0]
                if n.keys()[0].split("_")[0] == m or new_chr_name == m:
                    chr_value_x.append(float(n.values()[0]))
                    chr_value_y.append(int(n.keys()[0].split("_")[-1].split("-")[0]))
                    chr_value.append([round(float(n.values()[0]), 2), int(n.keys()[0].split("_")[-1].split("-")[0])])
                    left_categorie.append(n.keys()[0].split("_")[-1].split("-")[0])
                    left_value.append(float(n.keys()[0].split("_")[-1].split("-")[0]))
                    right_categorie.append(n.values()[0])
                    right_value.append(float(n.values()[0]))
            if len(chr_value_x) == 0:
                raise Exception("请检查染色体的命名是否正确，chr01必须写成chr1")
            x_max_ = self.get_max(chr_value_x)
            y_max_ = self.get_max(chr_value_y)
            self.sg_scatter_detail(scatter_id, m, chr_value, round(x_max_, 2), y_max_)
            right_value_max = self.get_max(right_value)
            left_value_max = self.get_max(left_value)
            collinearity_id = self.sg_collinearity(task_id, evalutaion_id, m, left_categorie, right_categorie,
                                                   data_type, "collinearity_single_chr")
            self.sg_collinearity_detail(collinearity_id, m, self.divide(left_value, left_value_max),
                                        self.divide(right_value, right_value_max))
        self.add_sg_collinearity_stat(evalutaion_id, spearman_path, data_type)
        print "相关性分析结果{} ok".format(data_type)

    def divide(self, list_, value_):
        """
        列表都除以同一个元素， value_ != 0
        :param list_:
        :param value_:
        :return:
        """
        new_value = []
        if value_ == 0:
            new_value = list_
        else:
            for m in list_:
                new_value.append(m / float(value_))
        return new_value

    def x_y_translation(self, list_, x, y, types="+"):
        """
        遍历列表，然后根据类型，每个元素进行加减
        add by hongdong@20180705
        :param list_:[[1,2], [3,4], [5,6]]
        :param x:
        :param y:
        :param types: “+” or “-”
        :return:
        """
        new_value = []
        for m in list_:
            if types == "+":
                new_value.append([m[0] + x, m[1] + y])
            else:
                new_value.append([m[0] - x, m[1] - y])
        return new_value

    def get_max(self, list_):
        """
        获取列表的最大值
        add by hongdong@20180705
        :param list_:
        :return:
        """
        return max(list_)

    def add_sg_collinearity_stat(self, evalutaion_id, spearman_path, data_type):
        """
        共线性分析表
        add by hongdong@20180705
        :param evalutaion_id:
        :param data_type:
        :param spearman_path::param spearman_path: total.phy.spearman.xls
        :return:
        """
        data_list = []
        with open(spearman_path, "r") as r:
            data = r.readlines()[1:]
            for line in data:
                temp = line.strip().split("\t")
                insert_data = {
                    "evalutaion_id": evalutaion_id,
                    "lg": int(temp[0]),
                    'pearman': abs(float(temp[3])),
                    'data_type': data_type
                }
                data_list.append(insert_data)
        if len(data_list) != 0:
            self.col_insert_data("sg_collinearity_stat", data_list)
        else:
            self.bind_object.logger.info("文件{}为空，不进行导表！".format(spearman_path))

    def add_eva_heatmap(self, file, type_groupping, task_id, origin_id):
        """
        file: 如果是非cp，则文件为total.csv
        如果是cp, 则为total.female.phase, total.male.phase 和 total.sexAver.phase。
        "iii此循环用于建立染色体和坐标对应的字典"
        """
        """
        file: /mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/yichuantupu/total.csv

        """
        starttime = datetime.datetime.now()
        self.check_exists(file)  # 检查文件是否存在
        location = "marker_map_evaluation"
        main_data = {
            "name": type_groupping,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "location": location,
            "type": 2,
            "task_id": task_id,
            "origin_id": origin_id
        }
        main_id1 = self.db["sg_heatmap"].insert_one(main_data).inserted_id
        insert_data_2 = []
        map1 = {}
        pos = []
        self.max1 = []
        with open(file)as f:
            lines = f.readlines()
            sample = None
            category = self.get_char()
            list3 = [1, 12]
            for iii in category:
                map2 = {}
                xmax = 0
                for num in range(len(lines)):
                    row = lines[num].strip().split(",")
                    calc = 0
                    if num == 0:
                        sample = row[3:]  # 储存横坐标样品ID
                        continue
                    for ii in range(len(sample)):
                        color = ii + 3
                        """
                            此处将表格中的内容转成统一格式，对于CP: “A”:1 ,"B":2 "H":3 "-":4
                            对于非CP："1":5 "2":6 "3":7 "4":8 "5":9
                        """
                        if type_groupping in ["total"]:
                            if row[color] not in ["A", "B", "H", "-"]:
                                pass
                                self.bind_object.logger.info("不是正确的CP格式")
                            else:
                                if row[color] == "A":
                                    row[color] = "1"
                                elif row[color] == "B":
                                    row[color] = "2"
                                elif row[color] == "H":
                                    row[color] = "3"
                                elif row[color] == "-":
                                    row[color] = "4"
                        elif type_groupping in ["female", "male", "sexaver"]:
                            if row[color] not in ["1", "2", "0", "-", "5"]:
                                pass
                                # self.bind_object.logger.info("不是正确的非CP格式")
                            else:
                                if row[color] == "1":
                                    row[color] = "5"
                                elif row[color] == "2":
                                    row[color] = "6"
                                elif row[color] == "0":
                                    row[color] = "7"
                                elif row[color] == "-":
                                    row[color] = "8"
                                elif row[color] == "5":
                                    row[color] = "9"
                        else:
                            pass
                            # self.bind_object.logger.info("请输入正确的CP类型")
                        axis = [num, ii, int(row[color])]
                        if row[1] == str(iii):
                            if iii not in map1.keys():
                                map1[iii] = []
                            map1[iii].append(axis)
                length = len(map1[iii])
                split = length / 40000
                if length % 40000 == 0:
                    list1 = [i * 40000 for i in range(1, split + 1)]
                else:
                    list1 = [i * 40000 for i in range(1, split + 1)]
                    list1.append(length)

                for factor in list1:
                    list2 = []
                    for m in range(list1.index(factor) * 40000, factor):
                        list2.append(map1[iii][m])
                    if factor not in map2.keys():
                        map2[factor] = []
                        map2[factor] = list2
                    pos = self.max1[(category.index(iii)) + 1]

                for chr1 in map2.keys():
                    data2 = {
                        "value": map2[chr1],
                        "seg": chr1,
                        "name": iii,
                        "pos": pos,
                        "heatmap_id": main_id1
                    }
                    insert_data_2.append(data2)
        if len(insert_data_2) == 0:
            # self.bind_object.logger.info("{}文件为空！".format(results_path))
            raise Exception("{}文件为空！".format(file))
        else:
            self.col_insert_data("sg_heatmap_detail", insert_data_2)
            self.update_db_record("sg_heatmap", {"_id": main_id1}, {"rows": category, "columns": sample})
        endtime = datetime.datetime.now()
        print (endtime - starttime).seconds

    def get_char(self, file):
        """
        获取染色体个数，同时生成热图中染色体最大值Pos。
        :return:
        """
        chr = []
        sum1 = 0
        with open(file)as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if i == 0:
                    pass
                else:
                    row = lines[i].strip().split(",")
                    if row[1] not in chr:
                        chr.append(row[1])
                        if i == (len(lines) - 1):
                            self.max1.append(sum1)
                            self.max1.append(1)
                        else:
                            self.max1.append(sum1)
                            sum1 = 1
                    else:
                        if i == (len(lines) - 1):
                            sum1 = sum1 + 1
                            self.max1.append(sum1)
                        else:
                            sum1 = sum1 + 1
        return chr

    def add_genetype_heatmap(self, file_path, data_type, task_id, origin_id):
        """
        导入基因型来源热图
         file: 如果是非cp，则文件为total.csv
        如果是cp, 则为total.female.phase, total.male.phase 和 total.sexAver.phase。
        add by hongdong@20180730
        :param file_path:
        :param data_type:
        :param task_id:
        :param origin_id:  # evalutaion_id
        :return:
        """
        origin_id = self.check_objectid(origin_id)
        is_cp = True
        if data_type in ["total"]:
            is_cp = False
        value_data = []
        chrs = []
        start = 3
        columns = []
        with open(file_path, "r") as r:
            for line in r:
                temp = line.strip().split(',')
                if re.match(r'Genotype.*', line):
                    columns = temp[3:]
                else:
                    if is_cp:
                        chrs.append(int(temp[2]))
                        start = 4
                    else:
                        chrs.append(int(temp[1]))
                    value_data.append(temp)
        chrs_new = list(set(chrs))
        chrs_new.sort()
        heatmap_id = self.add_sg_heatmap(task_id, origin_id, data_type, "2", "marker_map_evaluation")
        self.update_db_record("sg_heatmap", {"_id": heatmap_id}, {"rows": chrs_new, "columns": columns})
        for chr_ in chrs_new:
            chr_data = []
            n = 0
            for m in value_data:
                if m[0].split("_")[0] == "chr{}".format(chr_):
                    l = 0
                    for da in m[start:]:
                        color = self.set_color_value(da, is_cp)
                        chr_data.append([n, l, color])
                        l += 1
                    n += 1
            new_list = self.data_split(chr_data, 400000)
            insert_data = []
            x_max = chr_data[-1][0]
            for te in new_list:
                insert_data.append({
                    "heatmap_id": heatmap_id,
                    "name": "chr{}".format(chr_),
                    "value": te,
                    "x_max": x_max
                })
            self.col_insert_data("sg_heatmap_detail", insert_data, "false")

    def set_color_value(self, value, is_cp):
        """
        设定热图的染色
        :param value:
        :param is_cp:
        :return:
        """
        if not is_cp:
            if value == "A":
                new_value = 1
            elif value == "B":
                new_value = 2
            elif value == "H":
                new_value = 3
            elif value == "-":
                new_value = 4
            else:
                new_value = 0
        else:
            if value == "1":
                new_value = 5
            elif value == "2":
                new_value = 6
            elif value == "0":
                new_value = 7
            elif value == "-":
                new_value = 8
            elif value == "5":
                new_value = 9
            else:
                new_value = 0
        return new_value

    def add_sg_evalutaion_detail(self, file_path, evalutaion_id, data_type='total'):
        """
        图谱标记详情
        :param file_path:total.sexAver.info, total.male.info, total.female.info
        total.marker.info
        :param evalutaion_id:
        :param data_type:
        :return:
        """
        data_list = []
        with open(file_path, "r") as r:
            for line in r:
                temp = line.strip().split('\t')
                if re.match(r'^#.*', temp[0]):
                    pass
                else:
                    if len(temp) == 10:
                        insert_data = {
                            "evalutaion_id": evalutaion_id,
                            "marker_id": temp[0],
                            "lg": temp[1],
                            "genetic_pos": float(temp[2]),
                            "chr_id": temp[0].split("_")[0],
                            "physical_pos": int(temp[0].split("_")[-1].split('-')[0]),
                            "type": temp[3],
                            "nind": int(temp[4]),
                            "nmiss": int(temp[5]),
                            "miss_ratio": float("%.2f" % (1 - float(temp[5])/float(temp[4])))*100,
                            "geno": temp[6],
                            "sgeno": temp[7],
                            "signif": float(temp[8]),
                            "segretion": temp[9],
                            "data_type": data_type
                        }
                        data_list.append(insert_data)
        if len(data_list) != 0:
            self.col_insert_data("sg_evalutaion_detail", data_list)
        else:
            self.bind_object.logger.info("文件{}为空，不进行导表！".format(file_path))

    def add_sg_tree_detail(self, file_path, tree_id, marker_info_path, task_id):
        """
        连锁分群图细节表
        逻辑：根据每行的层级来判断当前层级是否在上一层中、与上一层平级、比上一层高级；记录每一层的上一层信息及母节点
        信息；为每一层标一个序号来区分同名的层级。
        :param file_path: 输入树图的dumper文件。
        :param marker_info_path: 如果不是CP：total.marker.info；如果是CP：total.sexAver.info
        :return:
        参数 last_num:上一层层级数
             last_name:上一层节点名
             mother_num:母节点层级数
             mother_name:母节点名
             dic:键为当前节点名，值为母节点名
             dic_level:键为层级数，值为节点名
        """
        level_list = []
        mother_num = None
        mother_name = None
        last_num = None
        last_name = None
        dic = {}
        dic_level = {}
        node_number = 0
        node_dic = {}
        marker_number = 0
        marker_dic = {}
        lg_detail_list = []
        data_list = []
        insert_data = {
            "name": None,
            "level": 0,
            "value": None,
            "tree_id": tree_id,
            "id": 0,
            "markers": None,
            "p_id": -1
        }
        data_list.append(insert_data)
        if len(data_list) != 0:
            self.col_insert_data("sg_tree_detail", data_list)
        else:
            self.bind_object.logger.info("文件{}为空，不进行导表！".format(file_path))
        with open(file_path, "r") as r:
            lines = r.readlines()
            for line in lines[::2]:
                temp = line.strip().split('\t')
                name = temp[0]
                level_nod = name.strip().split('/')
                if int(level_nod[0]) not in level_list:
                    level_list.append(int(level_nod[0]))
            level_list.sort()
            for line in lines[1::2]:
                temp = line.strip().split(',')
                marker_number += 1
                marker_dic[marker_number] = temp
            for line in lines[:1]:
                temp = line.strip().split('\t')
                name = temp[0]
                node_number += 1
                name_number = name + "-" + str(node_number)
                level_nod = name.strip().split('/')
                level_nod_num = int(level_nod[0])
                mother_num = 0
                mother_name = "0/0"
                last_num = level_nod_num
                last_name = name_number
                dic[name_number] = "0/0"
                node_dic[node_number] = name
            for line in lines[2::2]:
                temp = line.strip().split('\t')
                name = temp[0]
                level_nod = name.strip().split('/')
                level_nod_num = int(level_nod[0])
                node_number += 1
                node_dic[node_number] = name
                name_number = name + "-" + str(node_number)
                if level_nod_num > last_num:
                    dic[name_number] = last_name
                    mother_num = last_num
                    mother_name = last_name
                    last_num = level_nod_num
                    last_name = name_number
                    dic_level[level_nod_num] = name_number
                elif level_nod_num == last_num:
                    dic[name_number] = mother_name
                    last_num = level_nod_num
                    last_name = name_number
                    dic_level[level_nod_num] = name_number
                elif level_nod_num < last_num:
                    if level_nod_num < mother_num:
                        dic[name_number] = dic[dic_level[level_nod_num]]
                        mother_list = dic[dic_level[level_nod_num]].strip().split('/')
                        mother_num = int(mother_list[0])
                        mother_name = dic[dic_level[level_nod_num]]
                        last_num = level_nod_num
                        last_name = name_number
                        dic_level[level_nod_num] = name_number
                    elif level_nod_num == mother_num:
                        dic[name_number] = dic[mother_name]
                        mother_name = dic[mother_name]
                        mother_list = mother_name.strip().split('/')
                        mother_num = int(mother_list[0])
                        last_num = level_nod_num
                        last_name = name_number
                        dic_level[level_nod_num] = name_number
                    elif level_nod_num > mother_num:
                        dic[name_number] = mother_name
                        last_num = level_nod_num
                        last_name = name_number
                        dic_level[level_nod_num] = name_number
            for key, value in dic.items():
                data_list = []
                name_list = key.strip().split('-')
                value_list = value.strip().split('-')
                name = name_list[0]
                number = int(name_list[1])
                markers = marker_dic[number]
                level_nod = name.strip().split('/')
                level = level_list.index(int(level_nod[0])) + 1
                value_name = value_list[0]
                if value_name == "0/0":
                    father_number = 0
                else:
                    father_number = int(value_list[1])
                insert_data = {
                    "name": name,
                    "level": level,
                    "value": name,
                    "tree_id": tree_id,
                    "id": number,
                    "markers": markers,
                    "p_id": father_number
                }
                data_list.append(insert_data)
                if len(data_list) != 0:
                    tree_detail_id = self.col_insert_data("sg_tree_detail", data_list, is_show_log="false")
                    lg_detail_dic = {
                        "name": name,
                        "markers": markers,
                        "number": number,
                        "id": tree_detail_id
                    }
                    lg_detail_list.append(lg_detail_dic)
                else:
                    self.bind_object.logger.info("文件{}为空，不进行导表！".format(file_path))
        self.add_sg_grouping_detail(marker_info_path, lg_detail_list, task_id)

    def add_sg_tree(self, file_path, task_id, location, marker_info_path):
        """
        连锁分群图主表
        :param file_path: 输入树图的dumper文件；
        :param marker_info_path: 如果不是CP：total.marker.info；如果是CP：total.sexAver.info
        :param task_id:
        :param location: grouping_tree
        """
        data_list = []
        level_list = []
        with open(file_path, "r") as r:
            lines = r.readlines()
            for line in lines[::2]:
                temp = line.strip().split('\t')
                name = temp[0]
                level_nod = name.strip().split('/')
                if int(level_nod[0]) not in level_list:
                    level_list.append(int(level_nod[0]))
            level_list.sort()
            max_level = len(level_list)
            insert_data = {
                "level_num": max_level,
                "task_id": task_id,
                "location": location,
                "type": 1
            }
            data_list.append(insert_data)
        if len(data_list) != 0:
            tree_id = self.col_insert_data("sg_tree", data_list)
            self.update_db_record("sg_tree", {"_id": tree_id}, {"main_id": tree_id})
            self.add_sg_tree_detail(file_path, tree_id, marker_info_path, task_id)
        else:
            self.bind_object.logger.info("文件{}为空，不进行导表！".format(file_path))

    def add_sg_grouping_detail(self, file_path, lg_detail_list, task_id):
        """
        标记连锁分群结果表
        :param file_path: 如果不是CP：total.marker.info；如果是CP：total.sexAver.info
        :return:
        """
        data_list = []
        with open(file_path, "r") as r:
            lines = r.readlines()
            for line in lines[1:]:
                temp = line.strip().split('\t')
                markerid = temp[0]
                for i in lg_detail_list:
                    if markerid in i["markers"]:
                        marker_list = markerid.strip().split('_')
                        if len(temp) == 10:
                            chr_id = marker_list[0]
                            physical_pos = marker_list[-1]
                            type = temp[3]
                            nind = temp[4]
                            nmiss = temp[5]
                            geno = temp[6]
                            sgeno = temp[7]
                            signif = temp[8]
                            segretion = temp[9]
                            tree_id = i["id"]
                            insert_data = {
                                "marker_id": markerid,
                                "chr_id": chr_id,
                                "physical_pos": int(physical_pos.split('-')[0]),
                                "type": type,
                                "nind": int(nind),
                                "nmiss": int(nmiss),
                                "geno": geno,
                                "sgeno": sgeno,
                                "segretion": segretion,
                                "signif": float(signif),
                                "tree_detail_id": tree_id,
                                "task_id": task_id
                            }
                            data_list.append(insert_data)
        if len(data_list) != 0:
            self.col_insert_data("sg_grouping_detail", data_list)
        else:
            self.bind_object.logger.info("文件{}为空，不进行导表！".format(file_path))

    def add_sg_lg_stat(self, file_path, lg_id):
        """
        标记信息统计表
        :param file_path: Tool:get_grouping_result生成的total.marker.stat.xls文件。
        :param lg_id:主表ID
        :return:
        """
        data_list = []
        with open(file_path, "r") as r:
            lines = r.readlines()
            for line in lines[1:]:
                temp = line.strip().split('\t')
                lg = temp[0]
                marker_number = temp[1]
                snp_number = temp[2]
                indel_number = temp[3]
                insert_data = {
                    "lg": lg,
                    "marker_num": int(marker_number),
                    "snp_num": int(snp_number),
                    "indel_num": int(indel_number),
                    "lg_id": lg_id
                }
                data_list.append(insert_data)
            if data_list[0]:
                if re.match("^LG", data_list[0]["lg"]):
                    data_list.sort(key=lambda i: int(re.findall('\d+', i["lg"])[0]))
                else:
                    data_list.sort(key=lambda i: int(i["lg"]))
            if len(data_list) != 0:
                self.col_insert_data("sg_lg_stat", data_list)
            else:
                self.bind_object.logger.info("文件{}为空，不进行导表！".format(file_path))

    def add_sg_lg_detail(self, lg_path, map_path, marker_path, lg_id):
        """
        1.输入Total.lg，pop.filter.marker，total.map三个文件，将Total.lg和total.map整成字典，根据pop.filter.marker的每一行
        来取对应的数据。(如果是CP，就用中性的map文件；如果是bin，就用bin的结果文件来做，total.bin.marker。）
        2.动态加载samples：根据pop.filter.marker的第一行来取到samples，再每次添加一条记录时添加对应的sample和
        值到insert_data中。
        :param lg_path: Total.lg
        :param map_path: total.map
        :param marker_path: pop.filter.marker
        :param lg_id:
        :return:
        """
        data_list = []
        dict_lg = {}
        dict_map = {}
        value_list = []
        key_list = []
        lg_line_number = 0
        marker_path_row = 0
        sample_dict = {}
        sample_list = []
        lg_id = self.check_objectid(lg_id)
        with open(lg_path, "r") as r:
            lines = r.readlines()
            for line in lines:
                lg_line_number += 1
                if (lg_line_number % 2) == 0:
                    value = line.strip().split('\t')
                    value_list.append(value)
                else:
                    temp = line.strip().split('\t')
                    lg_id_list = temp[0].strip().split('>')
                    key = lg_id_list[1]
                    key_list.append(key)
            for i in range(len(value_list)):
                dict_lg[key_list[i]] = value_list[i]
        with open(map_path, "r")as m:
            lines = m.readlines()
            for line in lines[1:]:
                temp = line.strip().split('\t')
                key = temp[0]
                value = temp[1]
                dict_map[key] = value
        with open(marker_path, "r") as p:
            lines = p.readlines()
            for line in lines:
                marker_path_row += 1
                temp = line.strip().split('\t')
                marker_path_sample_num = len(temp) - 2
                if marker_path_row == 1:
                    for i in range(marker_path_sample_num):
                        sample_dict[(i + 2)] = temp[(i + 2)]
                else:
                    marker_id = temp[0]
                    count = 0
                    if marker_id in dict_map.keys():
                        lg_pos = dict_map[marker_id]
                    else:
                        lg_pos = "--"
                    for key, value in dict_lg.items():
                        if marker_id in value:
                            lg = key
                            count += 1
                            chr_list = marker_id.strip().split('_')
                            chr = chr_list[0]
                            chr_pos = chr_list[1]
                            insert_data = {
                                "marker_id": marker_id,
                                "lg": lg,
                                "lg_pos": str(lg_pos),
                                "chr": chr,
                                "chr_pos": str(chr_pos),
                                "type": temp[1],
                                "lg_id": lg_id
                            }
                            for x, y in sample_dict.items():
                                insert_data[y] = temp[x]
                            data_list.append(insert_data)
                    if count == 0:
                        lg = "--"
                        chr_list = marker_id.strip().split('_')
                        chr = chr_list[0]
                        chr_pos = chr_list[1]
                        insert_data = {
                            "marker_id": marker_id,
                            "lg": lg,
                            "lg_pos": str(lg_pos),
                            "chr": chr,
                            "chr_pos": str(chr_pos),
                            "type": temp[1],
                            "lg_id": lg_id
                        }
                        for x, y in sample_dict.items():
                            insert_data[y] = temp[x]
                        data_list.append(insert_data)
        for x, y in sample_dict.items():
            sample_list.append(y)
        if len(data_list) != 0:
            self.col_insert_data("sg_lg_detail", data_list)
            self.update_db_record("sg_lg", {"_id": lg_id}, {"specimen_ids": sample_list})
        else:
            self.bind_object.logger.info("文件{}为空，不进行导表！".format(marker_path))

    def update_sg_lg(self, lg_id, loc_path, map_path, types="cp", csv_path=None):
        lg_id = self.check_objectid(lg_id)
        if types == "cp":
            self.update_db_record("sg_lg", {"_id": lg_id}, {"sexAver_loc_path": loc_path, "sexAver_map_path": map_path})
        else:
            self.update_db_record("sg_lg", {"_id": lg_id}, {"total_loc_path": loc_path, "total_map_path": map_path,
                                                            "total_csv_path": csv_path})

    def add_sg_reorganization_heatmap(self, pic_path, target_path, evalutaion_id, types="nocp"):
        """
        pic_path:真实的热图数据路径
        target_path：远程磁盘路径
        双坐标热图
        :return:
        """
        evalutaion_id = self.check_objectid(evalutaion_id)
        data_list = []
        if types not in ['cp', "nocp"]:
            raise Exception("type:{}不合法！".format(types))
        insert_data = {
            "evalutaion_id": evalutaion_id,
            "created_ts": datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
            "heatmap_path": target_path,
            "genetype_path": os.path.join(os.path.dirname(target_path), "fig1")
        }
        if types == "cp":
            chrs_list = []
            for m in os.listdir(pic_path):
                n = re.match(r"(.*)\.heatMap\.sexaver\.png$", m)
                if n:
                    chrs_list.append(int(n.group(1)))
            chrs_list.sort()
            insert_data.update({"chrs_list": chrs_list})
            insert_data.update({"data_type": "female"})
            data_list.append(insert_data)
            insert_data1 = copy.deepcopy(insert_data)
            insert_data1.update({"data_type": "male"})
            data_list.append(insert_data1)
            insert_data2 = copy.deepcopy(insert_data)
            insert_data2.update({"data_type": "sexaver"})
            data_list.append(insert_data2)
        else:
            chrs_list = []
            for m in os.listdir(pic_path):
                n = re.match(r"(.*)\.heatMap\.png$", m)
                if n:
                    chrs_list.append(int(n.group(1)))
            chrs_list.sort()
            insert_data.update({"chrs_list": chrs_list})
            insert_data.update({"data_type": "total"})
            data_list.append(insert_data)
        if data_list:
            self.col_insert_data("sg_reorganization_heatmap", data_list)
        else:
            self.bind_object.logger.info("文件{}为空，不进行导表！")

if __name__ == "__main__":
    data = LinkageGroupping(None)
    # data.add_genetype_heatmap(
    #     "/mnt/ilustre/users/sanger-dev/workspace/20180727/Gmap_tsg_31236/output/04.grouping/evalutaion/total.csv",
    #     "total", "tsg_31236", "5b5adb31a4e1af018b2ccbad")
    #
    # data.add_yichuantupu("tsg_31236", "5b5adb31a4e1af018b2ccbad",
    #                      "/mnt/lustre/users/sanger/workspace/20181204/Gmap_sanger_143999/output/"
    #                      "04.grouping/evalutaion/total.map", "total", updata_chr_list=False)
    # data.add_yichuantupu("tsg_31236", "5b5adb31a4e1af018b2ccbad",
    #                      "/mnt/ilustre/users/sanger-dev/workspace/20180801/Gmap_tsg_31303/output/04.grouping/ev"
    #                      "alutaion/total.male.map", "male",
    #                      updata_chr_list=False)
    # data.add_yichuantupu("tsg_31236", "5b5adb31a4e1af018b2ccbad",
    #                      "/mnt/ilustre/users/sanger-dev/workspace/20180801/Gmap_tsg_31303/output/04.grouping/e"
    #                      "valutaion/total.sexAver.map", "sexaver",
    #                      updata_chr_list=False)
    # data.collinearity("test111",
    #                   "/mnt/lustre/users/sanger/workspace/20181204/Gmap_sanger_143999/output/04.grouping/evalutaion/total.map",
    #                   "5b5adb31a4e1af018b2ccbad",
    #                   "/mnt/lustre/users/sanger/workspace/20181204/Gmap_sanger_143999/output/04.grouping/evalutaion"
    #                   "/total.phy.spearman.xls", "total")
    # data.add_sg_evalutaion_detail("/mnt/lustre/users/sanger/workspace/20181204/Gmap_sanger_143999/output/04.grouping"
    #                               "/evalutaion/total.marker.info", "5b5adb31a4e1af018b2ccbad")
    # data.add_sg_lg_detail("/mnt/lustre/users/sanger/workspace/20181204/Gmap_sanger_143999/output/04.grouping/groupping/Total.lg",
    #                  "/mnt/lustre/users/sanger/workspace/20181204/Gmap_sanger_143999/output/04.grouping/evalutaion/total.map",
    #                       "/mnt/lustre/users/sanger/workspace/20181204/Gmap_sanger_143999/output/02.marker_filter/pop.filtered.marker", "5b5adb31a4e1af018b2ccbad")
    marker_info_path = "/mnt/lustre/users/sanger/workspace/20181204/Gmap_sanger_143999/output/04.grouping/evalutaion/total.marker.info"
    # data.add_sg_lg_stat("/mnt/lustre/users/sanger/workspace/20181204/Gmap_sanger_143999/output/04.grouping/marker_stat/total.marker.stat.xls", "5b5adb31a4e1af018b2ccbad")
    data.add_sg_tree("/mnt/lustre/users/sanger/workspace/20181204/Gmap_sanger_143999/output/temp/tree_data/Total.gTree.hash.dumper",
                        "111", "grouping_tree", marker_info_path)
    # data.add_sg_lg_stat(
    #     "/mnt/ilustre/tsanger-data/rerewrweset/files/m_188/188_5b7a78fa03634/tsanger_31644/interaction_results/LinkageGrouping/LinkageGrouping_by_hand_20180821_183315869665/marker_stat/total.marker.stat.xls",
    #     "test_1234567")
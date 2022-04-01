# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.api.database.base import Base, report_check
import re
import os
import json
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON


class AlphaCompare(Base):
    def __init__(self, bind_object):
        super(AlphaCompare, self).__init__(bind_object)
        self._project_type = "metaasv"

    def cat_files(self, statfile, cifiles=None):
        '''
        将组间差异比较的统计结果文件与posthoc结果文件cat成一个文件
        cifiles为posthoc文件路径的列表
        '''
        data = {}
        sort_list = []
        with open(statfile, 'rb') as s:
            head = s.readline().strip('\n').split()
            for line in s:
                line = line.strip('\n').split('\t')
                info = {}
                for i in head:
                    info[i] = line[head.index(i) + 1]
                data[line[0]] = info
                sort_list.append(line[0])
        if cifiles:
            for f in cifiles:
                with open(f, 'rb') as r:
                    head = r.readline().strip('\n').split('\t')
                    for line in r:
                        line = line.strip('\n').split('\t')
                        if line[0] in sort_list:
                            info = {}
                            for i in head[1:]:
                                info[i] = line[head.index(i)]
                            data[line[0]].update(info)
                        else:
                            pass
        if cifiles and len(cifiles) > 1:
            for name in data:
                keys = data[name].keys()
                mean = []
                low_ci = []
                ci = []
                for i in keys:
                    if re.search(r'mean$', i):
                        mean.append(float(data[name][i]))
                    elif re.search(r'lowerCI$', i):
                        group = i.split('_lowerCI')[0]
                        up = group + '_upperCI'
                        low_ci.append(float(data[name][i]))
                        ci.append(float(data[name][up]) - float(data[name][i]))
                max_mean = max(mean)
                max_ci = max(ci)
                min_ci = min(low_ci)
                if max_ci <= max_mean:
                    n = 1
                else:
                    n = round(max_ci / max_mean)
                l = round(abs(min_ci / n) + max_mean + 3)
                data[name].update({'n': n, 'l': l})
        return data, sort_list

    @report_check
    def add_species_difference_check_detail(self, statfile, cifiles, level=None, check_type=None, params=None, category_name=None, table_id=None, group_id=None, from_otu_table=None, major=False, posthoc=None, correlation_key=None, coll_name='multiple_group_detail', main_coll='multiple_group', type=None):
        if major:
            table_id = self.create_species_difference_check(level, check_type, params, category_name, group_id, from_otu_table)
        else:
            if table_id is None:
                self.bind_object.set_error("major为False时需提供table_id!")
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
                else:
                    self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!")
        stat_info, sort_list = self.cat_files(statfile, cifiles)
        self.bind_object.logger.info("stat_info: %s" % stat_info)
        self.bind_object.logger.info("sort_list: %s" % sort_list)
        with open(statfile, 'r') as f:
            table_name_list = f.readline().strip().split("\t")
        data_list = []
        # pvalue_list = []
        # qvalue_list = []
        posthoc_list = []
        for name in sort_list:
            data = [
                (correlation_key, table_id),
                ("estimators", name),
            ]
            for i in stat_info[name].keys():
                if i in ['pvalue', 'qvalue', 'corrected_pvalue', 'statistic', 'odds_ratio'] and stat_info[name][i] == 'NA':
                    if i == 'qvalue' or  i == 'corrected_pvalue':
                        data.append(('qvalue', stat_info[name][i]))
                    else:
                        data.append((i, stat_info[name][i]))
                elif re.search(r'_pvalue$', i) and posthoc in ['tukeykramer', 'gameshowell']:
                    data.append((i, stat_info[name][i]))
                    if i not in posthoc_list:
                        posthoc_list.append(i)
                elif i == 'qvalue':
                    data.append(('qvalue', float(stat_info[name][i])))
                else:
                    data.append((i, float(stat_info[name][i])))
                # if re.search(r'_pvalue$', i):
                #     if i not in pvalue_list:
                #         pvalue_list.append(i)
                # if re.search(r'_qvalue$', i):
                #     if i not in pvalue_list:
                #         qvalue_list.append(i)
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db[coll_name]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" %(coll_name,e))
            self.bind_object.set_error("导入%s信息出错"%(coll_name))
        else:
            self.bind_object.logger.info("导入%s信息成功!"%coll_name)
        try:
            main_collection = self.db[main_coll]
            settled_params = {"software" : "R-3.3.1 (stat)"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            if type in ["multiple"]:
                table_data = {"table_data": ["estimators"] + table_name_list + posthoc_list}
                table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id": table_id},{"$set": {"main_id": table_id,
                                                            "settled_params":settled_params_json,
                                                          "table_data":table_data_json}})
            elif type in ["two_group"]:
                table_data = {"table_data": ["estimators"] + table_name_list+ ["pvalue", "qvalue"]}
                table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id": table_id},{"$set": {"main_id": table_id,
                                                        "settled_params":settled_params_json,
                                                          "table_data":table_data_json}})
        except Exception, e:
            self.bind_object.logger.error("更新%s信息出错:%s" % (main_coll, e))
            self.bind_object.set_error("更新信息出错")
        else:
            self.bind_object.logger.info("更新%s信息成功!" % main_coll)
        return table_id

    @report_check
    def add_species_difference_check_barplot(self, file_path, table_id, correlation_key=None, coll_name='multiple_group_bar', main_coll='multiple_group'):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        index_list = []
        group_list = []
        with open(file_path, 'rb') as r:
            l = r.readline().strip('\n')
            l_all = l.split("\t")
            for i in l_all:
                if re.search(r"mean$", i):
                    if i not in group_list:
                        new_name = i.split("-mean")[0]
                        group_list.append(new_name)
            while True:
                line = r.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                if line_data[0] not in index_list:
                    index_list.append(line_data[0])
                i = 1
                for name in group_list:
                    data = [("specimen", name),(correlation_key, table_id),("type", line_data[0])]
                    mean = float(line_data[i])
                    sd_value = float(line_data[i +1 ])
                    upper = mean + sd_value
                    lower = mean - sd_value
                    data.append(("mean", float(line_data[i])))
                    data.append(("lower", lower))
                    data.append(("upper", upper))
                    i += 2
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db[coll_name]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入单物种柱状图信息出错:%s" % e)
            self.bind_object.set_error("导入单物种柱状图信息出错")
        else:
            self.bind_object.logger.info("导入单物种柱状图信息成功!")
        try:
            main_collection = self.db[main_coll]
            column_data = {
                "column_data": {"name":"specimen",
                            "data": "mean",
                            "condition": {"type": index_list}}}
            single_column_data_json = json.dumps(column_data, sort_keys=True, separators=(',', ':'))
            ishape_data = {
                "ishape_data": {"name":"specimen",
                            "data": ["mean", "lower", "upper"],
                            "condition": {"type": index_list}}}
            ishape_data_json = json.dumps(ishape_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": table_id},{"$set": {"main_id": table_id,
                                                               "column_data":single_column_data_json,
                                                               "ishape_data":ishape_data_json}})
        except Exception, e:
            self.bind_object.logger.error("更新%s信息出错:%s" % (main_coll, e))
            self.bind_object.set_error("更新信息出错")
        else:
            self.bind_object.logger.info("更新%s信息成功!" % main_coll)

    @report_check
    def add_species_difference_check_boxplot(self, file_path, table_id, correlation_key=None, coll_name='multiple_group_box', main_coll='multiple_group'):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        index_list = []
        with open(file_path, 'rb') as r:
            l = r.readline().strip('\n')
            group_list = re.findall(r'min\((.*?)\)', l)
            while True:
                line = r.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                if line_data[0] not in index_list:
                    index_list.append(line_data[0])
                i = 1
                for name in group_list:
                    data = [(correlation_key, table_id), ("type", line_data[0])]
                    data.append(("specimen", name))
                    data.append(("min", float(line_data[i])))
                    data.append(("q1", float(line_data[i + 1])))
                    data.append(("median", float(line_data[i + 2])))
                    data.append(("q3", float(line_data[i + 3])))
                    data.append(("max", float(line_data[i + 4])))
                    i += 5
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db[coll_name]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (coll_name, e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % coll_name)
        try:
            main_collection = self.db[main_coll]
            box_data = {
                "box_data": {"name":"specimen",
                            "data": ["min","q1","median","q3","max"],
                            "condition": {"type": index_list}}}
            box_data_json = json.dumps(box_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": table_id},{"$set": {"main_id": table_id,"box_data":box_data_json}})
        except Exception, e:
            self.bind_object.logger.error("更新%s信息出错:%s" % (main_coll, e))
            self.bind_object.set_error("更新信息出错")
        else:
            self.bind_object.logger.info("更新%s信息成功!" % main_coll)

    def add_est_detail(self,statfile, table_id, correlation_key=None, coll_name='multiple_group_detail',main_coll='multiple_group', type=None):
        """
        指数组间差异两组比较的结果进行导表
        :param statfile: 结果文件夹
        :param table_id: 主表id
        :param correlation_key:主表与详情表的关联字段
        :param coll_name: 导入详情表名称
        :param main_coll: 更新主表名称
        :param type: type类型
        :return:
        由于两组比较的数据，在不同的文件存在重复，需要注意去掉这部分的数据
        """
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!")
        dirs = os.listdir(statfile)
        box_file_list = []
        stat_file_list = []
        for file in dirs:
            file_path = os.path.join(statfile, file)
            if re.search(r'est_result', file):
                if file_path not in stat_file_list:
                    stat_file_list.append(file_path)
            elif re.search(r'est_barplot_result', file):
                if file_path not in box_file_list:
                    box_file_list.append(file_path)
        all_group_list = []
        data_list = []
        index_list = []
        self.bind_object.logger.info("box_file_list:{}".format(box_file_list))
        for file in box_file_list:
            with open(file, "r") as m:
                all_lines = len(m.readlines())
            j = 1
            file_group_list = [] # 用于记录每个文件的分组名称
            with open(file, "r") as f:
                self.bind_object.logger.info("file_path:{}".format(file))
                l = f.readline().strip()
                group_list = re.findall(r'min\((.*?)\)', l)
                while True:
                    line = f.readline().strip('\n')
                    j += 1
                    if not line:
                        break
                    line_data = line.split("\t")
                    if line_data[0] not in index_list:
                        index_list.append(line_data[0])
                    i = 1
                    for name in group_list:
                        data = [(correlation_key, table_id), ("type", line_data[0])]
                        data.append(("specimen", name))
                        data.append(("min", float(line_data[i])))
                        data.append(("q1", float(line_data[i + 1])))
                        data.append(("median", float(line_data[i + 2])))
                        data.append(("q3", float(line_data[i + 3])))
                        data.append(("max", float(line_data[i + 4])))
                        i += 5
                        if name not in all_group_list:
                            data_son = SON(data)
                            data_list.append(data_son)
                            if name not in file_group_list:
                                file_group_list.append(name)
                    if j == all_lines:
                        all_group_list += file_group_list
                        all_group_list = list(set(all_group_list))

        bar_list = []
        all_bar_group_list = []
        all_stat_list = []
        stat_list = []
        table_name_list = []
        for file in stat_file_list:
            new_group_list = []
            with open(file, "r") as m:
                all_lines2 = len(m.readlines())
            self.bind_object.logger.info("file: {}".format(file))
            all_group_list2 = []###用于记录文件本身的分组名称
            x = 1
            with open(file, "r") as f:
                l = f.readline().strip()
                l_all = l.split("\t")
                for i_line in l_all:
                    if re.search(r"mean$", i_line):
                        i_name = i_line.split("-mean")[0]
                        if i_name not in new_group_list:
                            new_group_list.append(i_name)
                pvalue = ""
                qvalue = ""
                while True:
                    line = f.readline().strip()
                    x += 1
                    if not line:
                        break
                    line_data = line.split("\t")
                    if line_data[0] not in index_list:
                        index_list.append(line_data[0])
                    for name in new_group_list:
                        # self.bind_object.logger.info("new_group_list: {}".format(new_group_list))
                        # self.bind_object.logger.info("all_bar_list: {}".format(all_group_list2))
                        data = [("specimen", name),(correlation_key, table_id),("type", line_data[0])]
                        group_name = name+"-mean"
                        # self.bind_object.logger.info("group_name: {}".format(group_name))
                        index = l_all.index(group_name)
                        mean = float(line_data[index + 1])
                        sd_value = float(line_data[index + 2])
                        upper = mean + sd_value
                        lower = mean - sd_value
                        data.append(("mean", float(mean)))
                        data.append(("lower", lower))
                        data.append(("upper", upper))
                        data_son2 = SON(data)
                        if name not in all_bar_group_list:
                            bar_list.append(data_son2)
                            if name not in all_group_list2:
                                all_group_list2.append(name)
                    if x == all_lines2:
                        all_bar_group_list += all_group_list2
                        all_bar_group_list = list(set(all_bar_group_list))

                    group_string = "-".join(new_group_list)
                    j = 1
                    stat_data = [("estimators", line_data[0]),(correlation_key, table_id),("type", type)]
                    for name in new_group_list:
                        if name not in all_stat_list:
                            group_name = name+"-mean"
                            index = l_all.index(group_name)
                            stat_data.append(("{}-mean".format(name), float(line_data[index + 1])))
                            stat_data.append(("{}-sd".format(name), float(line_data[index + 2])))
                            table_name_list.append(name + "-mean")
                            table_name_list.append(name + "-sd")
                        j += 2
                    pvalue = group_string + "_pvalue"
                    qvalue = group_string + "_qvalue"
                    stat_data.append((pvalue, line_data[-2]))
                    stat_data.append((qvalue, line_data[-1]))
                    data_son = SON(stat_data)
                    if name not in all_stat_list:
                        if name not in all_group_list2:
                            all_group_list2.append(name)
                        stat_list.append(data_son)
                    if x == all_lines2:
                        all_stat_list += all_group_list2
                        all_stat_list = list(set(all_stat_list))
                table_name_list.append(pvalue)
                table_name_list.append(qvalue)

        try:
            collection = self.db["alpha_diversity_diff_box"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % ("alpha_diversity_diff_box", e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % "alpha_diversity_diff_box")

        try:
            collection = self.db["alpha_diversity_diff_table"]
            collection.insert_many(stat_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % ("alpha_diversity_diff_bar", e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % "alpha_diversity_diff_bar")

        try:
            collection = self.db["alpha_diversity_diff_bar"]
            collection.insert_many(bar_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % ("alpha_diversity_diff_bar", e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % "alpha_diversity_diff_bar")

        try:
            main_collection = self.db["alpha_diversity_diff"]
            box_data = {
                "box_data": {"name":"specimen",
                            "data": ["min","q1","median","q3","max"],
                            "condition": {"type": index_list}}}
            box_data_json = json.dumps(box_data, sort_keys=True, separators=(',', ':'))
            column_data = {
                "column_data": {"name":"specimen",
                            "data": "mean",
                            "condition": {"type": index_list}}}
            single_column_data_json = json.dumps(column_data, sort_keys=True, separators=(',', ':'))
            ishape_data = {
                "ishape_data": {"name":"specimen",
                            "data": ["mean", "lower", "upper"],
                            "condition": {"type": index_list}}}
            ishape_data_json = json.dumps(ishape_data, sort_keys=True, separators=(',', ':'))
            settled_params = {"software" : "R-3.3.1 (stat)"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))

            table_data = {"table_data": ["estimators"] + table_name_list}
            table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))

            main_collection.update({"_id": table_id},{"$set": {"main_id": table_id,
                                                               "box_data":box_data_json,
                                                               "column_data":single_column_data_json,
                                                               "ishape_data":ishape_data_json,
                                                               "settled_params":settled_params_json,
                                                               "table_data":table_data_json}})
        except Exception, e:
            self.bind_object.logger.error("更新%s信息出错:%s" % (main_coll, e))
            self.bind_object.set_error("更新信息出错")
        else:
            self.bind_object.logger.info("更新%s信息成功!" % main_coll)


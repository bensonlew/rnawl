# -*- coding: utf-8 -*-

from biocluster.api.database.base import Base, report_check
import re
import json
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON



class BugbaseGroup(Base):
    def __init__(self, bind_object):
        super(BugbaseGroup, self).__init__(bind_object)
        # self._db_name = Config().MONGODB
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
            if head[0] == "name":
                head = head[1:]
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
    def add_species_difference_check_detail(self, statfile, cifiles, check_type=None, params=None, category_name=None, table_id=None, group_id=None, from_otu_table=None, major=False, posthoc=None, correlation_key=None, coll_name='multiple_group_detail', main_coll='multiple_group',bar_name="bugbase_group_bar",task_type=None):
        if major:
            table_id = self.create_species_difference_check("9",check_type, params, category_name, group_id, from_otu_table)
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
        all_function=[]
        with open(statfile, 'r') as f:
            table_name_list = f.readline().strip().split("\t")
            table_name_list.remove("statistic")
            if "name" in table_name_list:
                table_name_list.remove("name")
            data = f.readlines()
            for xx in data:
                all_function.append(xx.split('\t')[0])
        species_ishape = []
        for f in cifiles:
            with open(f, 'r') as inf:
                group_diff_list = inf.readline().strip().split()[0:3]
            if group_diff_list not in species_ishape:
                species_ishape.append(group_diff_list)
        data_list = []
        mean_list = []
        sd_list = []
        pvalue_list = []
        effective_list = []
        propotion_list = []
        for name in sort_list:
            all_mean = 0
            data = [
                (correlation_key, table_id),
                ("function_name", name)
            ]
            for i in stat_info[name].keys():
                if i in ['pvalue', 'qvalue', 'corrected_pvalue', 'statistic', 'odds_ratio'] and stat_info[name][i] == 'NA':
                    if i == 'qvalue' or  i == 'corrected_pvalue':
                        data.append(('qvalue', stat_info[name][i]))
                    else:
                        data.append((i, stat_info[name][i]))
                elif re.search(r'_pvalue$', i) and posthoc in ['tukeykramer', 'gameshowell']:
                    data.append((i, stat_info[name][i]))
                elif i == 'qvalue':
                    data.append(('qvalue', float(stat_info[name][i])))
                elif i == 'corrected_pvalue':
                    data.append(('qvalue', float(stat_info[name][i])))
                else:
                    if task_type in ["two_group"]:
                        data.append((i, float(stat_info[name][i])))
                    else:
                        data.append((i, float(stat_info[name][i])))
                if re.search(r'_pvalue$', i):
                    if i not in pvalue_list:
                        pvalue_list.append(i)
                elif re.search(r'-mean$', i):
                    if i not in mean_list:
                        mean_list.append(i)
                    all_mean += float(stat_info[name][i])
                elif re.search(r'-sd$', i):
                    if i not in sd_list:
                        sd_list.append(i)
                elif re.search(r'effectsize$', i):
                    if i not in effective_list:
                        effective_list.append(i)
                elif re.search(r'propotion$', i):
                    if i not in propotion_list:
                        propotion_list.append(i.strip("-propotion"))
            data.append(("all_value", all_mean))
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
            if task_type in ["multiple_group"]:
                table_data = {"table_data": ["function_name"] + table_name_list}
                table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
                column_data = {
                    "column_data": {"name":"function_name",
                                "data": mean_list,"category": ""}}
                column_data_json = json.dumps(column_data, sort_keys=True, separators=(',', ':'))
                ishape_data = {
                    "ishape_data": {"name":"function_name",
                                "data": sd_list,"category": ""}}
                ishape_data_json = json.dumps(ishape_data, sort_keys=True, separators=(',', ':'))
                species_column_data = {
                    "column_data": {"name":"function_name",
                                "data": effective_list,"category": ""}}
                species_column_data_json = json.dumps(species_column_data, sort_keys=True, separators=(',', ':'))
                species_text_data = {
                    "text_data": {"name":"function_name",
                                "data": pvalue_list,"category": ""}}
                species_text_data_json = json.dumps(species_text_data, sort_keys=True, separators=(',', ':'))

                multiple_text_data = {
                    "text_data": {"name":"function_name",
                                "data": ["pvalue", "qvalue"],"category": ""}}
                multiple_text_data_json = json.dumps(multiple_text_data, sort_keys=True, separators=(',', ':'))
                species_ishape_data = {
                    "ishape_data": {"name":"function_name",
                                "data": species_ishape,"category": ""}}
                species_ishape_data_json = json.dumps(species_ishape_data, sort_keys=True, separators=(',', ':'))

                main_collection.update({"_id": table_id},{"$set": {"main_id": table_id,
                                                            "species_column_data":species_column_data_json,
                                                          "table_data":table_data_json,
                                                           "column_data":column_data_json,
                                                           "settled_params":settled_params_json,
                                                           "ishape_data":ishape_data_json,
                                                           "species_text_data":species_text_data_json,
                                                           "species_ishape_data":species_ishape_data_json,
                                                           "multiple_text_data":multiple_text_data_json,
                                                           "all_function":",".join(sort_list)}})
            elif task_type in ["two_group"]:
                table_data = {"table_data": ["function_name"] + table_name_list[0:-1] + ["qvalue"]}
                table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
                column_data = {
                    "column_data": {"name":"function_name",
                                "data": mean_list,"category": ""}}
                column_data_json = json.dumps(column_data, sort_keys=True, separators=(',', ':'))
                left_ishape_data = {
                    "ishape_data": {"name":"function_name",
                                "data": sd_list,"category": ""}}
                ishape_data_json = json.dumps(left_ishape_data, sort_keys=True, separators=(',', ':'))
                species_column_data = {
                    "column_data": {"name":"function_name",
                                "data": effective_list,"category": ""}}
                species_column_data_json = json.dumps(species_column_data, sort_keys=True, separators=(',', ':'))

                multiple_text_data = {
                    "text_data": {"name":"function_name",
                                "data": ["pvalue", "qvalue"],"category": ""}}
                multiple_text_data_json = json.dumps(multiple_text_data, sort_keys=True, separators=(',', ':'))
                right_ishape_data = {
                    "ishape_data": {"name":"function_name",
                                "data": species_ishape,"category": ""}}
                species_ishape_data_json = json.dumps(right_ishape_data, sort_keys=True, separators=(',', ':'))

                main_collection.update({"_id": table_id},{"$set": {"main_id": table_id,
                                                            "species_column_data":species_column_data_json,
                                                          "table_data":table_data_json,
                                                           "column_data":column_data_json,
                                                           "settled_params":settled_params_json,
                                                           "left_ishape_data":ishape_data_json,
                                                           "right_ishape_data":species_ishape_data_json,
                                                           "multiple_text_data":multiple_text_data_json,
                                                           "all_function":",".join(sort_list)}})
        except Exception, e:
            self.bind_object.logger.error("更新%s信息出错:%s" % (main_coll, e))
            self.bind_object.set_error("更新信息出错")
        else:
            self.bind_object.logger.info("更新%s信息成功!" % main_coll)

        """
        导入bar的数据，主要为了pvalue
        """

        data_list2=[]
        if task_type in ["two_group"]:
            sample_list = []
            for i in stat_info[name].keys():
                if re.search(r'-mean$', i):
                    if i.split("-mean")[0] not in sample_list:
                        sample_list.append(i.split("-mean")[0])
            self.bind_object.logger.info("sample_list: %s" % sample_list)
            for name in sort_list:
                data = [
                    (correlation_key, table_id),
                    ("function_name", name)
                ]
                data.append(("pvalue", float(stat_info[name]["pvalue"])))
                data.append(("effectsize", float(stat_info[name]["effectsize"])))
                data.append(("upperCI", float(stat_info[name]["upperCI"])))
                data.append(("lowerCI", float(stat_info[name]["lowerCI"])))
                data.append(("name", sample_list[0] + "-" + sample_list[1]))
                data_son = SON(data)
                self.bind_object.logger.info("data_son: %s" % data_son)
                self.bind_object.logger.info("data_son: %s" % data_son)
                data_list2.append(data_son)
        else:
            sample_list = []
            for i in stat_info[name].keys():
                if re.search(r'_effectsize$', i):
                    sample_list.append(i.split("_effectsize")[0])
            for name in sort_list:
                data = [
                    (correlation_key, table_id),
                    ("function_name", name)
                ]
                for xxx in sample_list:
                    if xxx + "_effectsize" in stat_info[name].keys():
                        data.append(("effectsize", float(stat_info[name][xxx + "_effectsize"])))
                    if xxx + "_lowerCI" in stat_info[name].keys():
                        data.append(("lowerCI", float(stat_info[name][xxx + "_lowerCI"])))
                    if xxx + "_upperCI" in stat_info[name].keys():
                        data.append(("upperCI", float(stat_info[name][xxx + "_upperCI"])))
                    if xxx + "_pvalue" in stat_info[name].keys():
                        if " " in stat_info[name][xxx + "_pvalue"]:
                            pvalue = float(stat_info[name][xxx + "_pvalue"].split(" ")[-1])
                        else:
                            pvalue = float(stat_info[name][xxx + "_pvalue"])
                        data.append(("pvalue", pvalue))
                    data.append(("name", xxx))
                    data_son = SON(data)
                    data_list2.append(data_son)
        try:
            collection = self.db[bar_name]
            collection.insert_many(data_list2)
        except Exception, e:
            self.bind_object.logger.error("更新%s信息出错:%s" % (bar_name, e))
            self.bind_object.set_error("更新信息出错")
        else:
            self.bind_object.logger.info("更新%s信息成功!" % bar_name)
        return table_id

    @report_check
    def add_species_difference_check_boxplot(self, file_path, table_id, correlation_key=None,
                                             coll_name='multiple_group_box', coll_name2 = 'multiple_group_',main_coll='multiple_group',task_type=None):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        data_list2 = []
        with open(file_path, 'rb') as r:
            l = r.readline().strip('\n')
            group_list = re.findall(r'min\((.*?)\)', l)
            while True:
                line = r.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                i = 1
                for name in group_list:
                    data = [(correlation_key, table_id), ("function_name", line_data[0])]
                    data.append(("specimen", name))
                    data.append(("min", float(line_data[i])*100))
                    data.append(("q1", float(line_data[i + 1])*100))
                    data.append(("median", float(line_data[i + 2])*100))
                    data.append(("q3", float(line_data[i + 3])*100))
                    data.append(("max", float(line_data[i + 4])*100))
                    #data2 = [(correlation_key, table_id), ("function_name", line_data[0])]
                    #data.append(("name", float(line_data[i + 2])))
                    i += 5
                    data_son = SON(data)
                    data_list.append(data_son)
                    #data_son2 = SON(data2)
                    #data_list2.append(data_son2)
        try:
            collection = self.db[coll_name]
            collection.insert_many(data_list)
            #collection2 = self.db[coll_name2]
            #collection.insert_many(data_list2)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (coll_name, e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % coll_name)
        try:
            if task_type in ["multiple_group", "two_group"]:
                main_collection = self.db[main_coll]
                box_data = {
                    "box_data": {"name": "function_name",
                                 "data": ["min", "q1", "median", "q3", "max"],
                                 "category": ""}
                }
                box_data_json = json.dumps(box_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id": table_id}, {"$set": {"main_id": table_id, "box_data": box_data_json}})
                column_data = {
                    "column_data": {"name": "function_name",
                                 "data": group_list,
                                 "category": ""}
                    }
                colum_data_json = json.dumps(column_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id": table_id}, {"$set": {"main_id": table_id, "box_data": box_data_json}})
        except Exception, e:
            self.bind_object.logger.error("更新%s信息出错:%s" % (main_coll, e))
            self.bind_object.set_error("更新信息出错")
        else:
            self.bind_object.logger.info("更新%s信息成功!" % main_coll)

    @report_check
    def add_species_difference_check_barplot(self, file_path, table_id, correlation_key=None,
                                             coll_name='multiple_group_bar', main_coll='multiple_group',task_type=None):
        stat_info, sort_list = self.cat_files(statfile=file_path, cifiles=None)
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        specimen_list = []
        for name in sort_list:
            data = [
                ("function_name", name),
                (correlation_key, table_id)
            ]
            for i in stat_info[name].keys():
                data.append((i, float(stat_info[name][i])))
                if i not in specimen_list:
                    specimen_list.append(i)
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
            if task_type in ["multiple_group", "two_group"]:
                main_collection = self.db[main_coll]
                single_column_data = {
                    "column_data": {"name": "function_name",
                                    "data": specimen_list,
                                    "category": ""}
                }
                single_column_data_json = json.dumps(single_column_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id": table_id},
                                       {"$set": {"main_id": table_id, "single_column_data": single_column_data_json}})
        except Exception, e:
            self.bind_object.logger.error("更新%s信息出错:%s" % (main_coll, e))
            self.bind_object.set_error("更新信息出错")
        else:
            self.bind_object.logger.info("更新%s信息成功!" % main_coll)

    @report_check
    def update_species_difference_check(self, table_id, statfile, cifile, test, main_coll='two_group'):
        collection = self.db[main_coll]
        with open(statfile, 'rb') as s, open(cifile, 'rb') as c:
            sinfo = s.readlines()
            meanlist = []
            lowci = []
            length = len(sinfo)
            if test == 'twogroup':
                for i in range(1, length):
                    line = sinfo[i].strip('\n').split('\t')
                    meanlist.append(float(line[1]))
                    meanlist.append(float(line[3]))
            else:
                for i in range(1, length):
                    line = sinfo[i].strip('\n').split('\t')
                    meanlist.append(float(line[1]))
                    meanlist.append(float(line[2]))
            max_mean = max(meanlist)
            cinfo = c.readlines()
            ci = []
            len_ci = len(cinfo)
            for i in range(1, len_ci):
                line_ci = cinfo[i].strip('\n').split('\t')
                ci.append(float(line_ci[3]) - float(line_ci[2]))
                lowci.append(float(line_ci[2]))
            max_ci = max(ci)
            min_low = min(lowci)
            if max_ci <= max_mean:
                n = 1
            else:
                n = round(max_ci / max_mean)
            l = round(abs(min_low / n) + max_mean + 3)
            collection.update({"_id": ObjectId(table_id)}, {"$set": {"n": n, "l": l}})

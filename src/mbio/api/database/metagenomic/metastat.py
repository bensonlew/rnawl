# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# last_modified 20171120
from biocluster.api.database.base import Base, report_check
import re
from mainapp.libs.param_pack import group_detail_sort
from mbio.packages.metagenomic.id_convert import name2id
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
import gridfs
import datetime
import os
import json
import pandas as pd


class Metastat(Base):
    def __init__(self, bind_object):
        super(Metastat, self).__init__(bind_object)
        self.name_dic = {}
        self._project_type = "metagenomic"
        self.main_col = "metastat"
        self.detail_col = "metastat_detail"
        self.plot_col = "metastat_plot"

    def change_main_col(self, main_col):
        self.main_col = main_col
        self.detail_col = main_col + "_detail"
        self.plot_col = main_col + "_plot"

    @report_check
    def add_metastat(self, anno_type, check_type, name=None, params=None):
        task_id = self.bind_object.sheet.id
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "status": "end",
            "created_ts": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            "anno_type": anno_type,
            "test_type": check_type,
        }
        if check_type == 'two_sample':
            insert_data["desc"] = "两样本间差异性检验"
            insert_data["name"] = name
        elif check_type == "two_group":
            insert_data["desc"] = "两组间差异性检验"
            insert_data["name"] = name
        elif check_type == "multiple":
            insert_data["desc"] = "多组间差异性检验"
            insert_data["name"] = name
        else:
            self.bind_object.set_error("不存在此分析类型：%s", variables=(check_type), code="52801601")
        collection = self.db[self.main_col]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def update_metastat(self, table_id, statfile, cifile, test_type):
        with open(statfile, 'rb') as s, open(cifile, 'rb') as c:
            sinfo = s.readlines()
            meanlist = []
            lowci = []
            length = len(sinfo)
            if test_type == "twogroup":
                for i in range(1, length):
                    line = sinfo[i].strip('\n').split('\t')
                    meanlist.append(float(line[1]))
                    meanlist.append(float(line[3]))
            else:
                for i in range(1, length):
                    line = sinfo[i].strip("\n").split('\t')
                    meanlist.append(float(line[1]))
                    meanlist.append(float(line[2]))
            max_mean = max(meanlist)
            cinfo = c.readlines()
            ci = []
            len_ci = len(cinfo)
            for i in range(1, len_ci):
                line_ci = cinfo[i].strip("\n").split("\t")
                ci.append(float(line_ci[3]) - float(line_ci[2]))
                lowci.append(float(line_ci[2]))
            max_ci = max(ci)
            min_low = min(lowci)
            if max_ci <= max_mean:
                n = 1
            else:
                n = round(max_ci / max_mean)
            l = round(abs(min_low / n) + max_mean + 3)
            collection = self.db[self.main_col]
            collection.update({"_id": ObjectId(table_id)}, {"$set": {"n": n, "l": l}})
        if test_type == "two_group":
            self.update_detail(table_id)

    @report_check
    def update_detail(self, table_id):
        collection = self.db[self.detail_col]
        collection_plot = self.db[self.plot_col]
        result = collection_plot.find({"metastat_id": ObjectId(table_id), "type": "box"})
        for one in result:
            collection_id = collection.find_one({"metastat_id": ObjectId(table_id), "species_name": one["species_name"]})["_id"]
            collection.update({"metastat_id": ObjectId(table_id), "species_name": one["species_name"], "_id": collection_id}, {"$set": {one["category_name"] + "-median": one["median"]}})

    @report_check
    def add_metastat_detail(self, statfile, cifiles, main_col=None, check_type=None, metastat_id=None, posthoc=None):
        if main_col:
            self.change_main_col(main_col)
        if check_type == "two_sample":
            stat_info, sort_list = self.cat_files(statfile, cifiles, convert_sample="forward", metastat_id=metastat_id)
        else:
            stat_info, sort_list = self.cat_files(statfile, cifiles, metastat_id=metastat_id)
        data_list = []
        log_index = 0
        new_insert = 0
        self.bind_object.logger.info("run first for loop")
        for name in sort_list:
            data = [
                ("metastat_id", ObjectId(metastat_id)),
                ("species_name", name),
            ]
            # self.bind_object.logger.info("run second for loop")
            total_mean =0
            for i in stat_info.columns:
                stat_info_line = stat_info.loc[name,]
                if i in ['pvalue', 'qvalue', 'corrected_pvalue', 'statistic', 'odds_ratio'] and stat_info_line[i] == 'NA':
                    if i == 'qvalue' or  i == 'corrected_pvalue':
                        data.append(('corrected_pvalue', stat_info_line[i]))
                    else:
                        data.append((i, stat_info_line[i]))
                elif re.search(r'_pvalue$', i) and posthoc in ['tukeykramer', 'gameshowell']:
                    data.append((i, stat_info_line[i]))
                elif i == 'qvalue':
                    data.append(('corrected_pvalue', float(stat_info_line[i])))
                elif re.search(r'-mean$', i):
                    total_mean += float(stat_info_line[i])
                    data.append((i, float(stat_info_line[i])))
                else:
                    data.append((i, float(stat_info_line[i])))
            data.append(("total_mean", total_mean))
            # self.bind_object.logger.info("second loop end")
            data_son = SON(data)
            data_list.append(data_son)
            log_index += 1
            if log_index % 200000 == 0 and log_index > 0:
                new_insert = 1
                self.bind_object.logger.info("treating %s lines" % log_index)
            if new_insert == 1 and log_index > 0:
                self.bind_object.logger.info("开始导入数据库前%s lines" % log_index)
                try:
                    collection = self.db[self.detail_col]
                    collection.insert_many(data_list)
                    new_insert = 0
                    data_list = []
                except Exception, e:
                    self.bind_object.logger.error("导入前%slines失败:%s" % (log_index, e))
                    self.bind_object.set_error("导入数据失败")
        # self.bind_object.logger.info("first loop end")
        self.bind_object.logger.info("导入剩余数据")
        try:
            collection = self.db[self.detail_col]
            collection.insert_many(data_list)
            if total_mean != 0:
                main_collection = self.db[self.main_col]
                main_collection.update({"_id": ObjectId(metastat_id)}, {"$set": {"is_new": 1}})
        except Exception, e:
            self.bind_object.logger.error("导入sg_species_difference_check_detail信息出错:%s" % e)
            self.bind_object.set_error("导入sg_species_difference信息出错", code="52801602")
        else:
            self.bind_object.logger.info("导入sg_species_difference_check_detail信息成功!")
        return metastat_id

    @report_check
    def add_metastat_bar(self, file_path, metastat_id):
        metastat_id = check_id(metastat_id, "metastat_id")
        data_list = []
        log_index = 0
        new_insert = 0
        stat_info, sort_list = self.cat_files(statfile=file_path, cifiles=None, convert_sample="afterward", metastat_id=metastat_id)
        for name in sort_list:
            data = [
                ("species_name", name),
                ("type", "bar"),
                ("metastat_id", ObjectId(metastat_id)),
            ]
            for i in stat_info.columns:
                stat_info_line = stat_info.loc[name,]
                data.append((i, float(stat_info_line[i])))
            data_son = SON(data)
            data_list.append(data_son)
            log_index += 1
            if log_index % 200000 == 0 and log_index > 0:
                new_insert = 1
                self.bind_object.logger.info("treating %s lines" % log_index)
            if new_insert == 1 and log_index > 0:
                self.bind_object.logger.info("开始导入数据库前%s lines" % log_index)
                try:
                    collection = self.db[self.plot_col]
                    collection.insert_many(data_list)
                    new_insert = 0
                    data_list = []
                except Exception,e:
                    self.bind_object.logger.error("导入前%s")
        self.bind_object.logger.info("导入剩余数据")
        try:
            collection = self.db[self.plot_col]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入单物种/功能柱状图 %s信息出错：%s" % (file_path, e))
            self.bind_object.set_error("导入单物种/功能柱状图信息出错", code="52801603")
        else:
            self.bind_object.logger.info("导入单物种柱状图信息成功！")

    @report_check
    def add_metastat_box(self, file_path, metastat_id):
        metastat_id = check_id(metastat_id, "metastat_id")
        data_list = []
        log_index = 0
        new_insert = 0
        with open(file_path, 'rb') as r:
            l = r.readline().strip("\n")
            group_list = re.findall(r'min\((.*?)\)', l)
            while True:
                line = r.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                i = 1
                for name in group_list:
                    species_name = line_data[0].split(';')[-1]
                    data = [("metastat_id", ObjectId(metastat_id)), ("species_name", species_name), ("type", "box")]
                    data.append(("category_name", name))
                    data.append(("min", float(line_data[i])))
                    data.append(("q1", float(line_data[i + 1])))
                    data.append(("median", float(line_data[i + 2])))
                    data.append(("q3", float(line_data[i + 3])))
                    data.append(("max", float(line_data[i + 4])))
                    i += 5
                    data_son = SON(data)
                    data_list.append(data_son)
                    log_index += 1
                    if log_index % 200000 == 0 and log_index > 0:
                        new_insert = 1
                        self.bind_object.logger.info("treating %s lines" % log_index)
                    if new_insert == 1 and log_index > 0:
                        self.bind_object.logger.info("开始导入数据前%s lines" % log_index)
                        try:
                            collection = self.db[self.plot_col]
                            collection.insert_many(data_list)
                            new_insert = 0
                            data_list = []
                        except Exception, e:
                            self.bind_object.logger.error("导入前%s lines失败" % (log_index, e))
                            self.bind_object.set_error("导入数据失败")
        self.bind_object.logger.info("导入剩余数据")
        try:
            collection = self.db[self.plot_col]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错%s" % (file_path, e))
            self.bind_object.set_error("导入metastat_plot信息出错", code="52801604")
        else:
            self.bind_object.logger.info("导入%s信息成功！" % file_path)
        return data_list

    def cat_files(self, statfile, cifiles=None, convert_sample=False, metastat_id=None):
        """
        将组间差异比较的统计结果文件与posthoc结果文件cat成一个文件
        cifiles为posthoc文件路径的列表
        """
        if metastat_id == None:
            self.bind_object.set_error("必须输入metastat_id", code="52801605")
        main_task_id = self.db[self.main_col].find_one({'_id': ObjectId(metastat_id)})["task_id"]
        self.bind_object.logger.info("task_id is : %s" % main_task_id)
        self.name_dic = name2id(main_task_id, type='task')
        sdata = pd.read_table(statfile,index_col=0)
        self.bind_object.logger.info("read stat file over")
        rename_dict = dict()
        for i in sdata.columns:
            if '-' in i:
                str1, str2 = i.split('-')
                if convert_sample == "forward":
                    str1 = self.name_dic.get(str1, str1)
                elif convert_sample == "afterward":
                    str2 = self.name_dic.get(str2, str2)
                rename_dict[i] = '-'.join([str1,str2])
        sdata.rename(columns=rename_dict, inplace=True)
        self.bind_object.logger.info("rename over")
        sort_list = sdata.index.tolist()
        self.bind_object.logger.info("read cifiles...")
        if cifiles:
            for f in cifiles:
                self.bind_object.logger.info("%s" % f)
                cidata = pd.read_table(f, index_col=0)
                sdata = pd.merge(sdata, cidata, left_index=True, right_index=True, how="left")
        if cifiles and len(cifiles) > 1:
            self.bind_object.logger.info("add nl")
            tmp_data = sdata.apply(self.get_nl, axis=1)
            self.bind_object.logger.info(tmp_data.head())
            sdata["n"] = tmp_data["n"]
            sdata["l"] = tmp_data["l"]
        self.bind_object.logger.info("cat files over")
        return sdata, sort_list

    def get_nl(self, df):
        mean = []
        low_ci = []
        ci = []
        for i in df.index:
            if re.search(r'mean$', i):
                mean.append(float(df[i]))
            elif re.search(r'lowerCI$', i):
                group = i.split('_lowerCI')[0]
                up = group + '_upperCI'
                low_ci.append(float(df[i]))
                ci.append(float(df[up] - float(df[i])))
        max_mean = max(mean)
        max_ci = max(ci)
        min_ci = min(low_ci)
        if max_ci < max_mean:
            n = 1
        else:
            n = round(max_ci / max_mean)
        l = round(abs(min_ci / n) + max_mean + 3)
        return_value = pd.Series({"n":n, "l": l}, name=df.name)
        return return_value

def get_anno_info(anno_id):
    anno_id = int(anno_id)
    hash = {
        1: {"db": "NR", "level": "Domain"},
        2: {"db": "NR", "level": "Kingdom"},
        3: {"db": "NR", "level": "Phylum"},
        4: {"db": "NR", "level": "Class"},
        5: {"db": "NR", "level": "Order"},
        6: {"db": "NR", "level": "Family"},
        7: {"db": "NR", "level": "Genus"},
        8: {"db": "NR", "level": "Species"},
        9: {"db": "COG", "level": "Category"},
        10: {"db": "COG", "level": "Function"},
        11: {"db": "COG", "level": "NOG"},
        12: {"db": "KEGG", "level": "Pathway.Level1"},
        13: {"db": "KEGG", "level": "Pathway.Level2"},
        14: {"db": "KEGG", "level": "Pathway.Level3"},
        15: {"db": "KEGG", "level": "Module"},
        16: {"db": "KEGG", "level": "Enzyme"},
        17: {"db": "KEGG", "level": "KO"},
        18: {"db": "CAZY", "level": "Class"},
        19: {"db": "CAZY", "level": "Family"},
        20: {"db": "ARDB", "level": "Class"},
        21: {"db": "ARDB", "level": "Type"},
        22: {"db": "ARDB", "level": "Antibiotic.type"},
        23: {"db": "ARDB", "level": "ARG"},
        24: {"db": "CARD", "level": "Class"},
        25: {"db": "CARD", "level": "ARO"},
        26: {"db": "VFDB", "level": "Level1"},
        27: {"db": "VFDB", "level": "Level2"},
        28: {"db": "VFDB", "level": "Vfs"},
    }
    return hash[anno_id]["db"], hash[anno_id]["level"]


def check_id(objectid, id_name=""):
    if not isinstance(objectid, ObjectId):
        if isinstance(objectid, StringTypes):
            objectid = ObjectId(objectid)
        else:
            self.bind_object.set_error("%s必须为ObjectId对象或其对应的字符串!", variables=(id_name), code="52801606")
    return objectid

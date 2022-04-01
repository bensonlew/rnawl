# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.api.database.base import Base, report_check
import re
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
import gridfs
import datetime
import os
from mainapp.libs.param_pack import group_detail_sort, param_pack
import json
# from biocluster.config import Config


class StatTest(Base):
    def __init__(self, bind_object):
        super(StatTest, self).__init__(bind_object)
        # self._db_name = Config().MONGODB
        self._project_type = "meta"

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
    def add_species_difference_check_detail(self, statfile, cifiles, level=None, check_type=None, params=None, category_name=None, table_id=None, group_id=None, from_otu_table=None, major=False, posthoc=None):
        if major:
            table_id = self.create_species_difference_check(level, check_type, params, category_name, group_id, from_otu_table)
        else:
            if table_id is None:
                self.bind_object.set_error("major为False时需提供table_id!", code="51007401")
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
                else:
                    self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51007402")
        stat_info, sort_list = self.cat_files(statfile, cifiles)
        self.bind_object.logger.info("stat_info: %s" % stat_info)
        self.bind_object.logger.info("sort_list: %s" % sort_list)
        data_list = []
        for name in sort_list:
            data = [
                ("species_check_id", table_id),
                ("full_species_name", name),
                ("species_name", name.split('; ')[-1])
            ]
            for i in stat_info[name].keys():
                if i in ['pvalue', 'qvalue', 'corrected_pvalue', 'statistic', 'odds_ratio'] and stat_info[name][i] == 'NA':
                    if i == 'qvalue' or  i == 'corrected_pvalue':
                        data.append(('corrected_pvalue', stat_info[name][i]))
                    else:
                        data.append((i, stat_info[name][i]))
                elif re.search(r'_pvalue$', i) and posthoc in ['tukeykramer', 'gameshowell']:
                    data.append((i, stat_info[name][i]))
                elif i == 'qvalue':
                    data.append(('corrected_pvalue', float(stat_info[name][i])))
                else:
                    data.append((i, float(stat_info[name][i])))
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["sg_species_difference_check_detail"]
            collection.insert_many(data_list)

            main_collection = self.db["sg_species_difference_check"]
            #main_collection.update({"_id": table_id},{"$set": { "main_id": table_id}})

        except Exception, e:
            self.bind_object.logger.error("导入sg_species_difference_check_detail信息出错:%s" % e)
            self.bind_object.set_error("导入sg_species_difference_check_detail信息出错", code="51007403")
        else:
            self.bind_object.logger.info("导入sg_species_difference_check_detail信息成功!")
        return table_id

    @report_check
    def add_species_difference_check_barplot(self, file_path, table_id):
        stat_info, sort_list = self.cat_files(statfile=file_path, cifiles=None)
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51007402")
        data_list = []
        for name in sort_list:
            data = [
                ("species_name", name.split('; ')[-1]),
                ("type", "bar"),
                ("species_check_id", table_id)
            ]
            for i in stat_info[name].keys():
                data.append((i, float(stat_info[name][i])))
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["sg_species_difference_check_boxplot"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_species_difference_check"]
            #main_collection.update({"_id": table_id},{"$set": { "main_id": table_id}})
        except Exception, e:
            self.bind_object.logger.error("导入单物种柱状图信息出错:%s" % e)
            self.bind_object.set_error("导入单物种柱状图信息出错", code="51007404")
        else:
            self.bind_object.logger.info("导入单物种柱状图信息成功!")

    @report_check
    def add_species_difference_check_boxplot(self, file_path, table_id):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51007402")
        data_list = []
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
                    data = [("species_check_id", table_id), ("species_name", line_data[0].split('; ')[-1]), ("type", "box")]
                    data.append(("category_name", name))
                    data.append(("min", float(line_data[i])))
                    data.append(("q1", float(line_data[i + 1])))
                    data.append(("median", float(line_data[i + 2])))
                    data.append(("q3", float(line_data[i + 3])))
                    data.append(("max", float(line_data[i + 4])))
                    i += 5
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["sg_species_difference_check_boxplot"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错", code="51007405")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    @report_check
    def add_species_difference_lefse_detail(self, file_path, params=None, group_id=None, from_otu_table=None, table_id=None, major=False):
        if major:
            table_id = self.create_species_difference_lefse(params, group_id, from_otu_table)
        else:
            if table_id is None:
                self.bind_object.set_error("major为False时需提供table_id!", code="51007401")
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
                else:
                    self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51007402")
        data_list = []
        with open(file_path, 'rb') as r:
            r.readline()
            is_use = False
            for line in r:
                line = line.strip('\n')
                line_data = line.split('\t')
                if line_data[4] != '-':
                    line_data[4] = float(line_data[4])
                data = [
                    ("species_lefse_id", table_id),
                    ("species_name", line_data[0]),
                    ("category_name", line_data[2]),
                    ("median", float(line_data[1])), ("lda", line_data[3]),
                    ("pvalue", line_data[4])
                ]
                data_son = SON(data)
                data_list.append(data_son)
                if not line_data[2] == '-':
                    is_use = True
        if is_use:
            coll_main = self.db["sg_species_difference_lefse"]
            coll_main.update({"_id": ObjectId(table_id)}, {"$set": {"lda_png_id": 'is_useful', "lda_cladogram_id": 'is_useful'}})
        try:
            collection = self.db["sg_species_difference_lefse_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_species_difference_lefse"]
            #main_collection.update({"_id": table_id},{"$set": { "main_id": table_id}})

        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错", code="51007405")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list, table_id

    @report_check
    def update_species_difference_lefse(self, lda_png_path, lda_cladogram_path, table_id):
        size = os.path.getsize(lda_png_path)
        if size == 0:
            pass
        else:
            collection = self.db["sg_species_difference_lefse"]
            fs = gridfs.GridFS(self.db)
            ldaid = fs.put(open(lda_png_path, 'rb'))
            cladogramid = fs.put(open(lda_cladogram_path, 'rb'))
            try:
                collection.update({"_id": ObjectId(table_id)}, {"$set": {"lda_png_id": ldaid, "lda_cladogram_id": cladogramid}})
            except Exception, e:
                self.bind_object.logger.error("导入%s和%s信息出错:%s" % (lda_png_path, lda_cladogram_path, e))
                self.bind_object.set_error("导入信息出错", code="51007405")
            else:
                self.bind_object.logger.info("导入%s和%s信息成功!" % (lda_png_path, lda_cladogram_path))

    @report_check
    def update_species_difference_check(self, table_id, statfile, cifile, test):
        collection = self.db["sg_species_difference_check"]
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

    @report_check
    def create_species_difference_check(self, level, check_type, params, category_name, group_id=None, from_otu_table=None, name=None):
        if from_otu_table and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51007406")
        if group_id and not isinstance(group_id, ObjectId):
            if isinstance(group_id, StringTypes):
                group_id = ObjectId(group_id)
            else:
                self.bind_object.set_error("group_id必须为ObjectId对象或其对应的字符串!", code="51007407")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "组间差异性检验"
        if check_type == 'tow_sample':
            insert_data = {
                "type": check_type,
                "project_sn": project_sn,
                "task_id": task_id,
                "otu_id": from_otu_table,
                "name": self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else "DiffStatTwoSample_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
                "level_id": int(level),
                "params": params,
                "desc": desc,
                "status": "end",
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
        else:
            insert_data = {
                "type": check_type,
                "project_sn": project_sn,
                "task_id": task_id,
                "otu_id": from_otu_table,
                "group_id": group_id,
                "level_id": int(level),
                "params": params,
                "desc": desc,
                "status": "end",
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "category_name": category_name
            }
            if check_type == 'two_group':
                insert_data['name'] = self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else "DiffStatTwoGroup_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            else:
                insert_data['name'] = self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else "DiffStatMultiple_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        collection = self.db["sg_species_difference_check"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def create_species_difference_lefse(self, params, group_id=0, from_otu_table=0, name=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51007406")
        if group_id != 0 and not isinstance(group_id, ObjectId):
            if isinstance(group_id, StringTypes):
                group_id = ObjectId(group_id)
            else:
                self.bind_object.set_error("group_id必须为ObjectId对象或其对应的字符串!", code="51007407")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "lefse分析"
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "otu_id": from_otu_table,
            "group_id": group_id,
            "name": name if name else "LEfSe_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
            "params": params,
            "lda_cladogram_id": "",
            "lda_png_id": "",
            "desc": desc,
            "status": "end",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_species_difference_lefse"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def creat_multi_table(self,from_otu_table,group_id,params,level_id,category_name,name=None,status=None,spname_spid=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51004801")
        if group_id not in ["all", "All", "ALL"]:
            if not isinstance(group_id, ObjectId):
                if isinstance(group_id, StringTypes):
                    group_id = ObjectId(group_id)
                else:
                    self.bind_object.set_error("group_id必须为ObjectId对象或其对应的字符串!", code="51004802")
        if not status:
            status = "end"
        if spname_spid:
            my_params = json.loads(params)
            group_detail = {'All': [str(i) for i in spname_spid.values()]}
            my_params['group_detail'] = group_detail_sort(group_detail)
            my_params['otu_id'] = str(from_otu_table)
            params = param_pack(my_params)
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("sg_otu表找不到相应记录", code="51004803")
        project_sn = result['project_sn']
        task_id = result['task_id']
        insert_data = {
            "type": 'multiple',
            "project_sn": project_sn,
            "task_id": task_id,
            "otu_id": from_otu_table,
            "group_id": group_id,
            "name": name,
            "level_id": level_id,
            "params": params,
            "status": "end",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "category_name": category_name
        }

        collection = self.db["sg_species_difference_check"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def creat_two_table(self,from_otu_table,group_id,params,level_id,category_name,name=None,status=None,spname_spid=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51004801")
        if group_id not in ["all", "All", "ALL"]:
            if not isinstance(group_id, ObjectId):
                if isinstance(group_id, StringTypes):
                    group_id = ObjectId(group_id)
                else:
                    self.bind_object.set_error("group_id必须为ObjectId对象或其对应的字符串!", code="51004802")
        if not status:
            status = "end"
        if spname_spid:
            my_params = json.loads(params)
            group_detail = {'All': [str(i) for i in spname_spid.values()]}
            my_params['group_detail'] = group_detail_sort(group_detail)
            my_params['otu_id'] = str(from_otu_table)
            params = param_pack(my_params)
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("sg_otu表找不到相应记录", code="51004803")
        project_sn = result['project_sn']
        task_id = result['task_id']
        insert_data = {
            "type": 'two_group',
            "project_sn": project_sn,
            "task_id": task_id,
            "otu_id": from_otu_table,
            "group_id": group_id,
            "name": name,
            "level_id": level_id,
            "params": params,
            "status": "end",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "category_name": category_name
        }
        collection = self.db["sg_species_difference_check"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id



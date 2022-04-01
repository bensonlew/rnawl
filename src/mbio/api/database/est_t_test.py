# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# last modified guhaidong 20171115
from biocluster.api.database.base import Base, report_check
import re
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
import datetime
import json
import os
# from biocluster.config import Config


class EstTTest(Base):
    def __init__(self, bind_object):
        super(EstTTest, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB

    @report_check
    def get_another_name(self, name, group_list):
        another = ''
        for n in group_list:
            if n == name:
                pass
            else:
                another = n
        return another

    @report_check
    def add_est_t_test_detail(self, file_path, table_id, alpha_id, group_name=None):  #guanqing.zuo 20180510 add alpha_id
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51003001")
        # by houshuang 20191011 增加箱式图数据>>>
        box_file = re.sub("est_result", "est_barplot_result", file_path)
        self.bind_object.logger.info("boxfile_path:{}".format(box_file))
        dict = {}
        with open(box_file, 'rb') as r:
            l = r.readline().strip('\n')
            group_list = re.findall(r'min\((.*?)\)', l)
            while True:
                line = r.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                dict[line_data[0]] = {}
                i = 1
                for name in group_list:
                    value = [str(line_data[i]), str(line_data[i + 1]), str(line_data[i + 2]),
                             str(line_data[i + 3]), str(line_data[i + 4])]  # min,q1,median,q3,max
                    dict[line_data[0]][name] = ",".join(value)
                    i += 5
        self.bind_object.logger.info("name to box_value:{}".format(dict))
        # <<<
        data_list = []
        with open(file_path, 'rb') as r:
            l = r.readline().strip('\n')
            # group_list = re.findall(r'mean\((.*?)\)', l)
            group_list = re.findall(r'(\S+)-mean', l)  # by houshuang 20191011 更换成两组比较的脚本，结果文件不同
            while True:
                line = r.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                length = len(line_data)
                i = 1
                for name in group_list:
                    data = [("alpha_ttest_id", table_id), ("index_type", line_data[0]), ("qvalue", line_data[length-1]), ("pvalue", line_data[length-2]),("statistic", line_data[length-3])]
                    data.append(("category_name", name))
                    data.append(("compare_name", self.get_another_name(name, group_list)))
                    data.append(("mean", str('%0.5g' % float(line_data[i]))))
                    data.append(("sd", str('%0.5g' % float(line_data[i+1]))))
                    data.append(("box_value", dict[line_data[0]][name]))  # by houshuang 20191011，增加箱式图数据
                    i += 2
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["sg_alpha_ttest_detail"]
            collection.insert_many(data_list)
            collection_alpha = self.db["sg_alpha_diversity"]  # find alpha_diversity_id by guhaidong 20171115
            task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])  # get task_id by guhaidong 20171115
            result = collection_alpha.find_one({"task_id": task_id})
            alpha_diversity_id = result['_id']
            collection_main = self.db["sg_alpha_ttest"]
          #  collection_main.update({"_id": ObjectId(table_id), "alpha_diversity_id": alpha_diversity_id}, {"$set": {"compare_column": group_name}})
            collection_main.update({"_id": ObjectId(table_id), "alpha_diversity_id": ObjectId(alpha_id)}, {"$set": {"compare_column": group_name, "has_box": 1}})#guanqing.zou 20180510
            #collection_main.update({"_id": ObjectId(table_id)}, {"$set": {"main_id": ObjectId(table_id)}})
            self.bind_object.logger.info("is running add_est_t_test_detail")
            # add alpha_diversity_id by guhaidong 20171115
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)

    @report_check
    def add_est_t_test_collection(self, params, group_id, from_est_id=0, name=None, group_name=None):
        if isinstance(from_est_id, StringTypes):
            from_est_id = ObjectId(from_est_id)
        else:
            self.bind_object.set_error("est_id必须为ObjectId对象或其对应的字符串!", code="51003002")
        if group_id == "all":
            group_id = group_id
        else:
            group_id = ObjectId(group_id)
        task_id = self.bind_object.sheet.id  # add task_id by guhaidong 20171115
        task_id = "_".join(task_id.split('_')[0:2])  # get task_id by guhaidong 20171115
        collection = self.db["sg_alpha_diversity"]
        result = collection.find_one({"_id": from_est_id, "task_id": task_id})
        project_sn = result['project_sn']
        # task_id = result['task_id']
        otu_id = result['otu_id']
        level_id = result['level_id']
        desc = ""
        insert_data = {
                "project_sn": project_sn,
                "task_id": task_id,
                "otu_id": otu_id,
                "alpha_diversity_id": from_est_id,
                "name": self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else "多样性指数T检验结果表",
                "level_id": int(level_id),
                "group_id": group_id,
                "compare_column": group_name,
                "status": "end",
                "desc": desc,
                "params": params,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
        collection = self.db["sg_alpha_ttest"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        #collection.update({"_id": inserted_id}, {"$set": {"main_id": inserted_id}})
        self.bind_object.logger.info("is running add_est_t_test_collection")
        return inserted_id

    def add_est_table(self, est_path, level, est_id=None, group_name=None, otu_id=None, task_id=None, name=None, params=None,compare_column=None):
        if level not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level)
            self.bind_object.set_error("level参数不在允许的范围内", code="51005503")
        if task_id is None:
            task_id = self.bind_object.sheet.id
        if otu_id:
            if not isinstance(otu_id, ObjectId):
                otu_id = ObjectId(otu_id)
            params['otu_id'] = str(otu_id)  # otu_id在再metabase中不可用

        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "otu_id": otu_id,
            "alpha_diversity_id": est_id,
            "name": name if name else "EstTTest_Origin",
            "level_id": level,
            "status": "end",
            "desc": "Job has been finished",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        # if params is not None:
            # insert_data['params'] = json.dumps(params, sort_keys=True, separators=(',', ':'))
        collection = self.db["sg_alpha_ttest"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        for f in os.listdir(est_path):
            if re.search("est_result", f):
                self.add_est_t_test_detail(os.path.join(est_path, f), inserted_id, est_id, group_name=compare_column)
        return inserted_id
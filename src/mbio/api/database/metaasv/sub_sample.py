# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import re
import datetime
import json
from bson.objectid import ObjectId
from types import StringTypes


class SubSample(Base):
    def __init__(self, bind_object):
        super(SubSample, self).__init__(bind_object)
        self._project_type = 'metaasv'
        self.name_id = dict()  # otu表中样本名和id对照的字典
        self.otu_rep = dict()  # o
        self.main_task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])

    def add_sg_otu(self, params, from_otu_table=0, name=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!")

        task_id = self.main_task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            'task_id': task_id,
            'from_id': str(from_otu_table),
            'name': self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else name,
            "params": params,
            'status': 'end',
            "level_id": json.dumps([9]),
            'desc': 'ASV主表',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["asv"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def _get_name_id(self, from_otu_id):
        collection = self.db['asv_specimen']
        results = collection.find({"asv_id": from_otu_id})
        if not results.count():
            self.bind_object.logger.error("asv_id:{}未在asv_specimen表里找到相应的记录".format(from_otu_id))
            self.bind_object.set_error("asv_specimen表找不到相应记录")
        sp_ids = list()
        for result in results:
            sp_ids.append(result['specimen_id'])
        collection = self.db['specimen_detail']
        for id_ in sp_ids:
            result = collection.find_one({"_id": id_})
            if not result:
                self.bind_object.logger.error("意外错误， id: {}在asv_specimen表中找到，但未在specimen_detail表中出现")
                self.bind_object.set_error("asv_specimen表找不到相应记录")
            self.name_id[result["specimen"]] = id_

    def _get_sapmple_id(self, task_id):
        collection = self.db['specimen']
        results = collection.find_one({"task_id": task_id})
        # if not results.count():
        #     self.bind_object.logger.error("asv_id:{}未在specimen表里找到相应的记录".format(task_id))
        #     self.bind_object.set_error("asv_specimen表找不到相应记录")
        main_id = results["_id"]
        collection_detail = self.db['specimen_detail']
        results = collection_detail.find({"specimen_id": main_id})
        specimen_dic = {}
        for result in results:
            new_specimen = result['specimen']
            specimen_id = result["_id"]
            specimen_dic[new_specimen] = specimen_id
        return specimen_dic

    def _get_task_info(self, asv_id):
        collection = self.db['asv']
        result = collection.find_one({'_id': asv_id})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在asv表里找到相应的记录".format(asv_id))
            self.bind_object.set_error("asv_id表找不到相应记录")
        self.project_sn = result['project_sn']
        self.task_id = result['task_id']

    def _prepare_otu_rep(self, asv_id):
        """
        实际运行的时候发现对每一行(即每一个otu)都去数据库里查询一次，并获取otu_rep的时候，效率非常低，需要很长的时间， 因此，需要对mongo做一次查询， 将属于一个otu_id的otu_rep全部去读出来， 放到内存当中， 以提高效率
        """
        self.bind_object.logger.info("开始依据asv_id: {}查询所有的代表序列".format(asv_id))
        collection = self.db["asv_detail"]
        results = collection.find({"asv_id": asv_id})
        for result in results:
            self.otu_rep[result['asv']] = result["asv_rep"]
        self.bind_object.logger.info("代表序列查询完毕")

    @report_check
    def add_sg_otu_detail(self, file_path, from_otu_id, new_otu_id, new_samples=False):
        if not isinstance(from_otu_id, ObjectId):
            if isinstance(from_otu_id, StringTypes):
                from_otu_id = ObjectId(from_otu_id)
            else:
                self.bind_object.set_error("from_otu_id必须为ObjectId对象或其对应的字符串!")
        if not isinstance(new_otu_id, ObjectId):
            if isinstance(new_otu_id, StringTypes):
                new_otu_id = ObjectId(new_otu_id)
            else:
                self.bind_object.set_error("new_otu_id必须为ObjectId对象或其对应的字符串!")
        self._get_task_info(new_otu_id)
        # self._get_name_id(from_otu_id)
        sample_dict = self._get_sapmple_id(self.task_id)

        # 导入sg_otu_detail表
        self.bind_object.logger.info("开始导入asv_detail表")
        self.bind_object.logger.info("file_path:%s"%file_path)
        self._prepare_otu_rep(from_otu_id)
        insert_data = list()
        with open(file_path, 'rb') as r:
            head = r.next().strip('\r\n')
            head = re.split('\t', head)
            if head[-1] == "taxonomy":
                new_head = head[1:-1]
            else:
                new_head = head[1:]
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                otu_detail = dict()
                # OTU表可能有taxonomy列， 也可能没有， 需要适应
                if len(re.split(";", line[0])) > 1:
                    sample_num = line[1:]
                    otu = re.split(";", line[0])[-1]
                    classify_list = re.split(r"\s*;\s*", line[0])
                else:
                    sample_num = line[1:-1]
                    otu = line[0]
                    otu_detail['asv'] = line[0]
                    classify_list = re.split(r"\s*;\s*", line[-1])

                otu_detail['asv_id'] = new_otu_id
                for cf in classify_list:
                    if cf != "":
                        otu_detail[cf[0:3].lower()] = cf
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = int(sample_num[i])

                if otu not in self.otu_rep:
                    self.bind_object.logger.error("意外错误，asv_id: {}和asv: {}在asv_detail表里未找到".format(from_otu_id, line[0]))
                    self.bind_object.set_error("asv_detail表中找不到相应记录")
                otu_detail['asv_rep'] = self.otu_rep[otu]
                otu_detail['task_id'] = self.task_id
                insert_data.append(otu_detail)
        try:
            collection = self.db['asv_detail']
            collection.insert_many(insert_data)

            # main_collection = self.db["asv"]
            #main_collection.update({"_id": new_otu_id},{"$set": { "main_id": new_otu_id}})

        except Exception as e:
            self.bind_object.logger.error("导入asv_detail表格信息出错:{}".format(e))
            self.bind_object.set_error("导入asv_detail表格信息出错")
        else:
            self.bind_object.logger.info("导入asv_detail表格成功")
        # 导入sg_otu_specimen表
        self.bind_object.logger.info("开始导入asv_specimen表")
        # if not new_samples:
        #     insert_data = list()
        #     for sp in new_head:
        #         my_data = dict()
        #         my_data['asv_id'] = new_otu_id
        #         my_data["specimen_id"] = self.name_id[sp]
        #         insert_data.append(my_data)
        #     collection = self.db['asv_specimen']
        #     collection.insert_many(insert_data)
        if not new_samples:
            insert_data = list()
            for sp in new_head:
                my_data = dict()
                my_data['asv_id'] = new_otu_id
                if sp in sample_dict:
                    my_data["specimen_id"] = sample_dict[sp]
                    insert_data.append(my_data)
            collection = self.db['asv_specimen']
            collection.insert_many(insert_data)

    @report_check
    def add_sg_otu_detail_level(self, otu_path, from_otu_table, level, database=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!")
        collection = self.db["asv"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在asv表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("asv表找不到相应记录")
        project_sn = result['project_sn']
        task_id = result['task_id']
        # covered_level = list()
        # if "level_id" in result:
        #     covered_level = result["level_id"]
        #     covered_level.append(int(level))
        # else:
        #     covered_level.append(int(level))
        # covered_level = list(set(covered_level))
        # covered_level.sort()
        # result["level_id"] = json.dumps(covered_level)
        collection.update({"_id": from_otu_table}, {"$set": result}, upsert=False)
        insert_data = list()
        with open(otu_path, 'rb') as r:
            head = r.next().strip('\r\n')
            head = re.split('\t', head)
            if head[-1] == "taxonomy":
                new_head = head[1:-1]
            else:
                new_head = head[1:]
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                otu_detail = dict()

                if len(re.split("; ", line[0])) > 1:
                    sample_num = line[1:-2]
                    classify_list = re.split(r"\s*;\s*", line[0])
                elif len(re.split(";", line[0])) > 1:
                    sample_num = line[1:-2]
                    classify_list = re.split(r"\s*;\s*", line[0])
                else:
                    sample_num = line[1:-2]
                    otu_detail['asv'] = line[0]
                    classify_list = re.split(r"\s*;\s*", line[-1])

                otu_detail['asv_id'] = from_otu_table
                otu_detail['project_sn'] = project_sn
                otu_detail['task_id'] = task_id
                otu_detail["level_id"] = int(level)
                for cf in classify_list:
                    if cf != "":
                        otu_detail[cf[0:3].lower()] = cf
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = int(sample_num[i])
                    # count += int(sample_num[i])
                otu_detail["total"] = int(line[-2])
                otu_detail["percent"] = float(line[-1])
                insert_data.append(otu_detail)

        # 页面概览需要
        d_dict = []
        k_dict = []
        p_dict = []
        c_dict = []
        o_dict = []
        f_dict = []
        g_dict = []
        s_dict = []
        otu = []
        with open(otu_path, 'rb') as r:
            data = r.readlines()
            for i in data[1:]:
                classify_list = i.strip().split("\t")[0].split(";")
                for cf in classify_list:
                    if cf[0] == "d":
                        d_dict.append(cf)
                    elif cf[0] == "k":
                        k_dict.append(cf)
                    elif cf[0] == "p":
                        p_dict.append(cf)
                    elif cf[0] == "c":
                        c_dict.append(cf)
                    elif cf[0] == "o":
                        o_dict.append(cf)
                    elif cf[0] == "f":
                        f_dict.append(cf)
                    elif cf[0] == "g":
                        g_dict.append(cf)
                    elif cf[0] == "s":
                        s_dict.append(cf)
        table_dict = {"Domain": str(len(list(set(d_dict)))), "kingdom": str(len(list(set(k_dict)))), "Phylum": str(len(list(set(p_dict)))),
                      "Class": str(len(list(set(c_dict)))), "Order": str(len(list(set(o_dict)))), "Family": str(len(list(set(f_dict)))),
                      "Genus": str(len(list(set(g_dict)))), "Species": str(len(list(set(s_dict)))), "OTU": str(len(data) - 1)}
        table_info = json.dumps(table_dict, sort_keys=False, separators=(',', ':'))

        try:
            collection = self.db['asv_detail_level']
            collection.insert_many(insert_data)
            task = self.db["sg_task"]
            task.update({"task_id": task_id}, {"$set": {"species_stat": table_info}})
        except Exception as e:
            self.bind_object.logger.error("导入asv_detail_level表格失败：{}".format(e))
            self.bind_object.set_error("导入asv_detail_level表格失败")
        else:
            self.bind_object.logger.info("导入asv_detail_level表格成功")
        try:
            table_data = {
                "table_data": ["specimen","read_number"]
                }
            table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
            if self.bind_object.option("pipeline") == "functional_gene":
                if database:
                    origin_settled_params = self.get_task_info_from_task()
                    settled_params = {
                        "software": "Qiime2 (v2020.2) "+self.bind_object.option("denoise_method"),
                        "database": origin_settled_params["fungene_database"],
                        "anno_method": origin_settled_params["seq_style"],
                        "acid_length": self.bind_object.option("acid_length"),
                        "seq_identity": self.bind_object.option("seq_identity"),
                    }
                    if origin_settled_params["anno_method"] in ["vsearch", "multi_blast", "blast"]:
                        settled_params["identity"] = origin_settled_params["fungene_identity"]
                        settled_params["coverage"] = origin_settled_params["fungene_coverage"]
                    else:
                        settled_params["confidence"] = origin_settled_params["fungene_confidence"]
                    settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
                else:
                    settled_params = {
                        "software": "Qiime2 (v2020.2) "+self.bind_object.option("denoise_method"),
                        "database": self.bind_object.option("fungene_database"),
                        "anno_method": self.bind_object.option("seq_style"),
                        "acid_length": self.bind_object.option("acid_length"),
                        "seq_identity": self.bind_object.option("seq_identity"),
                    }
                    if self.bind_object.option("anno_method") in ["vsearch", "multi_blast", "blast"]:
                        settled_params["identity"] = self.bind_object.option("fungene_identity")
                        settled_params["coverage"] = self.bind_object.option("fungene_coverage")
                    else:
                        settled_params["confidence"] = self.bind_object.option("fungene_confidence")
                    settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            else:
                if database:
                    origin_settled_params = self.get_task_info_from_task()
                    settled_params = {
                        "software": "Qiime2 (v2020.2)",
                        "database": origin_settled_params["database"],
                        "anno_method": origin_settled_params["anno_method"]
                    }
                    if origin_settled_params["anno_method"] in ["vsearch", "multi_blast", "blast"]:
                        settled_params["identity"] = origin_settled_params["identity"]
                        settled_params["coverage"] = origin_settled_params["coverage"]
                    else:
                        settled_params["confidence"] = origin_settled_params["confidence"]
                    settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
                else:
                    settled_params = {
                        "software": "Qiime2 (v2020.2)",
                        "database": self.bind_object.option("database"),
                        "anno_method": self.bind_object.option("anno_method")
                    }
                    if self.bind_object.option("anno_method") in ["vsearch", "multi_blast", "blast"]:
                        settled_params["identity"] = self.bind_object.option("identity")
                        settled_params["coverage"] = self.bind_object.option("coverage")
                    else:
                        settled_params["confidence"] = self.bind_object.option("confidence")
                    settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            main_collection = self.db['asv']
            main_collection.update_one({"_id": from_otu_table}, {"$set": {"table_data": table_data_json,
                                                                          "main_id": from_otu_table,
                                                                          "settled_params": settled_params_json}})
        except Exception as e:
            self.bind_object.logger.error("导入asv_detail_level表格失败：{}".format(e))
            self.bind_object.set_error("导入asv_detail_level表格失败")

    def add_otu_detail(self, otu_path, from_otu_table, level):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                raise Exception("from_otu_table必须为ObjectId对象或其对应的字符串!")
        collection = self.db["asv"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            raise Exception("无法根据传入的_id:{}在asv表里找到相应的记录".format(str(from_otu_table)))
        project_sn = result['project_sn']
        task_id = result['task_id']
        covered_level = list()
        if "select_id" in result:
            covered_level = result["select_id"]
            covered_level = json.loads(covered_level)
            if level not in covered_level:
                covered_level.append(int(level))
        else:
            covered_level.append(9)
            if level not in covered_level:
                covered_level.append(int(level))
        covered_level = list(set(covered_level))
        covered_level.sort()
        result["select_id"] = json.dumps(covered_level)
        collection.update({"_id": from_otu_table}, {"$set": result}, upsert=False)
        insert_data = list()
        with open(otu_path, 'rb') as v1:
            data = v1.readlines()
            all = 0
            for line1 in data[1:]:
                line1 = line1.rstrip("\r\n")
                line1 = re.split('\t', line1)
                sample_num = line1[1:]
                for i in range(0, len(sample_num)):
                    all += int(float(sample_num[i]))
        with open(otu_path, 'rb') as r:
            head = r.next().strip('\r\n')
            head = re.split('\t', head)
            if head[-1] == "taxonomy":
                new_head = head[1:-1]
            else:
                new_head = head[1:]
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                sample_num = line[1:]
                classify_list = re.split(r"\s*;\s*", line[0])
                otu_detail = dict()
                otu_detail['asv_id'] = from_otu_table
                otu_detail['project_sn'] = project_sn
                otu_detail['task_id'] = task_id
                otu_detail["level_id"] = int(level)
                for cf in classify_list:
                    if cf != "":
                        otu_detail[cf[0:3].lower()] = cf
                count = 0
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = int(float(sample_num[i]))
                    count += int(float(sample_num[i]))
                otu_detail["total"] = count
                otu_detail["percent"] = float(count)/all
                insert_data.append(otu_detail)
        try:
            collection = self.db['asv_detail_level']
            collection.insert_many(insert_data)
        except Exception as e:
            raise Exception("导入asv_detail_level表格失败：{}".format(e))
        else:
            print "导入asv_detail_level表格成功"

    @report_check
    def add_sg_otu_seq_summary(self,otu_path, from_otu_table):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!")

        with open(otu_path, 'rb') as r:
            head = r.next().strip('\r\n')
            head = re.split('\t', head)
            if head[-1] == "taxonomy":
                new_head = head[1:-1]
            else:
                new_head = head[1:]
            sample_num = len(new_head)
            sample_count ={}
            for index in range(sample_num):    #guanqing.zou 20180525
                sample_count[new_head[index]] = 0
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                for index in range(sample_num):    #guanqing.zou 20180525  
                    sample_count[new_head[index]] += int(line[index+1])



        otu_sum = 0
        insert_data = []
        for k in sample_count.keys():
            insert_tmp = {"asv_id":from_otu_table}
            insert_tmp["specimen"] = k
            insert_tmp["read_number"] = sample_count[k]
            insert_data.append(insert_tmp)
            otu_sum += sample_count[k]

        try:
            collection = self.db['asv_specimen_stat']
            collection.insert_many(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入asv_specimen_stat表格失败：{}".format(e))
            self.bind_object.set_error("导入asv_specimen_stat表格失败")
        else:
            self.bind_object.logger.info("导入asv_specimen_stat表格成功")

        # find_task = self.db['asv'].find_one({"_id":from_otu_table})
        # task_id = find_task["task_id"]
        # sg = self.db['valid_sequence_info'].find_one({"task_id":task_id})
        # # add excution of empty result situation; by liulinmeng 20180611
        # if not sg:
        #     amplified_region = "--"
        # else:
        #     amplified_region = sg['amplified_region']
        # try:
        #     collection = self.db['asv_summary']
        #     collection.insert_one({"asv_id":from_otu_table,"samples":sample_num,"sequences":otu_sum,"amplified_region":amplified_region})
        # except Exception as e:
        #     self.bind_object.logger.error("导入asv_summary表格失败：{}".format(e))
        #     self.bind_object.set_error("导入asv_summary表格失败")
        # else:
        #     self.bind_object.logger.info("导入asv_summary表格成功")
    def get_task_info_from_task(self):
        """
        根据task读取origin数据表中的记录
        :return:
        """
        main_collection = self.db['asv']
        name = "ASV_Taxon_Origin"
        result = main_collection.find_one({"task_id": self.main_task_id, "name": name})
        if result:
            settle_params_json = result['settled_params']
            settle_params = json.loads(settle_params_json)
            return settle_params
        else:
            self.bind_object.set_error("未能成功查到asv_origin的任务！")


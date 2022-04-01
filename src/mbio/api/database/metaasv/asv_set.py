# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import re
import json
from bson.objectid import ObjectId
from types import StringTypes


class AsvSet(Base):

    def __init__(self, bind_object):
        super(AsvSet, self).__init__(bind_object)
        self._project_type = 'metaasv'
        self.task_id = ""
        self.name_id = dict()
        self.level_dict = {
            1: "d__",
            2: "k__",
            3: "p__",
            4: "c__",
            5: "o__",
            6: "f__",
            7: "g__",
            8: "s__",
            9: "asv"
        }

    @report_check
    def add_main_table(self, main_id, main_table_id, main_table, species_name=None, pvalue=None, qvalue=None, lda=None, top=None, label=None):
        """
        功能：查询详情表的物种名称，并更新asv集主表字段
        :param main_id: 主表id
        :param main_table_id: 查询主表的id
        :param main_table: 查询主表名称
        :param species_name: 物种名称
        :param pvalue: pvalue值
        :param qvalue: qvalue值
        :param lda: lda值
        :param top: top物种importance值
        :param label: venn标签值
        :return:
        """
        if main_table_id != 0 and not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                self.bind_object.set_error("main_table_id必须为ObjectId对象或其对应的字符串!")
        else:
            main_table_id = main_table_id

        if main_id != 0 and not isinstance(main_id, ObjectId):
            if isinstance(main_id, StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串!")
        else:
            main_id = main_id

        taxon_list = []
        if main_table in ["multiple_group", "two_group", "two_sample"]:
            if main_table in ["multiple_group"]:
                key_relation = "multiple_id"
            else:
                key_relation = "compare_id"
            detail_table = main_table + "_detail"
            if species_name:
                if pvalue:
                    if qvalue:
                        taxon_list = self.add_group_compare(main_table_id, detail_table, key_relation, species_name=species_name, pvalue=pvalue, qvalue=qvalue)
                    else:
                        taxon_list = self.add_group_compare(main_table_id, detail_table, key_relation, species_name=species_name, pvalue=pvalue)
                else:
                    if qvalue:
                        taxon_list = self.add_group_compare(main_table_id, detail_table, key_relation, species_name=species_name, qvalue=qvalue)
                    else:
                        taxon_list = self.add_group_compare(main_table_id, detail_table, key_relation, species_name=species_name)
            else:
                if pvalue:
                    if qvalue:
                        taxon_list = self.add_group_compare(main_table_id, detail_table, key_relation, pvalue=pvalue, qvalue=qvalue)
                    else:
                        taxon_list = self.add_group_compare(main_table_id, detail_table, key_relation, pvalue=pvalue)
                else:
                    if qvalue:
                        taxon_list = self.add_group_compare(main_table_id, detail_table, key_relation, qvalue=qvalue)
                    else:
                        self.bind_object.set_error("没有设置查询条件!")
        if main_table in ["lefse"]:
            detail_table = main_table + "_detail"
            key_relation = "lefse_id"
            if pvalue:
                if lda:
                    taxon_list = self.add_lefse(main_table_id, detail_table,key_relation, pvalue=pvalue, lda=lda)
                else:
                    taxon_list = self.add_lefse(main_table_id, detail_table,key_relation, pvalue=pvalue)
            else:
                if lda:
                    taxon_list = self.add_lefse(main_table_id, detail_table,key_relation, lda=lda)
                else:
                    self.bind_object.set_error("没有设置查询条件!")

        if main_table in ["randomforest"]:
            detail_table = "randomforest_bar"
            key_relation = "randomforest_id"
            taxon_list = self.add_randomforest(main_table_id, detail_table, key_relation, top)
        if main_table in ["venn"]:
            detail_table = "venn_detail"
            key_relation = "venn_id"
            taxon_list = self.add_venn(main_table_id, detail_table, key_relation, species_name, label)

        if len(taxon_list) == 0:
            self.bind_object.set_error("根据筛选条件得到的结果为空，请检查输入参数！")

        collection = self.db[main_table]
        result = collection.find_one({"_id": main_table_id})
        if not result:
            self.bind_object.set_error("在{}表中找不到相应记录".format(main_table))
        if main_table in ["lefse"]:
            params = json.loads(result["params"])
            asv_id = result["asv_id"] ## ObjectId 类型
            start_level = params["start_level"]
            end_level =  params["end_level"]
            asv_list = self.find_asv_from_level(taxon_list, asv_id, start_level)
            asv_list2 = self.find_asv_from_level(taxon_list, asv_id, end_level)
            asv_list = asv_list + asv_list2
        else:
            level_id = result["level_id"] ## int型
            asv_id = result["asv_id"] ## ObjectId 类型
            asv_list = self.find_asv_from_level(taxon_list, asv_id, level_id)
        asv_specimen = self.find_specimen_from_asv(asv_id)
        asv_size = len(asv_list)
        task_id = result["task_id"]
        project_sn = result["project_sn"]
        try:
            main_collection = self.db["asv_set"]
            main_collection.update({"_id": main_id}, {"$set": {"main_id": main_id,
                                                               "task_id": task_id,
                                                               "project_sn": project_sn,
                                                               "asv_list": asv_list,
                                                               "asv_size": asv_size,
                                                               "specimen": asv_specimen}})
        except:
            self.bind_object.set_error("更新主表asv_set主表失败！")
        self.bind_object.logger.info("更新主表asv_set成功！")

    @report_check
    def add_group_compare(self, main_id, detail_table, key_relation, species_name=None, pvalue=None, qvalue=None):
        """
        多组比较、两组比较、两样本比较
        :param main_id: 主表id
        :param detail_table: 详情表名称
        :param key_relation: 主表与详情表关联字段
        :param species_name: 物种名称
        :param pvalue: pvalue值
        :param qvalue: qvalue值
        :return:
        """
        self.bind_object.logger.info("开始查询表")
        collection = self.db[detail_table]
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, StringTypes):
                main_table_id = ObjectId(main_id)
        else:
            main_table_id = main_id
        species_list = []
        try:
            self.bind_object.logger.info("main_table_id: {}\ndetail_table: {}\n key_relation:{}\n".format(main_table_id, detail_table, key_relation))
            if species_name:
                if pvalue:
                    if qvalue:
                        results = collection.find({key_relation: main_table_id, "species_name": {"$regex": species_name}, "pvalue":{"$lt": pvalue}, "qvalue":{"$lt": qvalue}})
                    else:
                        results = collection.find({key_relation: main_table_id, "species_name": {"$regex": species_name}, "pvalue":{"$lt": pvalue}})
                else:
                    if qvalue:
                        results = collection.find({key_relation: main_table_id, "species_name": {"$regex": species_name}, "qvalue":{"$lt": qvalue}})
                    else:
                        results = collection.find({key_relation: main_table_id, "species_name": {"$regex": species_name}})
            else:
                if pvalue:
                    if qvalue:
                        results = collection.find({key_relation: main_table_id, "pvalue":{"$lt": pvalue}, "qvalue":{"$lt": qvalue}})
                    else:
                        results = collection.find({key_relation: main_table_id, "pvalue":{"$lt": pvalue}})
                else:
                    if qvalue:
                        results = collection.find({key_relation: main_table_id,"qvalue":{"$lt": qvalue}})
                    else:
                        self.bind_object.set_error("没有设置查询条件!")
            if results:
                for result in results:
                    species = result["species_name"]
                    if species not in species_list:
                        species_list.append(species)
                return species_list
            else:
                self.bind_object.set_error("在{}表中找不到相应记录".format(detail_table))
        except:
            self.bind_object.set_error("在{}表中找不到相应记录".format(detail_table))

    @report_check
    def add_lefse(self, main_id, detail_table, key_relation, lda=None, pvalue=None):
        """
        查询lefse分析的表
        :param main_id:主表id
        :param detail_table:详情表名称
        :param key_relation:主表与详情表关联字段
        :param lda: lda值
        :param pvalue:pvalue值
        :return:
        """
        self.bind_object.logger.info("开始查询表")
        lefse_collection = self.db[detail_table]
        main_table_id = main_id
        species_list = []
        try:
            if pvalue:
                if lda:
                    results = lefse_collection.find({key_relation: main_table_id, "lda": {"$gt": lda}, "pvalue":{"$lt": pvalue}})
                else:
                    results = lefse_collection.find({key_relation: main_table_id, "pvalue":{"$lt": pvalue}})
            else:
                if lda:
                    results = lefse_collection.find({key_relation: main_table_id, "lda": {"$gt": lda}})
                else:
                    self.bind_object.set_error("没有设置查询条件!")
            if results:
                for result in results:
                    species = result["species_name"]
                    all_species = species.split(".")
                    for sp in all_species:
                        if sp not in species_list:
                            species_list.append(sp)
                return species_list
            else:
                self.bind_object.set_error("在{}表中找不到相应记录".format(detail_table))

        except:
            self.bind_object.set_error("在{}表中找不到相应记录".format(detail_table))

    @report_check
    def add_randomforest(self, main_id, detail_table, key_relation, top):
        """
        随机森林查表
        :param main_id: 主表id
        :param detail_table: 详情表名称
        :param key_relation: 主表与详情表关联字段
        :param top: 丰度前多少的importance
        :return:
        """
        self.bind_object.logger.info("开始查询表")
        random_collection = self.db[detail_table]
        main_table_id = main_id
        accuracy_list = []
        all_species = {}
        species_list = []
        try:
            results = random_collection.find({key_relation: main_table_id})
            if results:
                for result in results:
                    species = result["species_name"]
                    accuracy = result["accuracy"]
                    if species not in all_species:
                        all_species[species] =accuracy
                    if accuracy not in accuracy_list:
                        accuracy_list.append(accuracy)
                accuracy_list = sorted(accuracy_list, reverse=True)
                new_accuracy_list = accuracy_list[0:top]
                for i in new_accuracy_list:
                    species_index = all_species.values().index(i)
                    species_name = all_species.keys()[species_index]
                    if species_name not in species_list:
                        species_list.append(species_name)
                return species_list
            else:
                self.bind_object.set_error("在{}表中找不到相应记录".format(detail_table))
        except:
            self.bind_object.set_error("在{}表中找不到相应记录".format(detail_table))

    @report_check
    def add_venn(self, main_id, detail_table, key_relation, species_name, label):
        """
        查询venn表
        :param main_id:主表id
        :param detail_table:详情表名称
        :param key_relation: 主表与详情表关联字段
        :param species_name: 物种名称
        :param label: label名称
        :return:
        """
        self.bind_object.logger.info("开始查询表")
        venn_collection = self.db[detail_table]
        main_table_id = main_id
        species_list = []
        try:
            results = venn_collection.find_one({key_relation: main_table_id, "specimen" : label})
            if results:
                species_name_display = results["species_name_display"]
                species_name_list = species_name_display.strip().split(",")
                for species in species_name_list:
                    if re.search(r"{}".format(species_name), species):
                        if species not in species_list:
                            species_list.append(species)
                return species_list

            else:
                self.bind_object.set_error("在{}表中找不到相应记录".format(detail_table))
        except:
            self.bind_object.set_error("在{}表中找不到相应记录".format(detail_table))

    def find_asv_from_level(self, taxon_list, asv_id, level_id):
        """
        根据level重asv_detail_level表中查到asv水平的asv
        :param taxon_list: 物种的level水平的名称
        :param asv_id: asv_id
        :param level_id: level水平
        :return:
        """
        self.bind_object.logger.info("开始asv查到对应level水平的asv")
        asv_collection = self.db["asv_detail"]
        main_table_id = asv_id
        level_name = self.level_dict[level_id]
        asv_list = []
        try:
            results = asv_collection.find({"asv_id": main_table_id})
            if results:
                for result in results:
                    taxon_name = result[level_name]
                    asv = result["asv"]
                    if taxon_name in taxon_list:
                        if asv not in asv_list:
                            asv_list.append(asv)
                if len(asv_list) != 0:
                    return asv_list
                else:
                    self.bind_object.set_error("在{}中查到的asv_id（{}）的asv_list为空，请重新筛选条件！".format("asv_detail", str(asv_id)))
            else:
                self.bind_object.set_error("在{}中没有查到对应的asv_id（{}）的详情表".format("asv_detail", str(asv_id)))
        except:
            self.bind_object.set_error("在{}表中找不到相应记录".format("asv_detail"))

    def find_specimen_from_asv(self, asv_id):
        """
        根据asv_id查到specimen
        :param asv_id: asv_id
        :return:
        """
        self.bind_object.logger.info("开始asv查到对应level水平的asv")
        asv_collection = self.db["asv_specimen"]
        main_table_id = asv_id
        specimen_list = []
        try:
            results = asv_collection.find({"asv_id": main_table_id})
            if results:
                for result in results:
                    specimen_id = str(result["specimen_id"])
                    if specimen_id not in specimen_list:
                        specimen_list.append(specimen_id)
                if len(specimen_list) != 0:
                    return ",".join(specimen_list)
                else:
                    self.bind_object.set_error("在{}中未查到对应的样本, 请检查查到的asv_id是否正确！".format("asv_specimen", str(asv_id)))

            else:
                self.bind_object.set_error("在{}中没有查到对应的asv_id（{}）的详情表".format("asv_specimen", str(asv_id)))
        except:
            self.bind_object.set_error("在{}表中找不到相应记录".format("asv_specimen"))

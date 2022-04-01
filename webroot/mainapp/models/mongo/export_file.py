# -*- coding: utf-8 -*-
# __author__ = 'xuting'
# last_modify by qiuping 20170113

import re
import json
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.models.mongo.core.base import Base
# from biocluster.config import Config


class ExportFile(Base):
    def __init__(self, bind_object=None):
        super(ExportFile, self).__init__(bind_object)
        self._project_type = 'meta'
        # self.client = Config().mongo_client
        # self.db = self.client[Config().MONGODB]
        self.level = {
            9: "otu", 8: "s__", 7: "g__", 6: "f__", 5: "o__",
            4: "c__", 3: "p__", 2: "k__", 1: "d__"
        }

    def export_otu_table(self, otuId, targetPath):
        self._try_write(targetPath)
        print "正在导出OTU表格:{}".format(targetPath)
        collection = self.db['sg_otu_specimen']
        my_collection = self.db['sg_specimen']
        results = collection.find({"otu_id": ObjectId(otuId)})
        samples = []
        for result in results:
            id_result = my_collection.find_one({"_id": result["specimen_id"]})
            if not id_result:
                raise Exception("意外错误，样本id:{}在sg_specimen中未找到！")
            samples.append(id_result["specimen_name"])
        collection = self.db['sg_otu_detail']
        with open(targetPath, "wb") as f:
            f.write("OTU ID\t%s\ttaxonomy\n" % "\t".join(samples))
            for col in collection.find({"otu_id": ObjectId(otuId)}):
                line = "%s\t" % col["otu"]
                for s in samples:
                    line += "%s\t" % col[s]
                for cls in ["d__", "k__", "p__", "c__", "o__", "f__", "g__"]:
                    if cls in col.keys():
                        line += "%s; " % col[cls]
                line += col["s__"]
                f.write("%s\n" % line)
        return targetPath

    def export_otu_table_by_level(self, otuId, targetPath, level=9):
        self._project_type = 'meta'
        # self.client = Config().mongo_client
        # self.db = self.client[Config().MONGODB]
        self.level = {
            9: "otu", 8: "s__", 7: "g__", 6: "f__", 5: "o__",
            4: "c__", 3: "p__", 2: "k__", 1: "d__"
        }
        self._try_write(targetPath)
        level = int(level)
        print "正在导出级别为{}:{}的OTU表格:{}".format(level, self.level[level], targetPath)
        collection = self.db['sg_otu_specimen']
        results = collection.find({"otu_id": ObjectId(otuId)})
        if not results.count():
            raise Exception("otu_id: {}在sg_otu_specimen表中未找到！".format(otuId))
        samples = list()
        for result in results:
            if "specimen_id" not in result:
                raise Exception("otu_id:{}错误，请使用新导入的OTU表的id".format(otuId))
            sp_id = result['specimen_id']
            my_collection = self.db['sg_specimen']
            my_result = my_collection.find_one({"_id": sp_id})
            if not my_result:
                raise Exception("意外错误，样本id:{}在sg_specimen表里未找到".format(sp_id))
            samples.append(my_result["specimen_name"])
        collection = self.db['sg_otu_detail']
        name_dic = dict()
        results = collection.find({"otu_id": ObjectId(otuId)})
        if not results.count():
            raise Exception("otu_id: {}在sg_otu_detail表中未找到！".format(otuId))
        for col in results:
            tmp = level + 1
            new_classify_name = self._create_classify_name(col, tmp)
            if new_classify_name not in name_dic:
                name_dic[new_classify_name] = dict()
                for sp in samples:
                    name_dic[new_classify_name][sp] = int(col[sp])
            else:
                for sp in samples:
                    name_dic[new_classify_name][sp] += int(col[sp])
        with open(targetPath, "wb") as f:
            f.write("OTU ID\t%s\n" % "\t".join(samples))
            for k in name_dic.iterkeys():
                line = k
                for s in samples:
                    line += "\t" + str(name_dic[k][s])
                line += "\n"
                f.write(line)
        return targetPath

    def _create_classify_name(self, col, tmp):
        for i in range(1, 10):
            if self.level[i] not in col:
                raise Exception("Mongo数据库中的taxonomy信息不完整, 缺少{}".format(self.level[i]))
        new_col = list()
        for i in range(1, tmp):
            new_col.append(col[self.level[i]])
        return "; ".join(new_col)

    def _try_write(self, targetPath):
        try:
            with open(targetPath, "wb"):
                pass
        except Exception as e:
            raise Exception(e)
        return True

    def _get_objectid(self, data):
        if not isinstance(data, ObjectId):
            if not isinstance(data, StringTypes):
                raise Exception("{}不为ObjectId类型或者其对应的字符串".format(data))
            else:
                try:
                    data = ObjectId(data)
                except:
                    raise Exception("{}不为ObjectId类型或者其对应的字符串".format(data))
        return data

    def export_group_table(self, groupId, category_name, targetPath):
        self._try_write(targetPath)
        print "正在按分组id和组名导出group表:{}".format(targetPath)
        if groupId == "all":
            with open(targetPath, "wb") as f:
                f.write("#sample\t" + "##empty_group##" + "\n")
            return targetPath
        group_table = self.db['sg_specimen_group']
        group_name_list = list()
        group_name = category_name
        group_name_list = re.split(',', group_name)
        groupId = self._get_objectid(groupId)
        group_schema = group_table.find_one({"_id": ObjectId(groupId)})
        if not group_schema:
            raise Exception("无法根据传入的group_id:{}在sg_specimen_group表里找到相应的记录".format(groupId))
        c_name = group_schema["category_names"]

        with open(targetPath, "wb") as f:
            f.write("#sample\t" + group_schema["group_name"] + "\n")

        sample_table_name = 'sg_specimen'
        sample_table = self.db[sample_table_name]
        index_list = self._get_index_list(group_name_list, c_name)
        sample_id = list()  # [[id名,组名], [id名, 组名]，...]
        sample_name = list()  # [[样本名, 组名], [样本名, 组名], ...]
        for i in index_list:
            for k in group_schema['specimen_names'][i]:
                # k 样品ID
                sample_id.append([k, c_name[i]])
        for pair in sample_id:
            result = sample_table.find_one({"_id": ObjectId(pair[0])})
            if not result:
                raise Exception("无法根据传入的group_id:{}在样本表{}里找到相应的样本名".format(groupId, sample_table_name))
            sample_name.append([result['specimen_name'], pair[1]])

        with open(targetPath, "ab") as f:
            for pair in sample_name:
                f.write("{}\t{}\n".format(pair[0], pair[1]))
        return targetPath

    def _get_index_list(self, group_name_list, c_name):
        length = len(c_name)
        index_list = list()
        for i in range(length):
            for g_name in group_name_list:
                if g_name == c_name[i]:
                    index_list.append(i)
                    break
        return index_list

    def export_group_table_by_detail(self, groupId, groupDetail, targetPath):
        """
        按分组的详细信息获取group表
        """
        self._try_write(targetPath)
        print "正在按分组详细信息导出group表，目标路径:{}".format(targetPath)
        groupId = self._get_objectid(groupId)
        group_table = self.db['sg_specimen_group']
        if not isinstance(groupDetail, dict):
            try:
                table_dict = json.loads(groupDetail)
            except Exception:
                raise Exception("生成group表失败，传入的groupDetail不是一个字典或者是字典对应的字符串")
        if not isinstance(table_dict, dict):
            raise Exception("生成group表失败，传入的groupDetail不是一个字典或者是字典对应的字符串")
        group_schema = group_table.find_one({"_id": ObjectId(groupId)})
        if not group_schema:
            raise Exception("无法根据传入的group_id:{}在sg_specimen_group表里找到相应的记录".format(groupId))
        with open(targetPath, "wb") as f:
            f.write("#sample\t" + group_schema["group_name"] + "\n")

        sample_table_name = 'sg_specimen'
        sample_table = self.db[sample_table_name]
        with open(targetPath, "ab") as f:
            for k in table_dict:
                for sp_id in table_dict[k]:
                    sp = sample_table.find_one({"_id": ObjectId(sp_id)})
                    if not sp:
                        raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                    else:
                        sp_name = sp["specimen_name"]
                    f.write("{}\t{}\n".format(sp_name, k))
        return targetPath

    def export_cascading_table_by_detail(self, groupId, groupDetail, targetPath):
        """
        根据group_detail生成group表或者二级group表
        """
        self._try_write(targetPath)
        print "正在按分组详细信息导出层级group表，目标路径:{}".format(targetPath)
        groupId = self._get_objectid(groupId)
        group_table = self.db['sg_specimen_group']
        if not isinstance(groupDetail, dict):
            try:
                table_list = json.loads(groupDetail)
            except Exception:
                raise Exception("生成group表失败，传入的groupDetail不是一个数组或者是数组对应的字符串")
        if not isinstance(table_list, dict):
            raise Exception("生成group表失败，传入的groupDetail不是一个数组或者是数组对应的字符串")
        group_schema = group_table.find_one({"_id": ObjectId(groupId)})
        if not group_schema:
            raise Exception("无法根据传入的group_id:{}在sg_specimen_group表里找到相应的记录".format(groupId))
        with open(targetPath, "wb") as f:
            f.write("#sample\t" + group_schema["group_name"])

        sample_table_name = 'sg_specimen'
        sample_table = self.db[sample_table_name]
        self._write_cascading_table(table_list, targetPath, sample_table, sample_table_name)
        return targetPath

    def _write_class1_table(self, table_list, targetPath, sample_table, sample_table_name):
        table_dict = table_list[0]
        with open(targetPath, "ab") as f:
            f.write("\n")
            for k in table_dict:
                for sp_id in table_dict[k]:
                    sp = sample_table.find_one({"_id": ObjectId(sp_id)})
                    if not sp:
                        raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                    else:
                        sp_name = sp["specimen_name"]
                    f.write("{}\t{}\n".format(sp_name, k))

    def _write_class2_table(self, table_list, targetPath, sample_table, sample_table_name):
        dict1 = dict()  # 样本id -> 分组名
        dict2 = dict()
        for k in table_list[0]:
            for id_ in table_list[0][k]:
                if id_ in dict1:
                    raise Exception("样本id:{}在分组{}，{}中都出现，一个样本在同一级别上只能属于一个分组！".format(id_, dict1[id_], k))
                dict1[id_] = k
        for k in table_list[1]:
            for id_ in table_list[1][k]:
                if id_ not in dict1:
                    raise Exception("样本id:{}在第一级的分组中未出现".format(id_))
                if id_ in dict2:
                    raise Exception("样本id:{}在分组{}，{}中都出现，一个样本在同一级别上只能属于一个分组！".format(id_, dict1[id_], k))
                dict2[id_] = k
        if len(dict1) != len(dict2):
            raise Exception("一级分组中的某些样本id在二级分组中未找到！")
        with open(targetPath, "ab") as f:
            f.write("\tsecond_group\n")
            for k in dict1:
                sp = sample_table.find_one({"_id": ObjectId(k)})
                if not sp:
                    raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(k, sample_table_name))
                else:
                    sp_name = sp["specimen_name"]
                f.write("{}\t{}\t{}\n".format(sp_name, dict1[k], dict2[k]))

    def _write_cascading_table(self, table_list, targetPath, sample_table, sample_table_name):
        length = len(table_list)
        if length == 1:
            self._write_class1_table(table_list, targetPath, sample_table, sample_table_name)
        elif length == 2:
            self._write_class2_table(table_list, targetPath, sample_table, sample_table_name)
        else:
            raise Exception("group_detal字段含有三个或以上的字典")

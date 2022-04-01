# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import re
import pandas as pd
import datetime
import json
from bson.objectid import ObjectId
from types import StringTypes


class Venn(Base):
    def __init__(self, bind_object):
        super(Venn, self).__init__(bind_object)
        self._project_type = 'metaasv'
        self.new_otu_id = list()
        self.single = dict()
        self.num = dict()

    @report_check
    def create_venn_table(self, params, group_id, level_id, from_otu_table=0, name=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!")
        if not isinstance(group_id, ObjectId):
            if isinstance(group_id, StringTypes):
                group_id = ObjectId(group_id)
            else:
                self.bind_object.set_error("group_id必须为ObjectId对象或其对应的字符串!")
        collection = self.db["asv"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("asv表里找不到相应的记录")
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "正在计算venn表格..."
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "level_id": level_id,
            "otu_id": from_otu_table,
            "group_id": group_id,
            "status": "end",
            "desc": desc,
            "name": self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else "venn表格",
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["venn"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_venn_detail(self, venn_table, venn_id):
        """
        需要计算每个标签所含有的物种数，
        :param venn_table:
        :param venn_id:
        :return:
        """
        data_list = []
        if not isinstance(venn_id, ObjectId):
            if isinstance(venn_id, StringTypes):
                venn_id = ObjectId(venn_id)
            else:
                self.bind_object.set_error("venn_id必须为ObjectId对象或其对应的字符串!")
        with open(venn_table, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                label = line[0].split(" only")[0].strip()
                number = line[1]
                try:
                    full_name = line[2]
                    tmp_name = re.split(',', line[2])
                except:
                    full_name = ""
                    tmp_name = []
                name_list = []
                for cla_info in tmp_name:
                    cla_info = re.split('; ', cla_info)
                    my_name = cla_info[-1]
                    name_list.append(my_name)
                display_name = ",".join(name_list)
                insert_data = {
                    'venn_id': venn_id,
                    'specimen': label,
                    'display_count': number,
                    'species_name': full_name,
                    'species_name_display':display_name,
                    'speceies_count': number
                }
                data_list.append(insert_data)
        try:
            collection = self.db["venn_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入Venn数据出错:%s" % e)
            self.bind_object.set_error("导入Venn数据出错")
        else:
            self.bind_object.logger.error("导入Venn数据成功")
        try:
            main_collection = self.db["venn"]
            settled_params = {"software" : "R-3.3.1 (stat)"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            venn_data = {
                "venn_data": {"name":"species_name",
                            "category": "specimen"}}
            venn_data_json = json.dumps(venn_data, sort_keys=True, separators=(',', ':'))
            table_data = {
                "table_data": {"specimen":"Group_label",
                            "display_count": "Species_num"}}
            table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
            pie_data = {
                "pie_data": {"name":"specimen",
                            "category": "species_name_display"}}
            pie_data_json = json.dumps(pie_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": venn_id},{"$set": {"main_id": venn_id,
                                                              "settled_params": settled_params_json,
                                                              "table_data": table_data_json,
                                                              "pie_data": pie_data_json,
                                                              "venn_data":venn_data_json}})
        except Exception, e:
            self.bind_object.logger.error("更新Venn数据出错:%s" % e)
            self.bind_object.set_error("更新Venn数据出错")
        else:
            self.bind_object.logger.error("更新Venn数据成功")

    @report_check
    def add_venn_pie(self, venn_table, venn_id, asv_table=None, group_table=None):
        """
        需要计算每个标签所含有的物种数，
        :param venn_table:
        :param venn_id:
        :param asv_table:
        :param group_table:
        :return:
        """
        data_list = []
        if not isinstance(venn_id, ObjectId):
            if isinstance(venn_id, StringTypes):
                venn_id = ObjectId(venn_id)
            else:
                self.bind_object.set_error("venn_id必须为ObjectId对象或其对应的字符串!")
        data = pd.read_table(asv_table, sep="\t", header=0)##分组名称
        columns = list(data.columns)
        data["level"] = (data["OTU ID"].str.split(";").str)[-1].str.strip()
        new_columns = ["level"] + columns[1:]
        all_data = data[new_columns]
        all_data = all_data.set_index("level")
        with open(venn_table, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                label = line[0].split(" only")[0].strip()
                label_list = label.split("&")
                label_list = [x.strip() for x in label_list]
                label_list = [x.split(" only")[0].strip() for x in label_list]
                try:
                    tmp_name = re.split(',', line[2])
                except:
                    tmp_name = ""
                name_list = []
                for cla_info in tmp_name:
                    cla_info = re.split('; ', cla_info)
                    my_name = cla_info[-1].strip()
                    if my_name not in name_list:
                        name_list.append(my_name)
                # all_species_number = 0
                # for sp_name in name_list:
                #     spe_number1 = 0
                #     for group in label_list:
                #         sp_num = all_data[group].loc[sp_name]
                #         spe_number1 += sp_num
                #     all_species_number += spe_number1
                if tmp_name != "":
                    for sp_name in name_list:
                        spe_number = 0
                        for group in label_list:
                            # self.bind_object.logger.info("group: {}".format(group))
                            # self.bind_object.logger.info("sp_name: {}".format(sp_name))
                            try:
                                sp_num = all_data[group].loc[sp_name] ###可能会出现某个水平上的值是有相同名称的（这种是bug，taxon可能是错误的，没有加上级水平的名称）
                                spe_number += int(sp_num)
                            except:
                                sp_num = all_data[group].loc[sp_name] ###可能会出现某个水平上的值是有相同名称的（这种是bug，taxon可能是错误的，没有加上级水平的名称）
                                spe_number += int(max(list(sp_num)))
                        # relative_sp_num = float(spe_number) / all_species_number
                        if int(spe_number) != 0:
                            insert_data = {
                                'venn_id': venn_id,
                                'specimen': label,
                                'species_name': sp_name,
                                'species_number':spe_number
                            }
                            data_list.append(insert_data)
                else:
                    insert_data = {
                            'venn_id': venn_id,
                            'specimen': label,
                            'species_name': "",
                            'species_number':0
                        }
                    data_list.append(insert_data)
        try:
            collection = self.db["venn_pie"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入Venn数据出错:%s" % e)
            self.bind_object.set_error("导入Venn数据出错")
        else:
            self.bind_object.logger.error("导入Venn数据成功")

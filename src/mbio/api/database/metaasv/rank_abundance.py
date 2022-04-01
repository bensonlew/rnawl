# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import re
import datetime
from random import choice
import json
import string
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.libs.param_pack import group_detail_sort, param_pack
from mbio.packages.metaasv.common_function import find_group_name
from bson.son import SON


class RankAbundance(Base):
    def __init__(self, bind_object):
        super(RankAbundance, self).__init__(bind_object)
        self._project_type = 'metaasv'
        self.main_task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])

    @report_check
    def add_rank(self, asv_id=None, params=None, name=None, spname_spid=None, group_id=None):
        if asv_id and not isinstance(asv_id, ObjectId):
            if isinstance(asv_id, StringTypes):
                origin_asv_id = ObjectId(asv_id)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!")
        else:
            origin_asv_id = asv_id
        project_sn = self.bind_object.sheet.project_sn
        if name:
            record_name = name
        else:
            record_name = self.bind_object.sheet.main_table_name
        task_id = self.main_task_id
        if spname_spid and params:
            if group_id not in [None, "All", "all", "ALL"]:
                ## 调用common模块，功能将导入的分组方案返回group_detail
                group_detail = find_group_name(task_id)
            else:
                group_detail = {'All': [str(i) for i in spname_spid.values()]}
            params['group_detail'] = group_detail_sort(group_detail)
        params = param_pack(params)
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "asv_id": origin_asv_id,
            "level_id": 7,
            "status": "end",
            "desc": "Rank_abundance主表",
            "submit_location": "rankabundace",
            "name": record_name if name else "Rank_abundance_Origin",
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["rank_abundance"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_rank_detail(self, file_path, rank_id):
        """
        导入详情表
        :param file_path: 文件路径
        :param rank_id: 主表id string类型
        :return:
        """
        if not isinstance(rank_id, ObjectId):
            if isinstance(rank_id, StringTypes):
                rank_id = ObjectId(rank_id)
            else:
                self.bind_object.set_error("rank_id必须为ObjectId对象或其对应的字符串!")
        all_data = {}
        with open(file_path, 'rb') as r:
            header = r.next().rstrip("\n")
            header = re.split('\t', header)
            asv_list = header[1:] ##导入MongoDB的key
            sp_list = []
            for line in r:
                line = line.strip().split("\t")
                sample_name = line[0]
                if sample_name not in sp_list:
                    sp_list.append(sample_name)
                all_data[sample_name] = line[1:]
        category_list = []
        rank_abundance_detail = []
        for sample in sp_list:
            for asv in asv_list:
                x = int(asv)
                y = float(all_data[sample][x-1])
                if x not in category_list:
                    category_list.append(x)
                insert_data = {
                    "rank_id": rank_id,
                    "specimen": sample,
                    "name" : "",
                }
                if y != 0.0:
                    insert_data["x"] = x
                    insert_data["y"] = y
                    all_insert_data = SON(insert_data)
                    rank_abundance_detail.append(all_insert_data)
        try:
            collection_detail = self.db['rank_abundance_detail']
            collection_detail.insert_many(rank_abundance_detail)
            self.bind_object.logger.info("导入rank_abundance_detail表格成功！")
            main_collection = self.db["rank_abundance"]
            settled_params = {"software" : "python-2.7"}
            line_data = {
                "line_data": {"name":"specimen",
                        "condition": {"type":["pan", "core"]}
                }
            }
            line_data_json = json.dumps(line_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id":rank_id},{"$set":{"main_id": rank_id,
                                                            "line_data": line_data_json,
                                                            "category_x": category_list,
                                                            "settled_params": settled_params}})
            self.bind_object.logger.info("更新rank_abundance表格成功！")
        except Exception as e:
            self.bind_object.logger.error("导入rank_abundance_detail表格{}信息出错:{}".format(file_path, e))
            self.bind_object.set_error("导入rank_abundance_detail表出错")
        else:
            self.bind_object.logger.info("导入rank_abundance_detail表格{}成功".format(file_path))


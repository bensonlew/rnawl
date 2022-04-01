# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify:20180927

import os
import json
import datetime
from api_base import ApiBase


class PopHistory(ApiBase):
    """
    种群历史导表
    """
    def __init__(self, bind_object):
        super(PopHistory, self).__init__(bind_object)
        self._project_type = "dna_evolution"
        self.project_sn = self.bind_object.sheet.project_sn
        self.task_id = self.bind_object.sheet.id.split("_PopHistory")[0]

    def add_sg_history(self, project_sn, task_id, params=None, name=None):
        """
        主表sg_history
        """
        if params:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        else:
            params = None
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "params": params if params else "",
            "name": name if name else "origin_pop_history",
            "desc": "种群历史",
        }
        main_id = self.db['sg_history'].insert_one(insert_data).inserted_id
        self.update_db_record("sg_history", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_history_detail(self, history_id, psmc_stat):
        """
        sg_history_detail
        psmc_stat: all_psmc_stat.xls
        """
        history_id = self.check_objectid(history_id)
        self.check_exists(psmc_stat)
        data_list = []
        with open(psmc_stat, "r") as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "history_id": history_id,
                    "generation": item[0]
                }
                for i in range(1, len(header)):
                    insert_data[header[i]] = item[i]
                data_list.append(insert_data)
        if "all" in header[1:]:
            header_ = header[1:]
            header_.remove("all")
        else:
            header_ = header[1:]
        self.update_db_record("sg_history", {"main_id": history_id}, {"header": header_})
        if len(data_list) == 0:
            self.bind_object.logger.info("{}的结果为空！".format(psmc_stat))
        else:
            self.col_insert_data("sg_history_detail", data_list)

    def add_sg_curve(self, history_id):
        """
        主表sg_curve
        """
        history_id = self.check_objectid(history_id)
        insert_data = {
            "task_id": self.task_id,
            "origin_id": history_id,
            "name": "",
            "type": 1,
            "location": "history_curve",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        main_id = self.db['sg_curve'].insert_one(insert_data).inserted_id
        return main_id

    def add_sg_curve_detail(self, curve_id, psmc_stat):
        """
        sg_curve_detail
        psmc_stat: all_psmc_stat.xls
        """
        curve_id = self.check_objectid(curve_id)
        self.check_exists(psmc_stat)
        data_list, categories = [], []
        data = {}
        with open(psmc_stat, "r") as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            for i in range(1, len(header)):
                data[header[i]] = []
            for line in lines[1:]:
                item = line.strip().split("\t")
                categories.append(item[0])
                for i in range(1, len(header)):
                    data[header[i]].append(round(float(item[i]), 5))
        for key in data.keys():
            insert_data = {
                "curve_id": curve_id,
                "name": key,
                "value": data[key]
            }
            data_list.append(insert_data)
        self.update_db_record("sg_curve", {"_id": curve_id}, {"categories": categories})
        if len(data_list) == 0:
            self.bind_object.logger.info("{}的结果为空！".format(psmc_stat))
        else:
            self.col_insert_data("sg_curve_detail", data_list)


if __name__ == "__main__":
    a = PopHistory(None)
    project_sn = "test_zj"
    task_id = "tsg_32120"
    # psmc_stat = "/mnt/ilustre/users/sanger-dev/workspace/20181009/Single_pop_history2/PopHistory/output/all_psmc_stat.xls"
    # history_id = a.add_sg_history(project_sn, task_id, params=None, name=None)
    # a.add_sg_history_detail(history_id, psmc_stat)
    # curve_id = a.add_sg_curve(task_id, history_id)
    # a.add_sg_curve_detail(curve_id, psmc_stat)
    curve_id = "5c072fd1a4e1af1f7ca4c573"
    psmc_stat = "/mnt/ilustre/users/sanger-dev/workspace/20181204/Evolution_tsg_32990/PopHistory/PsmcStat/output/all_psmc_stat.xls"
    a.add_sg_curve_detail(curve_id, psmc_stat)

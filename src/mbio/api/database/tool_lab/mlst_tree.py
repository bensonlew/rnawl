# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
#  20210409

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
from Bio import Phylo
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class MlstTree(ApiBase):
    def __init__(self, bind_object):
        super(MlstTree, self).__init__(bind_object)
        self._project_type = "tool_lab"

    def add_mlsttree(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "MlstTree" + "_"
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='MlstTree',
                params=params,
                status="start")
            main_id = self.create_db_table('mlst_tree', [main_info])
        else:
            main_id = ObjectId(main_id)
            try:
                self.update_db_record('mlst_tree', main_id, )
            except Exception as e:
                self.bind_object.logger.error("导入mlst_tree数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入mlst_tree数据成功")
        return main_id

    def add_mlsttree_detail(self, main_id, file1, file2=None):
        with open(file1, "r") as f:
            lines = f.readlines()
            list1 = [{"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample Name",
                      "type": "string"},
                     {"field": "st_id", "filter": "false", "sort": "false", "title": "ST",
                      "type": "string"}, ]
            funs = lines[0].strip().split("\t")
            for i in funs[2:]:
                list1.append({"field": i, "filter": "false", "sort": "false", "title": i, "type": "string"})
            table1_dict = {"column": list1, "condition": {}}
            data_list = []
            for line in lines[1:]:
                lin = line.strip().split("\t")
                insert_data = {
                    "mlst_id": ObjectId(main_id),
                    "sample_name": lin[0],
                    "st_id": lin[1],
                }
                for i in range(2, len(lin)):
                    insert_data[funs[i]] = lin[i]
                data_list.append(insert_data)
            self.create_db_table('mlst_tree_st', data_list)
        table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
        if file2:
            tree = Phylo.read(file2, "newick")
            list2 = []
            for i in tree.get_terminals():
                list2.append(str(i))
            with open(file2, "r") as f:
                line = f.readline()
            tree2 = line.strip()
            insert_data1 = {
                "mlst_id": ObjectId(main_id),
                "direction": "circle",
                "group": [{"groupname": "all", "value": list2}],
                "type": "tree",
                "data": tree2,
                "name": "",
            }
            self.create_db_table('mlst_tree_value', [insert_data1])
            table1_dict2 = {"name": "name", "condition": {"type": "tree"}}
            table1_info1 = json.dumps(table1_dict2, sort_keys=False, separators=(',', ':'))
            try:
                self.bind_object.logger.info(tree)
                self.update_db_record('mlst_tree', main_id, main_id=main_id, stat_table=table1_info,
                                      tree_data=table1_info1, status="end")
            except Exception as e:
                self.bind_object.logger.error("导入mlst_tree_st数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入mlst_tree_st数据成功")
        else:
            table1_dict2 = {"name": "name", "condition": {"type": "tree"}}
            table1_info1 = json.dumps(table1_dict2, sort_keys=False, separators=(',', ':'))
            try:
                self.update_db_record('mlst_tree', main_id, main_id=main_id, stat_table=table1_info,
                                      tree_data=table1_info1,status="end")
            except Exception as e:
                self.bind_object.logger.error("导入mlst_tree_st数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入mlst_tree_st数据成功")


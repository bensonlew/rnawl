# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last_modify: liubinxu 20181121
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import pandas as pd
import datetime
import json
import unittest
import os

from api_base import ApiBase


class Minet(ApiBase):
    def __init__(self, bind_object):
        super(Minet, self).__init__(bind_object)

    def import_net(self, result_dir=None, params=None, tf=None, target=None, nodes=None, edges=None):
        task_id = self.bind_object.sheet.id
        self.bind_object.logger.info("开始导入mirna调控网络")
        self.remove_table_by_main_record(main_table='minet', task_id=task_id, detail_table=['minet_tf', 'minet_target', 'minet_nodes', 'minet_edges'], detail_table_key='minet_id')
        minet_id = self.add_minet_main(result_dir, params)
        self.add_minet_tf(minet_id, tf)
        self.add_minet_target(minet_id, target)
        self.add_minet_nodes(minet_id, nodes)
        self.add_minet_edges(minet_id, edges)
        self.db['minet'].update({"_id": minet_id},
                                   {"$set": {"status": "end", "main_id": minet_id}})

    def import_net_webroot(self, main_id=None, result_dir=None, params=None, tf=None, target=None, nodes=None, edges=None):
        task_id = self.bind_object.sheet.id
        self.bind_object.logger.info("开始导入mirna调控网络")
        # self.remove_table_by_main_record(main_table='minet', task_id=task_id, detail_table=['minet_tf', 'minet_target', 'minet_nodes', 'minet_edges'], detail_table_key='minet_id')
        minet_id = main_id
        if os.path.exists(tf):
            self.add_minet_tf(minet_id, tf)
        self.add_minet_target(minet_id, target)
        self.add_minet_nodes(minet_id, nodes)
        self.add_minet_edges(minet_id, edges)
        self.db['minet'].update({"_id": minet_id},
                                   {"$set": {"status": "end", "main_id": minet_id}})

    @report_check
    def add_minet_main(self, result_dir, params):
        params_dict = params
        params_dict["submit_location"] = "mirna_mrna_net"
        params_dict["task_type"] = "workflow"
        params_dict["task_id"] = self.bind_object.sheet.id
        data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "status": "start",
            "name": 'miRNA_mRNA_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3],
            "desc": "",
            "created_ts": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            "params": json.dumps(params_dict, sort_keys=True, separators=(',', ':')),
            "result_dir": result_dir
        }
        col = self.db["minet"]
        main_id = col.insert_one(data).inserted_id
        self.bind_object.logger.info("主表创建成功")
        return main_id

    @report_check
    def add_minet_tf(self, minet_id, tf):
        if not isinstance(minet_id, ObjectId):
            if isinstance(minet_id, StringTypes):
                minet_id = ObjectId(minet_id)
            else:
                raise Exception("minet_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        with open(tf, "rb") as tf_f:
            head_line = tf_f.readline()
            header = head_line.strip("\n").split("\t")
            for line in tf_f:
                cols = line.strip("\n").split("\t")
                try:
                    cols[header.index("pvalue")] = float( cols[header.index("pvalue")] )
                except:
                    cols[header.index("pvalue")] = ""
                data = zip(header, cols)
                data.append(('minet_id', minet_id))
                data = SON(data)
                data_list.append(data)

        try:
            self.create_db_table('minet_tf', data_list)
        except Exception as e:
            self.bind_object.set_error("导入TF信息出错!")
        else:
            self.bind_object.logger.info("导入TF信息成功!")


    @report_check
    def add_minet_target(self, minet_id, target):
        if not isinstance(minet_id, ObjectId):
            if isinstance(minet_id, StringTypes):
                minet_id = ObjectId(minet_id)
            else:
                raise Exception("minet_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        with open(target, "rb") as target_f:
            head_line = target_f.readline()
            header = head_line.strip("\n").split("\t")
            for line in target_f:
                cols = line.strip("\n").split("\t")
                if len(cols) == 9:
                    cols[6] = float(cols[6])
                    cols[7] = float(cols[7])
                    cols[8] = float(cols[8])
                data = zip(header, cols)
                data.append(('minet_id', minet_id))
                data.append(('evidence1', "yes"))
                data = SON(data)
                data_list.append(data)

        try:
            self.create_db_table('minet_target', data_list)
        except Exception as e:
            self.bind_object.set_error("导入靶基因相关系性出错!")
        else:
            self.bind_object.logger.info("导入靶基因相关性信息成功!")


    @report_check
    def add_minet_nodes(self, minet_id, nodes):
        if not isinstance(minet_id, ObjectId):
            if isinstance(minet_id, StringTypes):
                minet_id = ObjectId(minet_id)
            else:
                raise Exception("minet_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        with open(nodes, "rb") as nodes_f:
            head_line = nodes_f.readline()
            header = head_line.strip("\n").split("\t")
            for line in nodes_f:
                cols = line.strip("\n").split("\t")
                data = zip(header, cols)
                data.append(('minet_id', minet_id))
                data = SON(data)
                data_list.append(data)

        try:
            self.create_db_table('minet_nodes', data_list)
        except Exception as e:
            self.bind_object.set_error("导入nodes信息出错!")
        else:
            self.bind_object.logger.info("导入nodes信息成功!")

    @report_check
    def add_minet_edges(self, minet_id, edges):
        if not isinstance(minet_id, ObjectId):
            if isinstance(minet_id, StringTypes):
                minet_id = ObjectId(minet_id)
            else:
                raise Exception("minet_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        with open(edges, "rb") as edges_f:
            head_line = edges_f.readline()
            header = head_line.strip("\n").split("\t")
            for line in edges_f:
                cols = line.strip("\n").split("\t")
                cols[7] = float(cols[7])
                cols[10] = abs(float(cols[10]))
                cols[11] = abs(float(cols[11]))
                cols[12] = abs(float(cols[12]))
                if cols[10] == 0:
                    cols[10] = 0.3
                data = zip(header, cols)
                data.append(('minet_id', minet_id))
                data.append(('evid1', "yes"))
                evid = "no"
                if cols[8] != "":
                    evid = "yes"
                data.append(('evid2', evid))
                data = SON(data)
                data_list.append(data)

        try:
            self.create_db_table('minet_edges', data_list)
        except Exception as e:
            self.bind_object.set_error("导入edges信息出错!")
        else:
            self.bind_object.logger.info("导入edges信息成功!")


if __name__ == '__main__':


    class TestFunction(unittest.TestCase):
        """
        测试导表函数
        """

        def test_mongo(test):
            from mbio.workflows.small_rna.small_rna_test_api import SmallRnaTestApiWorkflow
            from biocluster.wsheet import Sheet
            import random
            data = {
                "id": "small_rna",
                #+ str(random.randint(1,10000)),
                #"id": "denovo_rna_v2",
                "project_sn": "small_rna",
                #+ str(random.randint(1,10000)),
                "type": "workflow",
                "name": "small_rna.small_rna_test_api",
                "options": {

                }
            }
            wsheet = Sheet(data=data)
            wf = SmallRnaTestApiWorkflow(wsheet)

            wf.IMPORT_REPORT_DATA = True
            wf.IMPORT_REPORT_AFTER_END = False
            wf.test_api = wf.api.api("small_rna.minet")

            params = {
                "mirna_set": "5bc05fe4a4e1af0f2b2e0e5a",
                "tf_db": 0,
                "corr_method": "pearson",
                "corr_threshold": 0.5,
                "exp_matrix": "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA/exp",
                "group_dict": {"DMFJ":["DMFJ_2","DMFJ_3"],"FHZHQ":["FHZHQ_2","FHZHQ_3"],"LGJB":["LGJB_2","LGJB_3"]},
                "samples":["LGJB_2", "LGJB_3", "FHZHQ_2", "FHZHQ_3", "DMFJ_2", "DMFJ_3"]
            }

            test_dir = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA/"
            tf = test_dir + 'binding.xls'
            corr = test_dir + 'corr.xls'
            nodes = test_dir + 'node.xls'
            edges = test_dir + 'edge.xls'
            wf.test_api.import_net(test_dir, params, tf=tf, target=corr, nodes=nodes, edges=edges)
    unittest.main()

# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last_modify: liubinxu 20190305
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
from mbio.api.database.lnc_rna.api_base import ApiBase


class Cerna(ApiBase):
    def __init__(self, bind_object):
        super(Cerna, self).__init__(bind_object)

    def import_net(self, result_dir=None, params=None, mi2mrna=None, mi2lncrna=None, ce_corr=None, nodes=None, edges=None):
        task_id = self.bind_object.sheet.id
        self.bind_object.logger.info("开始导入cerna调控网络")
        self.remove_table_by_main_record(main_table='sg_cerna', task_id=task_id, detail_table=['sg_cerna_mi2mrna', 'sg_cerna_mi2lncrna', 'sg_cerna_cerna', 'sg_cerna_nodes', 'sg_cerna_edges'], detail_table_key='cerna_id')
        cerna_id = self.add_cerna_main(result_dir, params)
        self.add_cerna_mi2ce(cerna_id, mi2mrna, mi2lncrna)
        self.add_cerna(cerna_id, ce_corr)
        self.add_cerna_nodes(cerna_id, nodes)
        self.add_cerna_edges(cerna_id, edges)
        self.db['sg_cerna'].update({"_id": cerna_id},
                                   {"$set": {"status": "end", "main_id": cerna_id}})

    def import_net_webroot(self, main_id=None, result_dir=None, params=None, mi2mrna=None, mi2lncrna=None, ce_corr=None, nodes=None, edges=None):
        task_id = self.bind_object.sheet.id
        self.bind_object.logger.info("开始导入cerna调控网络")
        cerna_id = ObjectId(main_id)
        self.add_cerna_mi2ce(cerna_id, mi2mrna, mi2lncrna)
        self.add_cerna(cerna_id, ce_corr)
        self.add_cerna_nodes(cerna_id, nodes)
        self.add_cerna_edges(cerna_id, edges)
        self.db['sg_cerna'].update({"_id": cerna_id},
                                   {"$set": {"status": "end", "main_id": cerna_id}})

    @report_check
    def add_cerna_main(self, result_dir, params):
        params_dict = params
        params_dict["submit_location"] = "cerna"
        params_dict["task_type"] = "workflow"
        params_dict["task_id"] = self.bind_object.sheet.id
        data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "status": "start",
            "name": 'ceRNA' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3],
            "desc": "",
            "created_ts": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            "params": json.dumps(params_dict, sort_keys=True, separators=(',', ':')),
            "result_dir": result_dir
        }
        col = self.db["sg_cerna"]
        main_id = col.insert_one(data).inserted_id
        self.bind_object.logger.info("主表创建成功")
        return main_id

    @report_check
    def add_cerna_mi2ce(self, cerna_id, mi2mrna, mi2lncrna):
        '''
        导入miRNA与靶基因相关性信息
        '''
        if not isinstance(cerna_id, ObjectId):
            if isinstance(cerna_id, StringTypes):
                cerna_id = ObjectId(cerna_id)
            else:
                raise Exception("cerna_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []

        with open(mi2mrna, "rb") as tf_f:
            head_line = tf_f.readline()
            header = head_line.strip("\n").split("\t")
            if "energy" in header:
                energy_col = header.index('energy')
            for line in tf_f:
                cols = line.strip("\n").split("\t")
                cols[6] = float(cols[6])
                if (cols[7] == "None"):
                    cols[7] = 0
                else:
                    cols[7] = float(cols[7])

                if "energy" in header:
                    try:
                        cols[energy_col] = float(cols[energy_col])
                    except:
                        cols[energy_col] = 0
                if len(cols) == 11:
                    cols[8] = float(cols[8])
                    cols[9] = float(cols[9])
                    cols[10] = float(cols[10])
                data = zip(header, cols)
                data.append(('cerna_id', cerna_id))
                data = SON(data)
                data_list.append(data)

        try:
            self.create_db_table('sg_cerna_mi2mrna', data_list)
        except Exception as e:
            self.bind_object.set_error("导入miRNA mRNA 对应关系信息出错!")
        else:
            self.bind_object.logger.info("导入miRNA mRNA 对应关系信息成功!")

        data_list = []

        with open(mi2lncrna, "rb") as tf_f:
            head_line = tf_f.readline()
            header = head_line.strip("\n").split("\t")
            for line in tf_f:
                cols = line.strip("\n").split("\t")
                cols[6] = float(cols[6])
                if (cols[7] == "None"):
                    cols[7] = 0
                else:
                    cols[7] = float(cols[7])
                if len(cols) == 11:
                    cols[8] = float(cols[8])
                    cols[9] = float(cols[9])
                    cols[10] = float(cols[10])
                data = zip(header, cols)
                data.append(('cerna_id', cerna_id))
                data = SON(data)
                data_list.append(data)

        try:
            self.create_db_table('sg_cerna_mi2lncrna', data_list)
        except Exception as e:
            self.bind_object.set_error("导入miRNA lncRNA 对应关系信息出错!")
        else:
            self.bind_object.logger.info("导入miRNA lncRNA 对应关系信息成功!")



    @report_check
    def add_cerna(self, cerna_id, cerna):
        '''
        导入ceRNA信息
        '''
        if not isinstance(cerna_id, ObjectId):
            if isinstance(cerna_id, StringTypes):
                cerna_id = ObjectId(cerna_id)
            else:
                raise Exception("cerna_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        with open(cerna, "rb") as cerna_f:
            head_line = cerna_f.readline()
            header = head_line.strip("\n").split("\t")
            for line in cerna_f:
                cols = line.strip("\n").split("\t")
                if len(cols) >= 9:
                    cols[6] = int(cols[6])
                    cols[7] = float(cols[7])
                    cols[8] = float(cols[8])
                    cols[9] = float(cols[9])
                    cols[10] = float(cols[10])
                if len(cols) >= 16:
                    cols[14] = float(cols[14])
                    cols[15] = float(cols[15])
                data = zip(header, cols)
                data.append(('cerna_id', cerna_id))
                data = SON(data)
                data_list.append(data)

        try:
            self.create_db_table('sg_cerna_cerna', data_list)
        except Exception as e:
            self.bind_object.set_error("导入靶cerna出错!")
        else:
            self.bind_object.logger.info("导入靶cerna信息成功!")


    @report_check
    def add_cerna_nodes(self, cerna_id, nodes):
        if not isinstance(cerna_id, ObjectId):
            if isinstance(cerna_id, StringTypes):
                cerna_id = ObjectId(cerna_id)
            else:
                raise Exception("cerna_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        with open(nodes, "rb") as nodes_f:
            head_line = nodes_f.readline()
            header = head_line.strip("\n").split("\t")
            for line in nodes_f:
                cols = line.strip("\n").split("\t")
                data = zip(header, cols)
                data.append(('cerna_id', cerna_id))
                data = SON(data)
                data_list.append(data)

        try:
            self.create_db_table('sg_cerna_nodes', data_list)
        except Exception as e:
            self.bind_object.set_error("导入nodes信息出错!")
        else:
            self.bind_object.logger.info("导入nodes信息成功!")

    @report_check
    def add_cerna_edges(self, cerna_id, edges):
        if not isinstance(cerna_id, ObjectId):
            if isinstance(cerna_id, StringTypes):
                cerna_id = ObjectId(cerna_id)
            else:
                raise Exception("cerna_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        with open(edges, "rb") as edges_f:
            head_line = edges_f.readline()
            header = head_line.strip("\n").split("\t")
            for line in edges_f:
                cols = line.strip("\n").split("\t")
                cols[7] = int(cols[7])
                cols[10] = abs(float(cols[10]))
                cols[11] = abs(float(cols[11]))
                cols[12] = abs(float(cols[12]))
                cols[13] = cols[13]
                cols[14] = abs(float(cols[14]))
                cols[15] = cols[15]

                '''
                if cols[10] == 0:
                    cols[10] = 0.3
                '''
                data = zip(header, cols)
                data.append(('cerna_id', cerna_id))
                data = SON(data)
                data_list.append(data)

        try:
            self.create_db_table('sg_cerna_edges', data_list)
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
            from mbio.workflows.lnc_rna.lnc_rna_test_api import LncRnaTestApiWorkflow
            from biocluster.wsheet import Sheet
            import random
            data = {
                "id": "lnc_rna",
                #+ str(random.randint(1,10000)),
                #"id": "denovo_rna_v2",
                "project_sn": "lnc_rna",
                #+ str(random.randint(1,10000)),
                "type": "workflow",
                "name": "lnc_rna.lnc_rna_test_api",
                "options": {

                }
            }
            wsheet = Sheet(data=data)
            wf = LncRnaTestApiWorkflow(wsheet)

            wf.IMPORT_REPORT_DATA = True
            wf.IMPORT_REPORT_AFTER_END = False
            wf.test_api = wf.api.api("lnc_rna.cerna")

            params = {

            }

            test_dir = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_lnc_rna/ce_out/"
            mi2mrna = test_dir + 'mirna_mrna.corr.xls'
            mi2lncrna = test_dir + 'mirna_lncrna.corr.xls'
            ce_corr = test_dir + 'ceRNA.corr.xls'
            nodes = test_dir + 'nodes.xls'
            edges = test_dir + 'edges.xls'
            wf.test_api.import_net(test_dir, params, mi2mrna=mi2mrna, mi2lncrna=mi2lncrna, ce_corr=ce_corr,  nodes=nodes, edges=edges)
    unittest.main()

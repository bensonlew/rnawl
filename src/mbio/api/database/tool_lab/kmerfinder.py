# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
#  20210409

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re,os
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Kmerfinder(ApiBase):

    def __init__(self, bind_object):
        super(Kmerfinder, self).__init__(bind_object)
        self._project_type = "tool_lab"

    def add_kmerfinder(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "Kmerfinder" + "_"
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
                desc='Kmerfinder',
                params=params,
                status="start")
            main_id = self.create_db_table('pmlst', [main_info])
        else:
            main_id = ObjectId(main_id)
            try:
                self.update_db_record('kmerfinder', main_id)
            except Exception as e:
                self.bind_object.logger.error("导入Kmerfinder数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入Kmerfinder数据成功")
        return main_id


    def add_kmerfinder_detail(self, main_id, dir):
        data_list1 = []
        for i in os.listdir(dir):
            with open(dir + "/" + i, "r") as g:
                lines = g.readlines()

                for line in lines[1:]:
                    lin = line.strip().split("\t")
                    insert_data = {
                    "pmlst_id": ObjectId(main_id),
                    "sample_name": lin[0],
                    "ass_id": lin[1],
                    'score': lin[2],
                    "depth": lin[3],
                    "q_value": lin[4],
                    "p_value": lin[5],
                    "des": lin[6],
                    "acc_num": lin[7],
                    "taxid": lin[8],
                    "taxon": lin[9],
                    }
                    data_list1.append(insert_data)
        self.create_db_table('kmerfinder_detail', data_list1)
        table2_dict = {
            "column":[{"field":"sample_name","filter":"false","sort":"false","title":"sample name","type":"string"},
                      {"field":"ass_id","filter":"false","sort":"false","title":"AssemID","type":"string"},
                      {"field":"acc_num","filter":"false","sort":"false","title":"AccNum","type":"string"},
                      {"field":"des","filter":"false","sort":"false","title":"Description","type":"string"},
                      {"field":"taxon","filter":"false","sort":"false","title":"Taxonomy","type":"string"},
                      {"field":"taxid","filter":"false","sort":"false","title":"TaxID","type":"string"},
                      {"field":"score","filter":"false","sort":"false","title":"Score","type":"string"},
                      {"field":"depth","filter":"false","sort":"false","title":"Depth","type":"string"},
                      {"field":"q_value","filter":"false","sort":"false","title":"q_value","type":"string"},
                      {"field":"p_value","filter":"false","sort":"false","title":"p_value","type":"string"}],
            "condition": {}}
        table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
        try:
            self.update_db_record('kmerfinder', main_id, main_id=main_id, detail_table=table2_info,status="end")
        except Exception as e:
            self.bind_object.logger.error("导入kmerfinder_detail数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入kmerfinder_detail数据成功")
# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'xueqinwen'

from collections import OrderedDict
from bson import SON
from biocluster.config import Config
from io import StringIO
import types
import re
import datetime
import os
import json
from api_base import ApiBase


class AncestorAge(ApiBase):
    def __init__(self, bind_object):
        super(AncestorAge, self).__init__(bind_object)
        self._db_name = 'sanger_tool_lab'

    def add_ancestorage_detail(self, main_id, output_file):
        '''
        导入共祖时间文本数据
        '''
        ageout_file = output_file
        file_name = os.path.basename(ageout_file)
        insert_ageout_datas = []
        with open(ageout_file, 'r') as a:
            i = 0
            ageout = a.read()
            ancestor_age = ""

            all_sample = re.findall(r"All samples:\d{,}\.\d{,}", ageout)
            all_sample_time = re.findall(r"\d{,}\.\d{,}", str(all_sample[0]))
            # ancestor_age += "\n{}".format(all_sample[0])
            hp = re.findall(r't0 time:.{,}', ageout, re.M | re.I)
            ht = hp[0].decode('utf-8')
            # ft = re.finditer(r'YG\d{,}:after filter:\d{,}',ageout)
            outgroups = re.findall(r"outgroups:.{,}", ageout)
            # outgroup = re.finditer(r"YG\d{,}", outgroups[0])
            # og = ""
            # for i in outgroup:
            #     og += "{}\n".format(i.group())
            og = outgroups[0][10:]
            cp = ""
            sp = re.findall(r'.{,}-.{,}:.{,}',
                            ageout, re.M | re.I)[0].lstrip(">>>").split("-")[0]
            insert_sid_dict = {
                "type": "table",
                "ancestorage_id": self.check_objectid(main_id),
                "term": "sample_id",
                "result": "%s" % sp,
            }
            insert_sid = []
            insert_sid.append(insert_sid_dict)
            self.col_insert_data("ancestorage_detail", insert_sid)
            ag = re.finditer(r'.{,}-.{,}:.{,}', ageout, re.M | re.I)
            # ancestor_time = []
            # re.findall(r"YG\d{,}-YG\d{,}", a.group())[0],
            for a in ag:
                cp = cp + \
                    " {}".format(re.findall(r".{,}-.{,}:",
                                            a.group())[0].rstrip(":").split("-")[1])
            insert_cp_dict = {
                "type": "table",
                "ancestorage_id": self.check_objectid(main_id),
                "term": "compare_sample_id",
                "result": cp,
            }
            insert_hoplogroup_dict = {
                "type": "table",
                "ancestorage_id": self.check_objectid(main_id),
                "term": "hoplogroup",
                "result": "{}:{} years".format(ht.split(":")[1], round((float(ht.split(":")[-1])*10000),2)),
            }

            insert_outgroup_dict = {
                "type": "table",
                "ancestorage_id": self.check_objectid(main_id),
                "term": "outgroup",
                "result": og
            }

            insert_aaa_dict = {
                "type": "table",
                "ancestorage_id": self.check_objectid(main_id),
                "term": "all_ancestor_age_time",
                "result":  '{} years'.format(round((float(all_sample_time[0])*10000),2)),
            }
            insert_cp = []
            insert_cp.append(insert_cp_dict)
            insert_hoplogroup = []
            insert_hoplogroup.append(insert_hoplogroup_dict)
            insert_outgroup = []
            insert_outgroup.append(insert_outgroup_dict)
            insert_aaa = []
            insert_aaa.append(insert_aaa_dict)
            self.col_insert_data("ancestorage_detail", insert_hoplogroup)
            self.col_insert_data("ancestorage_detail", insert_cp)
            self.col_insert_data("ancestorage_detail", insert_outgroup)
            ag1 = re.finditer(r'.{,}-.{,}:.{,}', ageout, re.M | re.I)
            for a in ag1:
                insert_ancetortime = []
                insert_ancetortime_dict = {
                    "type": "table",
                    "ancestorage_id": self.check_objectid(main_id),
                }
                csample = '{}'.format(re.findall(
                    r".{,}-.{,}:", a.group())[0].rstrip(":").lstrip(">>>"))
                rt = '{} years'.format(
                    round((float(re.findall(r"\d{,}\.\d{,}", str(a.group()))[0])*10000),2))
                insert_ancetortime_dict["term"] = csample
                insert_ancetortime_dict["result"] = rt
                insert_ancetortime.append(insert_ancetortime_dict)
                self.col_insert_data("ancestorage_detail", insert_ancetortime)
            self.col_insert_data("ancestorage_detail", insert_aaa)

        # "ancestor_age": rt,
        # # "pos_for_calc": pos,
        # # "name": "{}".format(file_name)

        # ancestor_age += "\n{}".format(a.group())
        # sample1 = re.split("_","{}".format(file_name))[0]
        # insert_ageout_datas.append(insert_ageout_data)
        ancestorage_data_insert = []
        ancestorage_insert_dict = {
            "field": "term",
            "title": "",
            "sort": "false",
            "filter": "false",
            "type": "string",
        }
        ancestorage_insert_dict1 = {
            "field": "result",
            "title": "",
            "sort": "false",
            "filter": "false",
            "type": "string",
        }

        ancestorage_data_insert.append(ancestorage_insert_dict)
        ancestorage_data_insert.append(ancestorage_insert_dict1)

        ancestorage_table_dict = {
            "column": ancestorage_data_insert,
            "condition": {'type': "table"}
            # "text_data": json.dumps(ancestorage_data_insert)
        }
        update_dict = {
            "ancestor_detail_table": json.dumps(ancestorage_table_dict)
        }
        self.update_db_record(
            "ancestorage", {"_id": self.check_objectid(main_id)}, update_dict)


if __name__ == "__main__":
    a = AncestorAge(None)
    a.add_ancestorage_detail("5ea7b87917b2bf0f6184ce94",
                             "/mnt/ilustre/users/sanger-dev/workspace/20200923/AncestorAge_4g1eodvea7e5v4h1b1b062ie1b_0923172013569288_8938/output/ancestor_age/YG201901966_ancestor_age.txt")

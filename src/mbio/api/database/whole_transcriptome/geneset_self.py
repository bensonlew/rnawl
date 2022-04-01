# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.whole_transcriptome.api_base import ApiBase
import os
import pandas as pd
import unittest
from bson.objectid import ObjectId
from collections import OrderedDict
from biocluster.config import Config
from bson.objectid import ObjectId


class GenesetSelf(ApiBase):
    def __init__(self, bind_object):
        super(GenesetSelf, self).__init__(bind_object)
        self._project_type = 'whole_transcriptome'

    # @report_check
    def add_geneset(self, geneset_output_dir, name=None, level="G", source="diff", params=None,
                    project_sn='whole_transcriptome', task_id='whole_transcriptome', main_id=None):
        if params is None:
            params = dict(
                task_id=task_id,
                submit_location="geneset_upload",
                task_type="1",
                project_sn=project_sn,
                level=level,
            )
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))

        if main_id is None:
            # prepare main table info
            time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            desc = "the_geneset_upload_at_%s_by_the_customer" % time
            if name is None:
                name = desc

            main_info = dict(
                name=name,
                task_id=task_id,
                project_sn=project_sn,
                desc=desc,
                length=0,
                level=level,
                source=source,
                status="start",
                params=params,
                type="",
                is_use=0,
            )
            main_id = self.create_db_table('geneset', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        # prepare detail table info
        # target_file = os.path.join(geneset_output_dir, 'GenesetSelf.txt')
        with open(geneset_output_dir, 'r') as f:
            f.readline()
            seq_list = list()
            category_list = list()
            kind_list = list()
            for line in f:
                line = line.strip().split()
                seq_list.append(line[0])
                category_list.append(line[1])
                kind_list.append(line[2])
            length = len(seq_list)
            category = ""
            if int(category_list.count("mRNA")) >= 1:
                category += "mRNA:" + str(category_list.count("mRNA"))
            if int(category_list.count("lncRNA")) >= 1:
                if category != "":
                    category += " lncRNA:" + str(category_list.count("lncRNA"))
                else:
                    category += "lncRNA:" + str(category_list.count("lncRNA"))
            if int(category_list.count("miRNA")) >= 1:
                if category != "":
                    category += " miRNA:" + str(category_list.count("miRNA"))
                else:
                    category += "miRNA:" + str(category_list.count("miRNA"))
            if int(category_list.count("circRNA")) >= 1:
                if category != "":
                    category += " circRNA:" + str(category_list.count("circRNA"))
                else:
                    category += "circRNA:" + str(category_list.count("circRNA"))
            # category = "mRNA: " + str(category_list.count("mRNA")) + ", lncRNA: " + str(category_list.count("lncRNA"))\
            #            + ", miRNA: " + str(category_list.count("miRNA")) + ", circRNA: " + str(category_list.count("circRNA"))
        if length < 5:
            self.db['geneset'].update({'_id': main_id}, {'$set': {'params': None, 'status': 'failed'}})
            self.bind_object.set_error("预上传基因集中基因与该项目中基因匹配数目过少（＜5个），不予上传，请核查ID书写是否规范或者从交互页面直接创建")
        row_dict_list = [{"seq_list": seq_list, "category_list": category_list, "kind_list": kind_list}]
        # insert detail
        tag_dict = dict(geneset_id=main_id)
        self.create_db_table('geneset_detail', row_dict_list, tag_dict=tag_dict)
        self.update_db_record('geneset', main_id, status="end", main_id=main_id, length=length, type=category)


################################################
class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test(test):
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
            "id": "whole_transcriptome",
            "project_sn": "whole_transcriptome",
            "type": "workflow",
            "name": "whole_transcriptome_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.api_geneset_self = wf.api.api("whole_transcriptome.geneset_self")

        wf.api_geneset_self.add_geneset(
            geneset_output_dir="/mnt/ilustre/users/sanger-dev/workspace/20191012/GenesetSelf_whole_transcriptome_5272_2227/test.txt",
            level='T', main_id='5da13c4717b2bf14e012f752')


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
# __author__ = 'konghualei, 20170421'
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from bson import ObjectId
import datetime
import shutil
import re,os
import time
from biocluster.workflow import Workflow
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class CernaSankeyWorkflow(Workflow):
    """
    基因集venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CernaSankeyWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "file_id", "type": "string"},
            {"name": "task_type", "type": "int"},
            {"name": "submit_location", "type": "string"},
            {"name": "relation_type", "type": "string"},
            {"name": "relation_file", 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {"name": "relation_string", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.start_listener()
        self.fire("start")
        self.get_run_log()
        relations = []
        if self.option("relation_type") == "file":
            with open(self.option("relation_file").prop['path'], 'r') as f:
                f.readline()
                for line in f:
                    relation = line.strip().split("\t")
                    if len(relation) != 3:
                        self.set_error("文件格式错误需要为三列")
                    else:
                        relations.append(relation)
        elif self.option("relation_type") == "string":
            for line in self.option("relation_string").split("|"):
                relation = line.strip().split(",")
                if len(relation) != 3:
                    self.set_error("输入格式错误需要为三列")
                else:
                    relations.append(relation)

        sankey_file = self.relation2sankey(relations)

        self.set_db(sankey_file)
        self.end()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="cerna_sankey", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def relation2sankey(self, relations):
        # 转变输入表格为sankey图 js 元素
        mirna_w = dict()
        base_w = 1.5
        lnc2mi_w = dict()
        mi2lnc_num = dict()
        uniq2 = set()
        for relation in relations:
            if relation[1] + "|" + relation[2] in uniq2:
                pass
            else:
                if relation[1] in mirna_w:
                    mirna_w[relation[1]] = mirna_w[relation[1]] + base_w
                else:
                    mirna_w[relation[1]] = base_w
                uniq2.add(relation[1] + "|" + relation[2])

            if relation[1] in mi2lnc_num:
                mi2lnc_num[relation[1]].add(relation[0])
            else:
                mi2lnc_num[relation[1]] = set([relation[0]])
        print mirna_w

        uniq_w1 = set()
        uniq_w2 = set()
        with open(self.output_dir + "/cerna_sankey.xls", 'w') as fo:
            fo.write("type\tnode1\tnode2\twidth\n")
            for relation in relations:
                if relation[0] + "|" + relation[1] in uniq_w1:
                    pass
                else:
                    w1 = mirna_w[relation[1]]/len(mi2lnc_num[relation[1]])
                    fo.write("nc2mirna\t{}\t{}\t{}\n".format(relation[0], relation[1], w1))
                    uniq_w1.add(relation[0] + "|" + relation[1])
            for relation in relations:
                if relation[1] + "|" + relation[2] in uniq_w2:
                    pass
                else:
                    fo.write("mi2mrna\t{}\t{}\t{}\n".format(relation[1], relation[2], base_w))
                    uniq_w2.add(relation[1] + "|" + relation[2])
        return self.output_dir + "/cerna_sankey.xls"


    def end(self):
        super(CernaSankeyWorkflow, self).end()

    def set_db(self, sankey_file):
        cerna_api = self.api.api("whole_transcriptome.cerna")
        cerna_api.import_cerna_sankey(
            self.option('main_id'),
            sankey_file
        )

# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modifiy = modified 2019.02.21

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.workflows.denovo_rna_v2.denovo_test_api import DenovoTestApiWorkflow
from biocluster.wsheet import Sheet
import os
import re
import shutil
import json
import glob
from biocluster.config import Config
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class TargetCisWorkflow(Workflow):
    """
    交互分析进行筛选blast参数重注释时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(TargetCisWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "novol", "type": "infile", "format": "lnc_rna.gtf"},  # 输入文件
            {"name": "known", "type": "infile", "format": "lnc_rna.gtf"},  # 输入文件
            {'type': 'infile', 'name': 'mrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'annotation', 'format': 'lnc_rna.common'},
            {'type': 'string', 'name': 'geneset_id', 'default': 'All'},
            {'type': 'int', 'name': 'up_dis', 'default': 10},
            {'type': 'int', 'name': 'down_dis', 'default': 10},
            {"name": "stat_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "origin_result", "type": "string"},
            {"name": "target_cis_id", "type": "string"},
            {"name": "last_id_target", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "int"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.target_predict_known = self.add_tool("lnc_rna.lnc_target_cis")
        self.target_predict_novol = self.add_tool("lnc_rna.lnc_target_cis")
        self.on_rely([self.target_predict_known, self.target_predict_novol], self.set_output)
        # self.target_predict.on('end', self.set_db)

    def run_known(self):
        options = {
            "lncrna_gtf" : self.option("known"),
            "mrna_gtf" : self.option("mrna_gtf"),
            "annotation" : self.option("annotation"),
            "up_dis": self.option("up_dis"),
            "down_dis": self.option("down_dis"),
        }
        self.target_predict_known.set_options(options)
        self.target_predict_known.run()

    def run_novol(self):
        options = {
            "lncrna_gtf" : self.option("novol"),
            "mrna_gtf" : self.option("mrna_gtf"),
            "annotation" : self.option("annotation"),
            "up_dis": self.option("up_dis"),
            "down_dis": self.option("down_dis"),

        }
        self.target_predict_novol.set_options(options)
        self.target_predict_novol.run()


    def run(self):
        self.get_run_log()
        self.run_known()
        self.run_novol()
        super(TargetCisWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_target_cis", main_id=self.option('target_cis_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_output(self):
        fname = os.path.join(self.target_predict_known.output_dir, 'lnc_rna_cistarget.annot.xls')
        link = os.path.join(self.output_dir, "known_lncrna_cistarget.annot.xls")
        if os.path.exists(link):
            os.remove(link)
        os.link(fname, link)

        fname = os.path.join(self.target_predict_novol.output_dir, 'lnc_rna_cistarget.annot.xls')
        link = os.path.join(self.output_dir, "novol_lncrna_cistarget.annot.xls")
        if os.path.exists(link):
            os.remove(link)
        os.link(fname, link)

        self.set_db()


    def set_db(self):
        self.logger.info("保存结果到mongo")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.test_api = self.api.api("lnc_rna.lnc_target")

        self.test_api.import_target_cis_web(
            self.option("target_cis_id"),
            os.path.join(self.output_dir, "known_lncrna_cistarget.annot.xls"),
            os.path.join(self.output_dir, "novol_lncrna_cistarget.annot.xls")
        )
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        target_dir = self.output_dir
        repaths = [
            [".", "", "靶基因"],
            ["known_lncrna_cistarget.annot.xls", "", "已知miRNA 对应的靶基因详情表"],
            ["novol_lncrna_cistarget.annot.xls", "", "新miRNA 对应的靶基因详情表"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ]

        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        elif re.match(r'sanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        elif re.match(r'^\w+://\S+/.+$', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp
        else:
            self.set_error("json output wrong")

        result_dir = self.add_upload_dir(target_dir)
        result_dir.add_regexp_rules(repaths)
        '''
        db = Config().get_mongo_client(mtype="small_rna")[Config().get_mongo_dbname("small_rna")]
        # col1 = db["sg_annotation_stat"]
        # col1.update({"_id": ObjectId(self.option("stat_id"))}, {"$set": {"result_dir": self.workflow_output + "/Annotation"}}, upsert=True)
        col2 = db["sg_target"]
        col2.update({"_id": ObjectId(self.option("target_id"))}, {"$set": {"result_dir": self.workflow_output}}, upsert=True)
        '''

        super(TargetCisWorkflow, self).end()

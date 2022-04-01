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


class TargetTransWorkflow(Workflow):
    """
    交互分析进行筛选blast参数重注释时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(TargetTransWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "novol", "type": "infile", "format": "lnc_rna.fasta"},  # 输入文件
            {"name": "known", "type": "infile", "format": "lnc_rna.fasta"},  # 输入文件
            {'type': 'infile', 'name': 'mrna', 'format': 'lnc_rna.fasta'},
            {"name": "novol_gtf", "type": "infile", "format": "lnc_rna.gtf"},  # 输入文件
            {"name": "known_gtf", "type": "infile", "format": "lnc_rna.gtf"},  # 输入文件
            {'type': 'infile', 'name': 'mrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'annotation', 'format': 'lnc_rna.common'},

            {'type': 'string', 'name': 'lnc_geneset_id', 'default': 'All'},
            {'type': 'string', 'name': 'lnc_geneset_type', 'default': 'LG'},
            {'type': 'string', 'name': 'm_geneset_type', 'default': 'G'},
            {'type': 'string', 'name': 'm_geneset_id', 'default': 'All'},
            {'type': 'string', 'name': 'method', 'default': 'rnaplex'},
            {"name": "stat_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "origin_result", "type": "string"},
            {"name": "target_trans_id", "type": "string"},
            {"name": "last_id_target", "type": "string", "default": None},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "int"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.target_predict_known = self.add_module("lnc_rna.lnc_target_trans")
        self.target_predict_novol = self.add_module("lnc_rna.lnc_target_trans")

        # self.target_predict.on('end', self.set_db)
        self.choosed_known = True
        self.choosed_novol = True

    def choose_seq(self):
        '''
        根据列表筛选fasta
        '''
        if self.option("lnc_geneset_id") != "All":
            lnc_list = []
            with open(self.option("lnc_geneset_id"), 'r') as f:
                for line in f:
                    lnc_list.append(line.strip())

            if self.option("lnc_geneset_type") == "LG":
                lnc_tran2gene = self.option("known_gtf").get_txpt_gene_dic()
                lnc_tran2gene.update(self.option("novol_gtf").get_txpt_gene_dic())
                lnc_list = [k for k,v in lnc_tran2gene.items() if v in lnc_list]

            choosed = self.option("novol").choose_seq_by_list(lnc_list, self.work_dir + "/choose_novol_lnc.fa")
            if len(choosed) > 0:
                self.choosed_novol = True
                self.option("novol", self.work_dir + "/choose_novol_lnc.fa")
            else:
                self.choosed_novol = False

            # self.logger.info("lnc_list is {}".format(lnc_list))
            choosed = self.option("known").choose_seq_by_list(lnc_list, self.work_dir + "/choose_known_lnc.fa")
            if len(choosed) > 0:
                self.choosed_known = True
                self.option("known", self.work_dir + "/choose_known_lnc.fa")
            else:
                self.choosed_known = False

        if self.option("m_geneset_id") != "All":
            m_list = []
            with open(self.option("m_geneset_id"), 'r') as f:
                for line in f:
                    m_list.append(line.strip())
            if self.option("m_geneset_type") == "G":
                m_tran2gene = self.option("mrna_gtf").get_txpt_gene_dic()
                m_list = [k for k,v in m_tran2gene.items() if v in m_list]

            self.option("mrna").choose_seq_by_list(m_list, self.work_dir + "/choose_mrna.fa")
            self.option("mrna", self.work_dir + "/choose_mrna.fa")

    def run_known(self):
        options = {
            "query" : self.option("known"),
            "target": self.option("mrna"),
            "method": self.option("method"),
            "gtf": self.option("mrna_gtf"),
            "lnc_gtf": self.option("known_gtf"),
            "annotation": self.option("annotation"),
        }
        self.target_predict_known.set_options(options)
        self.target_predict_known.run()

    def run_novol(self):
        options = {
            "query" : self.option("novol"),
            "target": self.option("mrna"),
            "method": self.option("method"),
            "gtf": self.option("mrna_gtf"),
            "lnc_gtf": self.option("novol_gtf"),
            "annotation": self.option("annotation")
        }
        self.target_predict_novol.set_options(options)
        self.target_predict_novol.run()


    def run(self):
        self.get_run_log()
        self.choose_seq()
        if self.choosed_known and self.choosed_novol:
            self.on_rely([self.target_predict_known, self.target_predict_novol], self.set_output)
        elif self.choosed_known:
            self.target_predict_known.on('end', self.set_output)
        elif self.choosed_novol:
            self.target_predict_novol.on('end', self.set_output)
        if self.choosed_known:
            self.run_known()
        if self.choosed_novol:
            self.run_novol()
        super(TargetTransWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_target_trans", main_id=self.option('target_trans_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_output(self):
        if self.choosed_known:
            known_file = glob.glob(self.target_predict_known.output_dir + '/*merge_out.annot.xls')
            fname = known_file[0]
            link = os.path.join(self.output_dir, "known_lncrna_transtarget.annot.xls")
            if os.path.exists(link):
                os.remove(link)
            os.link(fname, link)
        if self.choosed_novol:
            novol_file = glob.glob(self.target_predict_novol.output_dir + '/*merge_out.annot.xls')
            fname = novol_file[0]
            link = os.path.join(self.output_dir, "novol_lncrna_transtarget.annot.xls")
            if os.path.exists(link):
                os.remove(link)
            if os.path.exists(fname):
                os.link(fname, link)
            else:
                print fname
        self.set_db()


    def set_db(self):
        self.logger.info("保存结果到mongo")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.test_api = self.api.api("lnc_rna.lnc_target")

        self.test_api.import_target_trans_web(
            self.option("target_trans_id"),
            self.option("last_id_target"),
            os.path.join(self.output_dir, "novol_lncrna_transtarget.annot.xls"),
            os.path.join(self.output_dir, "known_lncrna_transtarget.annot.xls")
        )
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        target_dir = self.output_dir

        repaths = [
            [".", "", "靶基因"],
            ["known_lncrna_transtarget.annot.xls", "", "已知miRNA 对应的靶基因详情表"],
            ["novol_lncrna_transtarget.annot.xls", "", "新miRNA 对应的靶基因详情表"],
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

        super(TargetTransWorkflow, self).end()

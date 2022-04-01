# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modifiy = modified 2019.09.15

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
from biocluster.core.function import filter_error_info, link, CJsonEncoder


class TargetCistransWorkflow(Workflow):
    """
    交互分析进行筛选blast参数重注释时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(TargetCistransWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "novol", "type": "infile", "format": "lnc_rna.gtf"},  # 输入文件
            {"name": "known", "type": "infile", "format": "lnc_rna.gtf"},  # 输入文件
            {'type': 'infile', 'name': 'mrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'annotation', 'format': 'lnc_rna.common'},
            {'name': 'exp_matrix_lnc', 'type': 'infile', 'format': 'lnc_rna.express_matrix'},
            {'name': 'exp_matrix_target', 'type': 'infile', 'format': 'lnc_rna.express_matrix'},
            {'name': 'pvalue_cutoff', 'type': 'float', 'default': 1},
            {'name': 'qvalue_cutoff', 'type': 'float', 'default': 1},
            {'name': 'cor_cutoff', 'type': 'float', 'default': 0.9},
            {'name': 'corr_way', 'type': 'string', 'default': "spearmanr"},
            {'name': 'padjust_way', 'type': 'string', 'default': "fdr_bh"},
            {'type': 'string', 'name': 'geneset_id', 'default': 'All'},
            {'type': 'int', 'name': 'up_dis', 'default': 10},
            {'type': 'int', 'name': 'down_dis', 'default': 10},
            {"name": "stat_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "origin_result", "type": "string"},
            {"name": "lncset_id", "type": "string"},
            {"name": "targetset_id", "type": "string"},
            {"name": "target_cis_id", "type": "string"},
            {"name": "last_id_target", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "int"},
            {"name": "update_info", "type": "string"},
            {"name": "group_dict", "type": "string"},
            {"name": "group_id", "type": "string"},

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.target_predict_known = self.add_tool("lnc_rna.lnc_target_cis")
        self.target_predict_novol = self.add_tool("lnc_rna.lnc_target_cis")
        self.merge_target = self.add_tool("lnc_rna.lnc_target_merge")
        self.expcorr_tool = self.add_tool("lnc_rna.expcorr_lnctarget")
        # self.target_predict.on('end', self.set_db)
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/03 LncRNA_Target')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(TargetCistransWorkflow, self).send_log(data)

    def run_exp_corr(self):
        options = dict(
            exp=self.option('exp_matrix_lnc'),
            exp_t=self.option('exp_matrix_target'),
            gt="G",
            anno=self.option('annotation').prop['path'],
            pvalue_cutoff=self.option('pvalue_cutoff'),
            qvalue_cutoff=self.option('qvalue_cutoff'),
            cor_cutoff=self.option('cor_cutoff'),
            corr_way=self.option('corr_way'),
            padjust_way=self.option('padjust_way'),
        )
        self.expcorr_tool.set_options(options)
        self.expcorr_tool.run()


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

    def run_merge_target(self):
        options = {
            "target_gtf" : self.option("mrna_gtf"),
            "known_lnc_gtf" : self.option("known"),
            "novel_lnc_gtf" : self.option("novol"),
            "cis_known" : self.target_predict_known.output_dir + "/lnc_rna_cistarget.annot.xls",
            "cis_novol" : self.target_predict_novol.output_dir + "/lnc_rna_cistarget.annot.xls",
            "annotation" : self.option("annotation"),
            "exp_corr": self.expcorr_tool.output_dir + "/corr.xls",

        }
        self.merge_target.set_options(options)
        self.merge_target.run()



    def run(self):
        self.on_rely([self.expcorr_tool, self.target_predict_known, self.target_predict_novol], self.run_merge_target)
        self.merge_target.on('end', self.set_output)
        self.get_run_log()
        self.run_exp_corr()
        self.run_known()
        self.run_novol()
        super(TargetCistransWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_target_cistrans", main_id=self.option('target_cis_id'),
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

        self.test_api.import_target_cistrans_web(
            self.option("target_cis_id"),
            os.path.join(self.output_dir, self.merge_target.output_dir + '/cis_annot.xls'),
            os.path.join(self.output_dir, self.merge_target.output_dir + '/trans_annot.xls')
        )
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.merge_target.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.merge_target.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.merge_target.output_dir, os.path.basename(self.run_log)))
        target_dir = self.merge_target.output_dir

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
        self.inter_dirs = [
            ["03 LncRNA_Target", "", "lncRNA靶基因预测结果目录", 0],
        ]
        result_dir.add_relpath_rules([
            [".", '', '靶基因预测文件', 0],
            ["cis_annot.xls", 'xls', 'cis作用靶基因预测详情表',0],
            ['trans_annot.xls', 'xls', 'trans作用靶基因预测详情表',0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])

        '''
        db = Config().get_mongo_client(mtype="small_rna")[Config().get_mongo_dbname("small_rna")]
        # col1 = db["sg_annotation_stat"]
        # col1.update({"_id": ObjectId(self.option("stat_id"))}, {"$set": {"result_dir": self.workflow_output + "/Annotation"}}, upsert=True)
        col2 = db["sg_target"]
        col2.update({"_id": ObjectId(self.option("target_id"))}, {"$set": {"result_dir": self.workflow_output}}, upsert=True)
        '''

        super(TargetCistransWorkflow, self).end()

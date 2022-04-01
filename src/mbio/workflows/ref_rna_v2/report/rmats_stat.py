# -*- coding: utf-8 -*-
# __author__ = fiona
# last modified by shicaiping at 20180515

import glob, os, Bio, shutil, argparse, sys, fileinput, urllib2
from bson import ObjectId
import re
from biocluster.workflow import Workflow
from biocluster.config import Config
from mbio.packages.ref_rna_v2.rmats_process_func import process_rmats_stat
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.ref_rna_v2.chart import Chart
import glob

class RmatsStatWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        self.new_options = {}
        super(RmatsStatWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "pvalue_fdr", "type": "string"},
            {"name": "rmats_out_root_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {"name": "fdr", "type": "string"},
            {"name": "psi", "type": "string"},
            {"name": "stat_id", 'type': "string"},
            {"name": "update_info", 'type': "string"},
            {"name": "task_id", 'type': "string"},
        ]
        self.logger.info(options)
        self.add_option(options)
        self.set_options(self._sheet.options())
    
    def run_rmats_stat(self):
        self.result_dir = self.work_dir + "/result_dir"
        if os.path.exists(self.result_dir):
            shutil.rmtree(self.result_dir)
        os.mkdir(self.result_dir)
        files = glob.glob(self.option("rmats_out_root_dir").prop['path'] + '/*.alter_id.txt')
        fdr = float(self.option("fdr"))
        psi = float(self.option("psi"))
        for file in files:
            os.link(file, self.result_dir + "/" + os.path.basename(file))
        process_rmats_stat(root=self.result_dir, pvalue_fdr=self.option('pvalue_fdr'), fdr=fdr, psi=psi)
        for file in os.listdir(self.result_dir):
            if file in ("psi_stats.file.txt", "event_stats.file.txt"):
                os.link(os.path.join(self.result_dir,file), os.path.join(self.output_dir, file))

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        splice_diff = self.output_dir + "/event_stats.file.txt"
        splice_psi = self.output_dir + "/psi_stats.file.txt"
        chart.chart_splice_diff_stat(splice_diff, splice_psi, cmp_name="")
        chart.to_pdf()

        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.output_dir +"/" + os.path.basename(p))

        # move pdf to result dir

    def end(self):
        # self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        relpath = [
            [".", "", "差异可变剪切事件统计结果目录", 0, "211224"],
            ["event_stats.file.txt", "", "差异可变剪切事件统计表", 0, "211225"],
            ["psi_stats.file.txt", "", "差异可变剪切模式变化统计表", 0, "211226"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ['*.pdf', 'txt', "差异可变剪切统计图", 0],
        ]
        result_dir.add_relpath_rules(relpath)
        super(RmatsStatWorkflow, self).end()

    def set_db(self):
        """
        保存结果表保存到mongo数据库中
        """
        self.logger.info("开始导表")
        self.workflow_output_tmp = self._sheet.output
        self.api_as = self.api.api("ref_rna_v2.splicing_rmats")
        self.api_as.add_sg_rmats_stat_for_controller(stat_id=ObjectId(self.option("stat_id")), outpath=self.output_dir)
        self.api_as.db['sg_splicing_rmats_stats'].update({'_id': ObjectId(self.option('stat_id'))},
                                                 {'$set': {'main_id': ObjectId(self.option('stat_id')), "status": "end"}}, upsert=True)

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_splicing_rmats_stats", main_id=self.option('stat_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run(self):
        self.start_listener()
        self.fire("start")
        self.get_run_log()
        self.run_rmats_stat()
        self.set_db()
        self.end()

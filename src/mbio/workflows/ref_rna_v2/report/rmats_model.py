# -*- coding: utf-8 -*-
# __author__ = fiona
# last modified by shicaiping at 20180515

import re, os, Bio, argparse, sys, fileinput, urllib2
from bson import ObjectId
import re,shutil
from biocluster.workflow import Workflow
from biocluster.config import Config
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.ref_rna_v2.chart_advance import ChartAdvance


class RmatsModelWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        self.new_options = {}
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        super(RmatsModelWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gene_id", "type": "string"},
            {"name": "rmats_out_root_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {"name": "label_a", "type": "string"},
            {"name": "label_b", "type": "string"},
            {"name": "event_type", 'type': "string"},
            {"name": "as_type", 'type': "string"},
            {"name": "intron_s", 'type': "int", "default": 1},
            {"name": "exon_s", 'type': "int", "default": 1},
            {"name": "splicing_id", 'type': "string"},
            {"name": "update_info", 'type': "string"},
            {"name": "task_id", 'type': "string"},
        ]
        self.logger.info(options)
        self.add_option(options)
        self.set_options(self._sheet.options())
    
    def run_rmats_model(self):
        a_bam_file = os.path.join(self.option('rmats_out_root_dir').prop['path'], "A_group_bam.txt")
        b_bam_file = os.path.join(self.option('rmats_out_root_dir').prop['path'], "B_group_bam.txt")
        if "interaction_results" in self.workflow_output:
            tmp = self.workflow_output.split("interaction_results")[0].replace("s3nb1","s3nb")
            bam_prefix = os.path.join(tmp, "workflow_results/Align/AlignBam")
        else:
            raise OptionError("output dir is not correct!")
        with open (a_bam_file, "r") as f1:
            a_bam_str = f1.readline().strip()
        with open (b_bam_file, "r") as f2:
            b_bam_str = f2.readline().strip()
        label_a = ""
        a_bam_str_new = ""
        b_bam_str_new = ""
        if os.path.exists(self.work_dir + "/bam_path"):
            shutil.rmtree(self.work_dir + "/bam_path")
        os.mkdir(self.work_dir + "/bam_path")
        transfer = MultiFileTransfer()
        for each in a_bam_str.split(","):
            if not os.path.exists(each):
                bam_basename = os.path.basename(each)
                transfer.add_download(bam_prefix + "/" + os.path.basename(each), self.work_dir + "/bam_path/")
                transfer.perform()
                if a_bam_str_new == "":
                    a_bam_str_new = self.work_dir + "/bam_path/" + bam_basename
                else:
                    a_bam_str_new += "," + self.work_dir + "/bam_path/" + bam_basename
            else:
                a_bam_str_new = a_bam_str
            if label_a:
                label_a = label_a + "," + each.split("/")[-1].replace(".bam","")
            else:
                label_a = each.split("/")[-1].replace(".bam","")
        label_b = ""
        for each in b_bam_str.split(","):
            if not os.path.exists(each):
                bam_basename = os.path.basename(each)
                transfer.add_download(bam_prefix + "/" + os.path.basename(each), self.work_dir + "/bam_path/")
                transfer.perform()
                if b_bam_str_new == "":
                    b_bam_str_new = self.work_dir + "/bam_path/" + bam_basename
                else:
                    b_bam_str_new += "," + self.work_dir + "/bam_path/" + bam_basename
            else:
                b_bam_str_new = b_bam_str
            if label_b:
                label_b = label_b + "," + each.split("/")[-1].replace(".bam","")
            else:
                label_b = each.split("/")[-1].replace(".bam","")
        rmats_event_file = self.option('rmats_out_root_dir').prop["path"] + "/" + self.option('event_type') + ".MATS." + self.option('as_type') + '.alter_id.txt'
        event_file = os.path.join(self.work_dir, "event_file.txt")
        event_num = 0
        with open (rmats_event_file, "r") as r, open(event_file, "w") as w:
            head = r.readline()
            w.write(head)
            for line in r:
                if line.split("\t")[1] == self.option('gene_id'):
                    w.write(line)
                    event_num += 1
        self.logger.info(event_file)

        if (event_num >= 1):
            self.new_options = {
                "a_bam_str": a_bam_str_new,
                "b_bam_str": b_bam_str_new,
                "label_a": label_a,
                "label_b": label_b,
                "event_type": self.option("event_type"),
                "event_file": event_file,
                "intron_s": self.option("intron_s"),
                "exon_s": self.option("exon_s")
            }
            self.rmats_model = self.add_tool("ref_rna_v2.rmats_model")
            self.rmats_model.set_options(self.new_options)
            self.rmats_model.on('end', self.set_db)
            self.rmats_model.run()
        else:
            self.set_error('没有可变剪切事件可画模式图')
            # db = Config().get_mongo_client(mtype="ref_rna_v2")[Config().get_mongo_dbname("ref_rna_v2")]
            # conn = db["sg_splicing_rmats_model"]
            # conn.update({"task_id": self.option('task_id')}, {"$set": {'graph_dir': "", 'pdf_files' : "", 'png_files': "", 'status': 'end'}}, upsert=True)
            # self.end()
    
    def end(self):
        if os.path.exists(os.path.join(self.rmats_model.output_dir + '/Sashimi_plot', os.path.basename(self.run_log))):
            os.remove(os.path.join(self.rmats_model.output_dir + '/Sashimi_plot', os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.rmats_model.output_dir + '/Sashimi_plot', os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.rmats_model.output_dir + '/Sashimi_plot')
        relpath = [
            [".", "", "可变剪切模式图", 0, "211223"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ]
        result_dir.add_relpath_rules(relpath)
        super(RmatsModelWorkflow, self).end()
    
    def set_db(self):
        """
        保存结果表保存到mongo数据库中
        """
        db = Config().get_mongo_client(mtype="ref_rna_v2")[Config().get_mongo_dbname("ref_rna_v2")]
        pdf_files = [f for f in os.listdir(self.rmats_model.output_dir + '/Sashimi_plot') if re.match('^\S+\.(pdf)$', f.strip())]
        png_files = [f for f in os.listdir(self.rmats_model.output_dir + '/Sashimi_plot') if re.match('^\S+\.(png)$', f.strip())]
        graph_dir = self.workflow_output
        self.logger.info("rmats model图文件夹为： %s" % graph_dir)
        conn = db["sg_splicing_rmats_model"]
        self.main_id = conn.find_one({"task_id": self.option('task_id')})["_id"]
        conn.update({"task_id": self.option('task_id')}, {"$set": {'graph_dir': graph_dir, 'pdf_files' : pdf_files, 'png_files': png_files, 'status': 'end'}}, upsert=True)
        self.end()

    def get_run_log(self):
        db = Config().get_mongo_client(mtype="ref_rna_v2")[Config().get_mongo_dbname("ref_rna_v2")]
        conn = db["sg_splicing_rmats_model"]
        self.main_id = conn.find_one({"task_id": self.option('task_id')})["_id"]
        get_run_log = GetRunLog("ref_rna_v2", table="sg_splicing_rmats_model", main_id=self.main_id,
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()
    
    def run(self):
        self.get_run_log()
        self.run_rmats_model()
        super(RmatsModelWorkflow, self).run()

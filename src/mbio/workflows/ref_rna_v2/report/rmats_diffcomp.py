# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from bson.objectid import ObjectId
import os
import re
from biocluster.file import getsize, exists
from biocluster.file import download
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.ref_rna_v2.chart import Chart


class RmatsDiffcompWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RmatsDiffcompWorkflow,self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type":"string", "default": None},
            {"name": "task_id", "type":"string", "default": None},
            {"name": "submit_location", "type":"string", "default": None},
            {"name": "task_type", "type":"string", "default": None},
            {"name": "main_id", "type":"string", "default": None},
            {"name": "rmats_detail", "type": "string", "default": None},
            {"name": "significant_value", "type": "string", "default": "0.05"},
            {"name": "delta_PSI", "type": "string", "default": "0.0"},
            {"name": "significant_diff", "type": "string", "default": "FDR"},
            {"name": "last_id", "type": "string", "default": None},
            {"name": "compare_plan", "type": "string", "default": None}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.rmasts_diffcomp = self.add_tool("ref_rna_v2.rmats_uniqid")
        self.dump_tool = self.api.api("ref_rna_v2.rmats_diffcomp")

    def run(self):
        self.rmasts_diffcomp.on("end", self.set_db)
        self.get_run_log()
        self.download_rmat_detail()
        self.run_rmats_diffcomp()
        super(RmatsDiffcompWorkflow,self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_splicing_rmats_diffcomp", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        elif os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find')
        return to_path

    def download_rmat_detail(self):
        rmat_download = list()
        compares = self.option("compare_plan").split(",")
        i = 0
        for rmat in self.option("rmats_detail").split(","):
            rmat_new = self.download_s3_file(rmat.strip(), os.path.basename(rmat.strip()) + "_" + compares[i])
            rmat_download.append(rmat_new)
            i += 1
        self.option("rmats_detail", ",".join(rmat_download))

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"

        a.chart_geneset_enrich_circ(sys.argv[1], sys.argv[2])
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")[0]
        os.link(pdf_file, self.tool.output_dir + "/sample_correlation.pdf")

    def set_db(self):
        workflow_output = self.get_workflow_output_dir()
        self.dump_tool.add_rmats_diffcomp_detail(self.rmasts_diffcomp.output_dir,self.option("main_id"))
        self.dump_tool.update_db_record('sg_splicing_rmats_diffcomp', self.option('main_id'), output_dir=workflow_output, status="end")
        if self.option('last_id') == "":
            pass
        else:
            self.dump_tool.remove_table_by_main_record("sg_splicing_rmats_diffcomp", _id =ObjectId(self.option('last_id')), detail_table="sg_splicing_rmats_diffcomp_detail", detail_table_key="rmats_diffcomp_id")
        self.end()

    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.rmasts_diffcomp.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.rmasts_diffcomp.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.rmasts_diffcomp.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.rmasts_diffcomp.output_dir)
        result_dir.add_relpath_rules([
            [".","","RmastDiffcomp", 0, "211222"]
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(RmatsDiffcompWorkflow,self).end()

    def run_rmats_diffcomp(self):
        options = {
            "rmats_detail" : self.option('rmats_detail'),
            "significant_value" : float(self.option('significant_value')),
            "delta_PSI" : float(self.option('delta_PSI')),
            "significant_diff" : self.option('significant_diff'),
            "compare_plan" : self.option('compare_plan')
        }
        self.rmasts_diffcomp.set_options(options)
        self.rmasts_diffcomp.run()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output

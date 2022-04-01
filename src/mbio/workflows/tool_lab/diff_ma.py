# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import glob
import unittest
import types
import json
from biocluster.config import Config
from biocluster.file import getsize, exists
from biocluster.file import download
from bson.objectid import ObjectId
import time


class DiffMaWorkflow(Workflow):
    """
    Used for cds to protein code
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffMaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "raw_file", "type": "infile", "format": "ref_rna_v2.common"},  # fasta文件
            # {"name": "raw_file_info", "type": "string"},  # fasta文件
            {"name": "project_type", "type": "string","default": None},  # fasta文件
            {"name": "pvalue", "type": "float", "default": 0.05},
            {"name": "fc", "type": "float", "default": 2.0},
            {"name": "x_axis_name", "type": "string", "default": "log10(TPM)"},
            {"name": "y_axis_name", "type": "string", "default": "-log10(Padjust)"},
            {"name": "title_name", "type": "string","default": "Volcano Plot"},
            {"name": "color", "type": "string", "default": "red_blue_grey"},
            {"name": "source", "type": "string", "default": "tool_lab"},
            {"name": "main_id", "type": "string"},
            {'name': 'project_task_id', 'type': 'string'},
            {'name': 'relate_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.diff_ma")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(DiffMaWorkflow, self).run()

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
            self.set_error('file can not find %s', variables=(path), code = '13700502')
        return to_path

    def check_file_path(self):
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        conn_upset = db[collection_name]
        status = 'start'
        count_time = 0
        while status == 'start':
            if count_time > 600:
                self.set_error('超过十分钟还没有结果文件生成，请检查是否生成文件时报错')
                break
            time.sleep(10)
            print('sleep 10s')
            try:
                upset = conn_upset.find_one(
                    {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
                status = upset['status']
            except:
                pass
            count_time += 10
        upset = conn_upset.find_one(
            {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
        file_path = upset['file_path']
        return file_path

    def get_raw_file_info(self):
        if not self.option("project_type"):
            raw_file_path = self.option("raw_file").prop["path"]
            project_type = "custom"
            return project_type, raw_file_path
        elif  self.option("project_type") == "refrna" :
            raw_file_path = self.option("raw_file").prop["path"]
            project_type = "ref_rna_v2"
            return project_type, raw_file_path
        else:
            self.file_path = self.check_file_path()
            raw_file_path = self.download_s3_file(self.file_path, 'diff_ma.txt')
            project_type = "custom"
            return project_type, raw_file_path

            # if self.option('source') == 'project':
            #     self.file_path = self.check_file_path()
            #     raw_file_path = self.download_s3_file(self.file_path, 'diff_ma.txt')
            #     project_type = "custom"
            #     return project_type,raw_file_path
            # else:
            #     raw_file_path = self.option("raw_file").prop["path"]
            #     project_type = "custom"
            #     return project_type,raw_file_path


        # try:
        #     raw_file_info = json.loads(self.option("raw_file_info"))
        #     project_type = raw_file_info["project_type"]
        #     raw_file_path = raw_file_info["file_path"]
        #     return project_type,raw_file_path
        # except:
        #     project_type = "custom"
        #     raw_file_path = self.option("raw_file_info")
        #     return project_type,raw_file_path


    def run_tool(self):
        project_type,file_path = self.get_raw_file_info()
        opts = {
            'raw_file': file_path,
            'project_type' :project_type,
            'pvalue': self.option('pvalue'),
            'fc': self.option('fc'),
            'x_axis_name' : self.option("x_axis_name"),
            'y_axis_name': self.option("y_axis_name"),
            'title_name': self.option("title_name"),
            'color': self.option("color"),
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        output_dir = self.tool.output_dir
        pdf_path = glob.glob(os.path.join(output_dir, '*.pdf'))
        pdf_s3_path = os.path.join(self._sheet.output,os.path.basename(pdf_path[0]))
        diff_ma = self.api.api("tool_lab.diff_ma")
        diff_ma.add_s3_result(
                                 main_id=self.option('main_id'),
                                 s3_output=pdf_s3_path,
                                 )
        self.set_output()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达量差异ma图",0],
            [r'diff_ma.pdf', 'pdf', '表达量差异ma图', 0],
            [r'diff_ma.png', 'png', '表达量差异ma图', 0],
            [r'diff_ma.svg', 'svg', '表达量差异ma图', 0],
        ])

        super(DiffMaWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.diff_ma import DiffMaWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "diff_volcano" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.diff_ma",
            "options": dict(
                raw_file='/mnt/ilustre/users/sanger-dev/workspace/20210322/DiffexpBatch_jgo0_k3guobgpi1lrvksgs4pu61_6989_3353/DiffexpBatch/output/ma.xls',
                pvalue=0.05,
                fc=2,
                x_axis_name="log10(TPM)",
                y_axis_name="log2(FC)",
                title_name="MA Plot",
                color="ref_blue_grey"
            )
        }
        wsheet = Sheet(data=data)
        wf =DiffMaWorkflow(wsheet)
        wf.sheet.id = 'diff_ma'
        wf.sheet.project_sn = 'diff_ma'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

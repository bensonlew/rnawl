# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'

import os
import gevent
import datetime
import unittest
import types
import shutil
from bson.objectid import ObjectId
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import HTMLParser
from mbio.packages.tool_lab.common import down_gbk_files


class AntismashWorkflow(Workflow):
    """
    antismash

    """
    def __init__(self, wsheet_object):
        self.samples = {}
        self._sheet = wsheet_object
        super(AntismashWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "input_file", "type": "infile", "format": "tool_lab.no_empty"},
            {"name": "input_format", "type": "string", "default": "gbk"},
            {"name": "taxon", "type": "string", "default": "bacteria"},  # 选择物种为真菌fungi时不能输入fasta文件
            {"name": "main_id", "type": "string"},
            {'name': "update_info", 'type': 'string'},
            {'name': 'project_data', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        if self.option('project_data'):
            self.project_data = eval(HTMLParser.HTMLParser().unescape(self._sheet.option('project_data')))
            for i in self.project_data['specimens']:
                self.samples[i['id']] = i['name']
            assemble_dir = os.path.join(self._sheet.work_dir, "gbk_dir")
            if os.path.exists(assemble_dir):
                shutil.rmtree(assemble_dir)
            (self.assemble_dir, self.analysis_type) = down_gbk_files(self.project_data['my_type'],
                                                                     self.project_data['db_version'], assemble_dir,
                                                                     self._sheet.option("project_task_id"),
                                                                     self.samples)

    def run(self):
        self.run_antismash()
        super(AntismashWorkflow, self).run()

    def run_antismash(self):
        self.antismash = self.add_tool("tool_lab.antismash")
        if self.option('project_data'):
            sample = self.samples.values()[0]
            gbk = os.path.join(self.work_dir, "gbk_dir", sample + ".gbk")
            self.logger.info(self._sheet.output)
            options = {
                "input_file": gbk,
                "input_format": self.option('input_format'),
                "taxon": self.option('taxon')
            }
        else:
            options = {
                "input_file": self.option('input_file'),
                "input_format": self.option('input_format'),
                "taxon": self.option('taxon')
            }
        self.antismash.set_options(options)
        self.antismash.on('end', self.set_output)
        self.antismash.run()

    """
    def set_output(self):
        if os.path.exists(self.antismash.output_dir + "/core_structures"):
            if not os.path.exists(self.output_dir + "/core_structures"):
                os.mkdir(self.output_dir + "/core_structures")
            for file in os.listdir(self.antismash.output_dir + "/core_structures"):
                os.link(os.path.join(self.antismash.output_dir, "core_structures", file), os.path.join(self.output_dir, "core_structures", file))
        os.link(os.path.join(self.antismash.output_dir, "antismash_anno.xls"), os.path.join(self.output_dir, "antismash_anno.xls"))
        os.link(os.path.join(self.antismash.output_dir, "gene_antismash.xls"), os.path.join(self.output_dir, "gene_antismash.xls"))
        os.link(os.path.join(self.antismash.output_dir, "gene_list.xls"), os.path.join(self.output_dir, "gene_list.xls"))
        self.set_db()
    """

    def set_db(self):
        remote_dir = self._sheet.output + '/' + 'core_structures'
        print remote_dir
        self.logger.info(self._sheet.output)
        #file_path = os.path.join(remote_dir, self.option('output'))
        self.logger.info("开始导表")
        api_antismash = self.api.api("tool_lab.antismash")
        api_antismash.add_Antismash(self.antismash.output_dir,remote_dir,main_id = self.option('main_id'))
        self.logger.info("导表结束")
        self.end()

    def download_file(self):
        """
        download file from s3
        :return:
        """
        assemble_dir = os.path.join(self.work_dir, "gbk_dir")
        for sample in self.samples.keys():
            sample_path = os.path.join(assemble_dir, self.samples[sample] + ".gbk")
            assemble_dir3 = os.path.join(assemble_dir, self.samples[sample])
            dir_list = os.listdir(assemble_dir3)
            for file2 in dir_list:
                os.system("cat {} >> {}".format(os.path.join(assemble_dir3, file2), sample_path))

    def set_output(self):
        if os.path.exists(self.output_dir + "/" + "antismash_anno.xls"):
            os.remove(self.output_dir + "/" + "antismash_anno.xls")
        if os.path.exists(self.output_dir + "/" + "gene_antismash.xls"):
            os.remove(self.output_dir + "/" + "gene_antismash.xls")
        if os.path.exists(self.output_dir + "/" + "gene_list.xls"):
            os.remove(self.output_dir + "/" + "gene_list.xls")
        if os.path.exists(self.antismash.output_dir + "/" + "antismash_anno.xls"):
            os.link(self.antismash.output_dir + "/"+ "antismash_anno.xls", self.output_dir + "/" + "antismash_anno.xls")
            os.link(self.antismash.output_dir + "/" + "gene_antismash.xls", self.output_dir + "/" + "gene_antismash.xls")
            os.link(self.antismash.output_dir + "/" + "gene_list.xls", self.output_dir + "/" + "gene_list.xls")
            allfiles = os.listdir(self.antismash.output_dir + "/" + "core_structures")
            #oldfiles = [os.path.join(self.antismash.output_dir + "/" + "core_structures", i) for i in allfiles]
            #newfiles = [os.path.join(self.output_dir + "/" + "core_structures", i) for i in allfiles]
            print self.antismash.output_dir + "/" + "core_structures"
            print self.output_dir + "/" + "core_structures"
            if not os.path.exists(self.output_dir + "/" + "core_structures"):
                os.mkdir(self.output_dir + "/" + "core_structures")
            for newfile in allfiles:
                if os.path.exists(os.path.join(self.output_dir + "/" + "core_structures", newfile)):
                    os.remove(newfile)
            for i in allfiles:
                os.link(os.path.join(self.antismash.output_dir + "/" + "core_structures", i), os.path.join(self.output_dir + "/" + "core_structures", i))
            self.set_db()
        else:
            api_antismash = self.api.api("tool_lab.antismash")
            api_antismash.add_Antismash2(main_id=self.option('main_id'))
            self.end()
    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AntismashWorkflow, self).end()
# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

import pandas as pd
import os
from biocluster.workflow import Workflow
import datetime
import types
from bson.objectid import ObjectId
import re
import shutil
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file
import HTMLParser
from mbio.packages.tool_lab.common import down_seq_files

class MlstTreeWorkflow(Workflow):
    """
    多个基因组的MLST分析，构建进化树
    """
    def __init__(self, wsheet_object):
        self.samples = {}
        self._sheet = wsheet_object
        super(MlstTreeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "tool_lab.fasta_dir"},  # 序列文件夹
            {"name": "species", "type": "string", },  # 物种名称
            {'name': 'task_id', 'type': 'string'},
            {'name': 'project_data', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        if self.option('project_data'):
            print(self._sheet)
            self.project_data = eval(HTMLParser.HTMLParser().unescape(self._sheet.option('project_data')))
            for i in self.project_data['specimens']:
                self.samples[i['id']] = i['name']
            assemble_dir = os.path.join(self._sheet.work_dir, "assemble_dir")
            if os.path.exists(assemble_dir):
                shutil.rmtree(assemble_dir)
            (self.assemble_dir, self.analysis_type) = down_seq_files(self.project_data['my_type'],
                                                                     self.project_data['db_version'], assemble_dir,
                                                                     self._sheet.option("project_task_id"),
                                                                     self.samples)

    def run(self):
        self.run_mlst()
        super(MlstTreeWorkflow, self).run()

    def run_mlst(self):
        dir = ''
        if self.option('project_data'):
            dir = self.download_file()
        elif self.option('fasta_dir').is_set:
            dir = self.option('fasta_dir')
        self.mlst = self.add_module("tool_lab.mlst")
        self.mlst.set_options({
            'fasta_dir': dir,
            'species': self.option('species')
        })
        self.mlst.on('end', self.set_db)
        self.mlst.run()

    def download_file(self):
        """
        download file from s3
        :return:
        """
        path = ''
        if self.analysis_type in ['complete']:
            assemble_dir2 = os.path.join(self.work_dir, "assemble_dir2")
            if os.path.exists(assemble_dir2):
                shutil.rmtree(assemble_dir2)
            os.mkdir(assemble_dir2)
            for sample in self.samples.keys():
                sample_path = os.path.join(assemble_dir2, self.samples[sample] + ".fna")
                assemble_dir3 = os.path.join(self.assemble_dir, self.samples[sample])
                dir_list = os.listdir(assemble_dir3)
                for file2 in dir_list:
                    os.system("cat {} >> {}".format(os.path.join(assemble_dir3, file2), sample_path))
            with open(assemble_dir2 + "/list.txt", "w") as f:
                f.write("sample\tfile\n")
                for sample in self.samples.keys():
                    f.write("{}\t{}\n".format(self.samples[sample], self.samples[sample] + ".fna"))
            path = assemble_dir2
        elif self.analysis_type in ['uncomplete']:
            with open(self.assemble_dir + "/list.txt", "w") as f:
                f.write("sample\tfile\n")
                for sample in self.samples.keys():
                    f.write("{}\t{}\n".format(self.samples[sample], self.samples[sample] + ".fna"))
            path = self.assemble_dir
        return path

    def set_db(self):

        """
        保存结果标准化数据到mongo数据库中
        """
        link_dir(self.mlst.output_dir, self.output_dir)
        mlst_tree = self.api.api('tool_lab.mlst_tree')
        path1 = self.output_dir + "/all.ST.xls"
        path2 = self.output_dir + "/mlst_tree.nwk"
        if os.path.exists(path2):
            mlst_tree.add_mlsttree_detail(ObjectId(self.option("main_id")), path1, file2=path2)
        else:
            mlst_tree.add_mlsttree_detail(ObjectId(self.option("main_id")), path1)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.mlst.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Mlst_Tree结果输出目录"],
            ["./mlst_tree.nwk", "xls", "MLST的进化树文件"],
            ["./all.ST.xls", "xls", "MLST的分型统计表"],
        ])
        result_dir.add_regexp_rules([
            ["", "", "Mlst_Tree结果输出目录"],
            ["./mlst_tree.nwk", "nwk", "MLST的进化树文件"],
            ["./all.ST.xls", "xls", "MLST的分型统计表"],
        ])
        super(MlstTreeWorkflow, self).end()

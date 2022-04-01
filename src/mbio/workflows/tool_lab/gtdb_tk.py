# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
from biocluster.workflow import Workflow
import datetime
from biocluster.core.exceptions import OptionError
import shutil
from bson.objectid import ObjectId
from mbio.packages.metagbin.common_function import link_dir,link_file
import HTMLParser
from mbio.packages.tool_lab.common import down_seq_files
from biocluster.file import download,exists


class GtdbTkWorkflow(Workflow):
    """
    GTDB-TK的物种注释分类
    """
    def __init__(self, wsheet_object):
        self.samples = {}
        self._sheet = wsheet_object
        if self._sheet.option('project_data'):
            print(self._sheet)
            self.project_data = eval(HTMLParser.HTMLParser().unescape(self._sheet.option('project_data')))
            if self.project_data['my_type'] in ['bacgenome']:
                for i in self.project_data['specimens']:
                    self.samples[i['id']] = i['name']
                assemble_dir = os.path.join(self._sheet.work_dir, "assemble_dir")
                if os.path.exists(assemble_dir):
                    shutil.rmtree(assemble_dir)
                (self.assemble_dir, self.analysis_type) = down_seq_files(self.project_data['my_type'],
                                                                         self.project_data['db_version'], assemble_dir,
                                                                         self._sheet.option("project_task_id"),
                                                                         self.samples)
            elif self.project_data['my_type'] in ['metagbin']:
                assemble_dir = os.path.join(self._sheet.work_dir, "assemble_dir")
                if os.path.exists(assemble_dir):
                    shutil.rmtree(assemble_dir)
                os.mkdir(assemble_dir)
                download(self.project_data['genome_path'], assemble_dir+"/"+self.project_data['genome_id']+".fna")
                if os.path.exists(assemble_dir+"/"+self.project_data['genome_id']+".fna"):
                    pass
                else:
                    raise OptionError("下载序列失败或有问题！")

        super(GtdbTkWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {'name': "update_info", 'type': 'string'},
            {'name': "input_format", 'type': 'string'},
            {'name': 'input_dir', 'type': 'infile', 'format': 'toolapps.fasta_dir'},
            {'name': 'input_file', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'project_data', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.gtdb_tk = self.add_tool("toolapps.gtdb_tk")

    def run(self):
        self.run_gtdbtk()
        super(GtdbTkWorkflow, self).run()

    def get_input(self):
        dir = ''
        if self.option('project_data'):
            if self.project_data['my_type'] in ['metagbin']:
                dir = self.work_dir+"/assemble_dir"
            elif self.project_data['my_type'] in ['bacgenome']:
                dir = self.download_file()
        else:
            if os.path.exists(self.work_dir + "/input"):
                shutil.rmtree(self.work_dir + "/input")
            os.mkdir(self.work_dir + "/input")
            if self.option("input_file").is_set:
                name = ".".join(os.path.basename(self.option("input_file").prop['path']).split(".")[0:-1])
                link_file(self.option("input_file").prop['path'], self.work_dir + "/input/" + name + ".fna")
            elif self.option("input_dir").is_set:
                with open(self.option("input_dir").prop['path'] + "/list.txt", "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        lin = line.strip().split("\t")
                        os.link(lin[0], self.work_dir + "/input/" + lin[1] + ".fna")
            else:
                raise OptionError("必须输入input_file文件! 或 input_dir文件夹！")
            dir = self.work_dir+"/input"
        return dir

    def run_gtdbtk(self):

        dir = self.get_input()
        options = {
            "genome_dir": dir,
        }
        self.gtdb_tk.set_options(options)
        self.gtdb_tk.on("end", self.set_output)
        self.gtdb_tk.run()

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
            path = assemble_dir2
        elif self.analysis_type in ['uncomplete']:
            path = self.assemble_dir
        return path

    def set_output(self):
        link_dir(self.gtdb_tk.output_dir, self.output_dir)
        self.set_db()

    def set_db(self):
        self.logger.info("导表开始")
        file_taxon = self.output_dir + "/all.taxon.xls"
        api_gtdbtk = self.api.api("tool_lab.gtdb_tk")
        api_gtdbtk.add_gtdb_detail(file_taxon, ObjectId(self.option("main_id")))
        self.end()
        self.logger.info("导表结束")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "GTDB结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(GtdbTkWorkflow, self).end()
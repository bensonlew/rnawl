# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.tool_lab.common_function import fasta_dir_sort
from mbio.packages.tool_lab.common import down_seq_files
import os
import re, shutil
import HTMLParser
import gevent

class QuastWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self.samples = {}
        self._sheet = wsheet_object
        super(QuastWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "file_style", "type": "string"},
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},
            {'name': 'project_data', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.tools = []
        if self.option('project_data'):
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
        self.run_quast()
        super(QuastWorkflow, self).run()

    def run_quast(self):
        """
        QUAST基因组评估
        :return:
        """
        path2 =''
        if os.path.exists(os.path.join(self.work_dir, "fasta_dir_sort")):
            shutil.rmtree(os.path.join(self.work_dir, "fasta_dir_sort"))
        os.mkdir(os.path.join(self.work_dir, "fasta_dir_sort"))
        if self.option("fasta_dir").is_set:
            fasta_dir_sort(self.option("fasta_dir").prop['path'], os.path.join(self.work_dir, "fasta_dir_sort"))
            path2 = os.path.join(self.work_dir,"fasta_dir_sort")
        elif self.option("fasta").is_set:
            name = self.option("fasta").prop["path"].split("/")[-1]
            if name.endswith(".gz"):
                os.system("gzip -d {}".format(self.option("fasta").prop["path"]))
            os.link(self.option("fasta").prop["path"].strip(".gz"),self.work_dir+"/"+"fasta_dir_sort"+"/"+name.split(".")[0] + ".fasta")
            path2 =os.path.join(self.work_dir,"fasta_dir_sort")
        elif self.option("project_data"):
            path2 =  self.download_file()
        self.quast = self.add_tool("tool_lab.quast")
        opts = {
                "fasta_dir": path2,
                }
        if self.option("ref_fasta").is_set:
            opts["ref_fasta"] =self.option("ref_fasta")
        self.quast.set_options(opts)
        self.quast.on('end', self.set_output)
        self.quast.run()

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
        link_dir(self.quast.output_dir,self.output_dir)
        self.set_db()

    def set_db(self):
        api_main = self.api.api("tool_lab.quast")
        if self.option('main_id'):
            main_id =  self.option('main_id')
        else:
            main_id = api_main.add_quast_main()
        picture = self._sheet.output + '/' + 'cumulative_plot.pdf'
        self.logger.info(picture)
        if self.option("ref_fasta").is_set:
            ref = True
        else:
            ref = False
        api_main.add_quast_detail(main_id,file_path=self.output_dir + "/QUAST_Evaluation.xls",picture_path=picture,is_ref=ref)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "样品流程分析结果目录", 0],
            ["./cumulative_plot.pdf", "pdf", "基因组分布图", 0],
            ["./QUAST_Evaluation.xls", "xls", "QUAST评估结果表", 0],
        ])
        super(QuastWorkflow, self).end()
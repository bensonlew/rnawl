# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.tool_lab.common_function import fasta_dir_sort
import os
import re, shutil
import subprocess
import gevent
from mbio.packages.tool_lab.common import down_seq_files
import HTMLParser

class TagseqWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self.samples = {}
        self._sheet = wsheet_object
        super(TagseqWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "file_style", "type": "string"},
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},
            {"name": "ref_style", "type": "string"},
            {"name": "ref", "type": "string"}, # 目的片段
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

    def check_options(self):
        if self.option("project_data"):
            pass
        else:
            if (not self.option("fasta").is_set) and (not self.option("fasta_dir").is_set):
                raise OptionError("请传入序列文件或序列文件夹！", code="")
            if (not self.option("ref")) and (not self.option("ref_fasta").is_set):
                raise OptionError("请传入目的片段！", code="")

    def run(self):
        if self.option("fasta_dir").is_set:
            fasta_dir_sort(self.option("fasta_dir").prop['path'],self.work_dir + "/fasta_dir_sort")
        elif self.option("project_data"):
            self.download_file()
        self.run_tagseq()
        super(TagseqWorkflow, self).run()

    def run_tagseq(self):
        """
        目的片段查找
        :return:
        """
        if self.option("ref"):
            seq = self.option("ref")
            with open(self.work_dir + "/query.fasta","w") as t:
                if seq.startswith(">"):
                    t.write(">" + self.option("ref").split("\n")[0].lstrip(">") + "\n")
                    t.write("\n".join(self.option("ref").split("\n")[1:]) + "\n")
                elif seq.startswith("&gt;"):
                    t.write(">" + self.option("ref").split("\n")[0].lstrip("&gt;") + "\n")
                    t.write("\n".join(self.option("ref").split("\n")[1:]) + "\n")
                else:
                    t.write(">seq1" + "\n")
                    t.write("\n".join(self.option("ref").split("\n")) + "\n")
            query = self.work_dir + "/query.fasta"
        else:
            query = self.option("ref_fasta")

        if self.option("fasta").is_set:
            self.blast = self.add_tool("align.blast")
            opts = {
                "query": query,
                "blast": "blastn",
                "query_type": "nucl",
                "reference_type": "nucl",
                "reference": self.option("fasta"),
                "database": "customer_mode",
                "outfmt" : 6,
                "num_alignment": 1,
            }
            self.blast.set_options(opts)
            self.blast.on('end', self.set_db)
            self.blast.run()
        elif self.option("project_data"):
            self.blast = self.add_tool("align.blast")
            sample = self.samples.values()[0]
            fasta = os.path.join(self.work_dir, "assemble_dir", sample + ".fna")
            opts = {
                "query": query,
                "blast": "blastn",
                "query_type": "nucl",
                "reference_type": "nucl",
                "reference": fasta,
                "database": "customer_mode",
                "outfmt": 6,
                "num_alignment": 1,
            }
            self.blast.set_options(opts)
            self.blast.on('end', self.set_db)
            self.blast.run()
        else:
            files = os.listdir(os.path.join(self.work_dir, "fasta_dir_sort"))
            for file in files:
                self.blast = self.add_tool("align.blast")
                opts = {
                    "query": query,
                    "blast": "blastn",
                    "query_type": "nucl",
                    "reference_type": "nucl",
                    "reference": os.path.join(self.work_dir, "fasta_dir_sort") + "/" + file,
                    "database": "customer_mode",
                    "outfmt": 6,
                    "num_alignment": 1
                }
                self.blast.set_options(opts)
                self.tools.append(self.blast)
                self.logger.info(self.tools)
                if self.tools:
                    if len(self.tools) > 1:
                        self.on_rely(self.tools, self.set_db)
                    elif len(self.tools) == 1:
                        self.tools[0].on('end', self.set_db)
                    for tool in self.tools:
                        self.logger.info(tool)
                        gevent.sleep(1)
                        tool.run()

    def download_file(self):
        """
        download file from s3
        :return:
        """
        if self.analysis_type in ['complete']:
            assemble_dir = os.path.join(self.work_dir, "assemble_dir")
            for sample in self.samples.keys():
                sample_path = os.path.join(assemble_dir, self.samples[sample] + ".fna")
                assemble_dir3 = os.path.join(assemble_dir, self.samples[sample])
                dir_list = os.listdir(assemble_dir3)
                for file2 in dir_list:
                    os.system("cat {} >> {}".format(os.path.join(assemble_dir3, file2), sample_path))

    def set_db(self):
        api_main = self.api.api("tool_lab.tagseq")
        if self.option('main_id'):
            main_id =  self.option('main_id')
        else:
            main_id = api_main.add_tagseq_main()
        if self.option("fasta").is_set or self.option("project_data"):
            api_main.add_tagseq_detail(main_id=main_id, file_path=self.blast.output_dir)
            self.logger.info(self.blast.output_dir)
        else:
            for tool in self.tools:
                if os.listdir(tool.output_dir):
                    api_main.add_tagseq_detail(main_id=main_id, file_path=tool.output_dir)
        api_main.add_tagseq_main(main_id)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(TagseqWorkflow, self).end()
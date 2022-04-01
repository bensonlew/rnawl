# -*- coding: utf-8 -*-
# __author__: binbin.zhao
# modified: 20200622

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.tool_lab.common import down_seq_files
import HTMLParser
import shutil

class SsrAnalysisWorkflow(Workflow):
    """
    交互分析：所有样本的SSR设计
    """
    def __init__(self, wsheet_object):
        self.samples = {}
        self._sheet = wsheet_object
        if self._sheet.option('project_data'):
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
        super(SsrAnalysisWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "reffa", "type": "infile", "format": "sequence.fasta"},  # 参考基因组ref.fa
            {"name": "rept", "type": "string", "default": "10,6,5,5,5,5"},
            # {"name": "rept_1", "type": "int", "default": 10},
            # {"name": "rept_2", "type": "int", "default": 6},
            # {"name": "rept_3", "type": "int", "default": 5},
            # {"name": "rept_4", "type": "int", "default": 5},
            # {"name": "rept_5", "type": "int", "default": 5},
            # {"name": "rept_6", "type": "int", "default": 5},
            {"name": "ssr_distance", "type": "int", "default": 100},
            {"name": "update_info", "type": "string"},
            {'name': 'project_data', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("reffa"):
            raise OptionError("请设置参考基因组的ref.fa")

    def run_ssr_ref_primer(self):
        rept_list = self.option("rept").split(",")
        if rept_list[0] != "":
            rept_1 = int(rept_list[0])
        else:
            rept_1 = 10
        if rept_list[1] != "":
            rept_2 = int(rept_list[1])
        else:
            rept_2 = 6
        if rept_list[2] != "":
            rept_3 = int(rept_list[2])
        else:
            rept_3 = 5
        if rept_list[3] != "":
            rept_4 = int(rept_list[3])
        else:
            rept_4 = 5
        if rept_list[4] != "":
            rept_5 = int(rept_list[4])
        else:
            rept_5 = 5
        if rept_list[5] != "":
            rept_6 = int(rept_list[5])
        else:
            rept_6 = 5
        options = {
            "needini": True,
            "needprimer": False,
            "rept_1": rept_1,
            "rept_2": rept_2,
            "rept_3": rept_3,
            "rept_4": rept_4,
            "rept_5": rept_5,
            "rept_6": rept_6,
            # "rept_1": self.option("rept_1"),
            # "rept_2": self.option("rept_2"),
            # "rept_3": self.option("rept_3"),
            # "rept_4": self.option("rept_4"),
            # "rept_5": self.option("rept_5"),
            # "rept_6": self.option("rept_6"),
            "ssr_distance": self.option("ssr_distance"),
        }
        if self.option("project_data"):
            sample = self.samples.values()[0]
            path = self.download_file()
            fasta_file = os.path.join(self.work_dir, path, sample + ".fna")
            options["reffa"] = fasta_file
        else:
            options["reffa"] = self.option("reffa")
        self.ssr_ref_primer = self.add_module("tool_lab.ssr_ref_primer")
        self.ssr_ref_primer.set_options(options)
        self.ssr_ref_primer.on("end", self.set_output, "ssr_merge")
        self.ssr_ref_primer.run()

    def set_output(self, event):
        for f in os.listdir(event["bind_object"].output_dir):
            new = os.path.join(self.output_dir, f)
            if os.path.exists(new):
                os.remove(new)
            os.link(os.path.join(event["bind_object"].output_dir, f), new)
        self.set_db()

    def set_db(self):
        self.logger.info("开始进行导表")
        if self.option("main_id"):
            ssr_api = self.api.api("tool_lab.ssr_specimen")
            ssr_id = self.option("main_id")
            if self.option("project_data"):
                sample = self.samples.values()[0]
                stat_file = self.output_dir + "/" + sample + ".ssr.stat"
            else:
                stat_file = self.output_dir + "/ref.ssr.stat"
            ssr_api.add_sg_ssr_marker_stat(self.option("task_id"), ssr_id, stat_file)
        self.end()

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

    def run(self):
        self.run_ssr_ref_primer()
        super(SsrAnalysisWorkflow, self).run()

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SsrAnalysisWorkflow, self).end()

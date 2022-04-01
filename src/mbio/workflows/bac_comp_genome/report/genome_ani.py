# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modifies 20191104

import os,shutil
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
from mbio.packages.bac_comp_genome.common_function import get_fasta, get_sample_from_tree,link_dir
from biocluster.api.file.lib.transfer import MultiFileTransfer

class GenomeAniWorkflow(Workflow):
    """
    基因组ani计算
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenomeAniWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "sample_list", "type": "string"},
            {"name": "seq_dir", "type": "string"},  # 输入参考序列文件夹
            {"name": "method", "type": "string", "default": "ANIm"}, #ANIm,ANIb,ANIblastall,TETRA
            {"name": "linkage", "type": "string", "default": "average"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.ani = self.add_module("bac_comp_genome.genome_ani")

    def run_ani(self):
        transfer = MultiFileTransfer()
        transfer.add_download(self.option("seq_dir"), self.work_dir + "/seq_dir/")
        transfer.perform()
        if os.path.exists(self.work_dir + "/data"):
            shutil.rmtree(self.work_dir + "/data")
        os.mkdir(self.work_dir + "/data")
        for i in self.option("sample_list").split(","):
            os.link(self.work_dir + "/seq_dir/" + i + ".fna", self.work_dir + "/data/" + i + ".fna")
        self.ani.set_options({
            "seq_dir": self.work_dir + "/data",
            "method":self.option("method"),
            'linkage': self.option("linkage")
        })
        self.ani.on("end", self.set_output)
        self.ani.run()

    def run(self):
        self.run_ani()
        super(GenomeAniWorkflow, self).run()

    def set_output(self):
        link_dir(self.ani.output_dir, self.output_dir)
        self.set_db()

    def change_samples(self):
        samples_dict = {}
        files =os.listdir(self.work_dir + "/data/")
        for file in files:
            if file.endswith(".fna"):
                name = file.split('.fna')[0]
                new_name = name.replace(".", "_")
                samples_dict[name] = new_name
        return samples_dict

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.samples_dict = self.change_samples()
        self.logger.info(self.samples_dict)
        api_path = self.api.api("bac_comp_genome.common_api")
        main_id = self.option("main_id")
        samples = get_sample_from_tree(self.output_dir + "/hcluster.tre")
        api_path.add_main_detail2(self.output_dir + '/all.ANI_summary.xls', 'ani_detail', main_id, "specimen_id", main_name='ani_id', names_convert=self.samples_dict, names_co="true")
        api_path.add_main_tree2(self.output_dir + "/hcluster.tre", "ani", main_id, update_dic={"samples_order": samples, "status": "end", "main_id": ObjectId(main_id)})
        self.end()

    def end(self):
        repaths = [
            [".", "", "",0],
        ]
        regexps = [
            [r'.*\.xls$', 'xls', '',0],
            [r'.*\.stat$', 'stat', '',0]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(GenomeAniWorkflow, self).end()
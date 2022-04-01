# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_file,link_dir
import os

class SangerIdentityWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SangerIdentityWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "file_dir", "type": "infile", 'format': "denovo_rna_v2.common_dir"},
            {"name": "method", "type": "string","default": "false"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.identity = self.add_tool('tool_lab.sanger_identity')

    def check_options(self):
        if not self.option("file_dir").is_set:
            raise OptionError("必须输入文件！")

    def run(self):
        self.run_identity()
        super(SangerIdentityWorkflow, self).run()

    def run_identity(self):
        self.identity.set_options({
            "file_dir": self.option("file_dir"),
            "method": self.option("method")
        })
        self.identity.on("end", self.set_output)
        self.identity.run()

    def set_output(self):
        self.logger.info("set_output")
        stat_file = self.output_dir + "/sangerSeq_stat.xls"
        nt_file = self.output_dir + "/sangerSeq_NTblast.xls"
        if os.path.exists(stat_file):
            os.remove(stat_file)
        if os.path.exists(nt_file):
            os.remove(nt_file)
        os.link(self.identity.output_dir + "/sangerSeq_stat.xls",stat_file)
        if os.path.exists(self.identity.output_dir+ "/sangerSeq_NTblast.xls"):
            os.link(self.identity.output_dir+ "/sangerSeq_NTblast.xls",nt_file)
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_main = self.api.api("tool_lab.sanger_identity")
        if os.path.exists(self.output_dir + "/sangerSeq_NTblast.xls"):
            self.logger.info("有")
            api_main.add_main(file1 = self.output_dir + "/sangerSeq_stat.xls", file2=self.output_dir + "/sangerSeq_NTblast.xls",method =self.option("method"), main_id = self.option('main_id'))
        else:
            api_main.add_main(file1=self.output_dir + "/sangerSeq_stat.xls", main_id=self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "一代测序菌种鉴定结果目录", 0],
            ["./sangerSeq_stat.xls", "xls", "一代测序菌鉴概览", 0],
            ["./sangerSeq_NTblast.xls", "xls", "一代测序菌鉴NT比对详情表", 0],
        ])
        super(SangerIdentityWorkflow, self).end()

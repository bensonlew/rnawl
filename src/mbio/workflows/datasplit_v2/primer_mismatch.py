# -*- coding: utf-8 -*-
# __author__: zengjing
# last_modify: 20210223

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class PrimerMismatchWorkflow(Workflow):
    """
    primer错配：每个barcode对应的前n的primer
    """
    def __init__(self, wsheet_object):
        super(PrimerMismatchWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "lib_list", "type": "infile", "format": "datasplit.path"},  # 文库extendedFrags.fastq和对应的barcode文件
            {"name": "top_row", "type": "int", "default": 10},  # 取前多少的引物
            {"name": "update_info", "type": "string"},
            {"name": "verify_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("lib_list").is_set:
            raise OptionError("必须设置lib_list")

    def run_primer_mismatch(self):
        self.start_times, self.end_times = 0, 0
        with open(self.option("lib_list").prop["path"], "rb") as f:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                options = {
                    "fq": item[0],
                    "barcode_primer_info": item[1],
                    "top_row": self.option("top_row")
                }
                self.primer_mismatch = self.add_tool("datasplit_v2.primer_mismatch")
                self.primer_mismatch.set_options(options)
                self.primer_mismatch.on("end", self.set_output)
                self.primer_mismatch.run()
                self.start_times += 1

    def set_output(self, event):
        self.logger.info("设置结果目录")
        obj = event["bind_object"]
        info = os.path.join(obj.output_dir, "primer_reads.txt")
        with open(info, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                self.info_list.append(line)
        self.end_times += 1
        if self.start_times == self.end_times:
            all_primer = os.path.join(self.output_dir, "primer_reads.txt")
            with open(all_primer, "wb") as w:
                w.write("Sample\tPrimer\tReads\n")
                for line in self.info_list:
                    w.write(line)
            self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        if self.option("verify_id"):
            datasplit_api = self.api.api("datasplit.datasplit_new")
            verify_id = self.option("verify_id")
            info_path = os.path.join(self.output_dir, "primer_reads.txt")
            datasplit_api.add_sg_meta_verify_primer_detail(verify_id, info_path)
        self.end()

    def run(self):
        self.info_list = []
        self.run_primer_mismatch()
        super(PrimerMismatchWorkflow, self).run()

    def end(self):
        super(PrimerMismatchWorkflow, self).end()

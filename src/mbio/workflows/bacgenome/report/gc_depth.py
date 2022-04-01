# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# __modify__ = '2020/03/31'

from biocluster.workflow import Workflow
from mbio.packages.bacgenome.common import fasta_cutoff
import os


class GcDepthWorkflow(Workflow):
    """
    用于扫描图的gc_depth分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GcDepthWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},
            {"name": "samp", "type": "string"},
            {"name": "min_len", "type": "string", "default":"1000"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.gc_depth = self.add_tool('bacgenome.gc_depth')
        self.file_path = self._sheet.output

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.get_info()
        self.run_gcdepth()
        super(GcDepthWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_path = self.api.api('bac_assem.assess')
        self.logger.info("开始进行导表")
        gc_path = self.gc_depth.output_dir
        gc_remote = self.file_path
        api_path.add_assess_gc("gc_depth_detail", self.option("main_id"), self.option("samp"), gc_path, gc_remote)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.gc_depth.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果目录"],
        ])
        super(GcDepthWorkflow, self).end()

    def run_gcdepth(self):
        fasta_cutoff(self.option('fa').prop['path'], self.option("min_len"), self.work_dir + "/" + "seq.fa")
        opts = {
            "seq": self.work_dir + "/" + "seq.fa",
            "fastq_list": self.path,
        }
        self.gc_depth.set_options(opts)
        self.gc_depth.on("end", self.set_db)
        self.gc_depth.run()

    def get_info(self):
        self.path = os.path.join(self.work_dir, 'list.txt')
        with open(self.path,'w') as file:
            name = self.option("samp")
            lib =str(self.option("read1").prop['path'].split("/")[-1].split(".clean")[0].split("_PE_")[1])
            file.write("{}\t{}\t{}\n".format(name, "PE" + lib, self.option("read1").prop['path']+";"+self.option("read2").prop['path']))

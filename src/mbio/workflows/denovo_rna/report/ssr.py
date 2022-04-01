# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from biocluster.workflow import Workflow
from mbio.api.to_file.denovo import *


class SsrWorkflow(Workflow):
    """
    报告中计算比对质量评估时使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SsrWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gene_fasta", "type": "string"},  # gene转录本，fasta格式
            {"name": "bed_file", "type": "string"},  # bed格式文件
            {"name": "update_info", "type": "string"},
            {"name": "insert_id", "type": "string"},  # 插入的主表的ID
            {"name": "primer", "type": "string"}  # 是否设置引物
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ssr = self.add_tool('denovo_rna.gene_structure.ssr')

    def run(self):
        self.ssr.set_options({
            "fasta": self.option("gene_fasta"),
            "bed": self.option("bed_file"),
            "primer": True
        })
        print(self.option("bed_file"))
        self.ssr.on("end", self.set_db)
        self.ssr.run()
        super(SsrWorkflow, self).run()

    def set_db(self):
        fasta_file_name = os.path.basename(self.option("gene_fasta"))
        api_ssr = self.api.denovo_gene_structure
        ssr_path = self.ssr.output_dir + "/" + fasta_file_name + ".misa"
        stat_path = self.ssr.work_dir + "/" + fasta_file_name + ".statistics"
        if not os.path.isfile(ssr_path):
            raise Exception("找不到报告文件:{}".format(ssr_path))
        api_ssr.add_ssr_detail(ssr_path, self.option("insert_id"))
        api_ssr.add_ssr_stat(stat_path, self.option("insert_id"))
        if self.option("primer") == "true":
            print("lllllllllllllll")
            primer_path = self.ssr.work_dir + "/" + fasta_file_name + ".misa.results"
            print(primer_path)
            api_ssr.add_ssr_primer(primer_path, self.option("insert_id"))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.ssr.output_dir)
        result_dir.add_regexp_rules([
            [r"misa$", "misa", "ssr结果"]
        ])
        super(SsrWorkflow, self).end()

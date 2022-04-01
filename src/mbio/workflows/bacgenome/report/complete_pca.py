# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# __modify__ = '2020/04/02'

from biocluster.workflow import Workflow
import os


class CompletePcaWorkflow(Workflow):
    """
    完成图的scaffold kmer交互分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CompletePcaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "genome_id", "type": "string"},
            {"name": "samp", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.pca = self.add_tool('bacgenome.kmer_pca')

    def run(self):
        self.run_pca()
        super(CompletePcaWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_path = self.api.api('bac_assem.blast')
        self.logger.info("开始进行导表")
        table = os.path.join(self.pca.output_dir, "4mer_pca_sites.xls")
        table2 = os.path.join(self.pca.output_dir, "4mer_pca_importance.xls")
        api_path.add_kmer_pca_detail("complete_kmer_pca_detail", self.option("main_id"), self.option("samp"), table,
                                     table2, self.option("genome_id"))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.pca.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果目录"],
        ])
        super(CompletePcaWorkflow, self).end()

    def run_pca(self):
        opts = {
            "scaf": self.option("fa")
        }
        self.pca.set_options(opts)
        self.pca.on("end", self.set_db)
        self.pca.run()

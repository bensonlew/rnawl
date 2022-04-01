# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# __modify__ = '20190216'

from biocluster.workflow import Workflow
import os
import json

class TrfPcoaWorkflow(Workflow):
    """
    四核苷酸频率pcoa
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TrfPcoaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "bin_id", "type": "string"},
            {'name': 'input_genome', 'type': 'infile', 'format': 'sequence.fasta'},  # bin的序列
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.pcoa = self.add_module("metagbin.tre_pcoa")

    def run(self):
        self.run_seq()
        super(TrfPcoaWorkflow, self).run()

    def run_seq(self):
        opts =({
            "input_genome":self.option("input_genome"),
        })
        self.pcoa.set_options(opts)
        self.pcoa.on('end',self.set_db)
        self.pcoa.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        main_id =self.option('main_id')
        self.api_path = self.api.api('metagbin.tetra_pcoa')
        self.api_path.add_pcoa_detail(main_id, self.pcoa.output_dir + '/pcoa_sites.xls',self.pcoa.output_dir + '/pcoa_eigenvaluespre.xls')
        self.end()

    def end(self):
        bin_name = self.option("bin_id")
        result_dir = self.add_upload_dir(self.pcoa.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "%s四核苷酸频率pcoa结果目录" % bin_name],
            ["pcoa_eigenvalues.xls", "", "%s四核苷酸频率的pcoa矩阵特征值" % bin_name],
            ["tetra.summary.xls  ", "", "%s四核苷酸频率结果表" % bin_name],
            ["pcoa_eigenvaluespre.xls", "", "%s四核苷酸频率pcoa特征解释度百分比" % bin_name],
            ["pcoa_sites.xls", "", "%s四核苷酸频率scafollds坐标表" % bin_name],
            ["bray_curtis_tetra.pcoa.xls", "", "%s基于四核苷酸计算的bray_curtis距离的矩阵" % bin_name],
        ])
        super(TrfPcoaWorkflow, self).end()
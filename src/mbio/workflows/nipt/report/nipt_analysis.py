# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# create_time：20170518
from biocluster.workflow import Workflow
import os


class NiptAnalysisWorkflow(Workflow):
    """
    该workflow用于nipt的交互分析，实现对bed文件进行处理，获得百分比数据，用于网页端的展示
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(NiptAnalysisWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "bed_file", "type": "string"},  # 输入样本的id
            {"name": "bw", "type": "int", "default": 10},
            {"name": "bs", "type": "int", "default": 1},
            {"name": "ref_group", "type": "int", "default": 2},
            {"name": "update_info", "type": "string"},
            {"name": "nipt_task_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.niptanalysis = self.add_tool('nipt.bed_analysis')


    def run_niptanalysis(self):
        print "test", self.option("bed_file")
        # self.api.nipt_analysis.add_bed_file(self.option("bed_file").prop['path'])
        bed_file = self.api.nipt_analysis.export_bed_file(sample=self.option("bed_file"), dir=self.work_dir)
        options = {
            'bed_file': bed_file,
            'bw': self.option('bw'),
            'bs': self.option('bs'),
            'ref_group': self.option('ref_group'),
            'single_chr': "false"
        }
        self.niptanalysis.set_options(options)
        self.niptanalysis.on('end', self.set_db)
        self.output_dir = self.niptanalysis.output_dir
        self.niptanalysis.run()


    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "nipt分析结果输出目录"],
            ["./z.xls", "xls", "百分比统计数据（Z值）"],
            ["./zz.xls", "xls", "各个染色体上不同bin的统计数值"]
        ])
        super(NiptAnalysisWorkflow, self).end()

    def set_db(self):
        """
        报存分析结果到mongo数据库中
        """
        api_niptanalysis = self.api.nipt_analysis
        file_name = str(self.option('bed_file')) + '_' + str(self.option("bw")) + "_" + str(self.option('bs'))
        z_file_path = self.output_dir + '/' + file_name + '_z.xls'
        zz_file_path = self.output_dir + '/' + file_name + '_zz.xls'
        if not os.path.isfile(z_file_path):
            raise Exception("找不到报告文件:{}".format(z_file_path))
        if not os.path.isfile(zz_file_path):
            raise Exception("找不到报告文件:{}".format(zz_file_path))
        print 'start insert'
        api_niptanalysis.add_z_result(file_path=z_file_path, table_id=self.option("nipt_task_id"))
        api_niptanalysis.add_zz_result(file_path=zz_file_path, table_id=self.option("nipt_task_id"))
        print 'end insert'
        self.end()

    def run(self):
        self.run_niptanalysis()
        super(NiptAnalysisWorkflow, self).run()

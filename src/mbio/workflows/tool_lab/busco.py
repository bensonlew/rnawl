# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
import pandas as pd
import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import matplotlib.pyplot as plt
from mbio.packages.ref_rna_v2.chart import Chart
import glob



class BuscoWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BuscoWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'fa', 'type': 'infile', 'format': 'denovo_rna_v2.trinity_fasta'},
            {'name': 'classification', 'type': 'string'},
            {'name': 'fa_dir', 'type': 'infile', 'format': 'denovo_rna_v2.common_dir'},
            {'name': 'odb9', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.run_busco()
        super(BuscoWorkflow, self).run()

    def run_busco(self):
        self.busco = self.add_module("tool_lab.busco")
        self.busco.set_options({
            'fa_dir': self.option('fa_dir'),
            'odb9': self.option('odb9')
        })
        self.busco.on('end', self.set_db)
        self.busco.run()


    # def run_busco(self):
    #     self.busco = self.add_tool("tool_lab.busco")
    #     self.busco.set_options({
    #         'fa': self.option('fa'),
    #         'odb9': self.option('odb9')
    #     })
    #     self.busco.on('end', self.set_db)
    #     self.busco.run()

    # def run_picture_bar(self):
    #     # sample_name = os.path.basename(self.option('fa').path).split('.')[0]
    #     # summary_result = self.busco.option('summary_result').path
    #     # busco = pd.read_table(summary_result, header=0, index_col=0, sep='\t')
    #     # persent = busco['persent'].tolist()
    #     # persent = persent[1:-1]
    #     # color_list = ['lightblue', 'deepskyblue', 'yellow', 'red']
    #     # label_list = busco.index.tolist()[1:-1]
    #     # for i in range(len(persent)):
    #     #     plt.bar(sample_name, persent[i], width=0.2, bottom=sum(persent[:i]), color=color_list[i % len(color_list)],
    #     #             alpha=0.5, label=label_list[i])
    #     # plt.legend( bbox_to_anchor=(1, 1), loc='upper center')
    #     # plt.savefig(os.path.join(self.busco.output_dir, 'busco.jpg'))
    #     # self.set_db()
    #     with open(os.path.join(self.busco.output_dir, 'short_summary_busco_result.txt'), 'r') as s:
    #         for line in s.readlines():
    #             if line.startswith('#') or line.startswith('\n'):
    #                 continue
    #             else:
    #                 if line.strip().startswith('C:'):
    #                     text = line.strip()
    #     sample_name = os.path.basename(self.option('fa').path).split('.')[0]
    #     self.bar = self.add_tool("tool_lab.busco_bar")
    #     self.bar.set_options({
    #         'summary_result': self.busco.option('summary_result').path,
    #         'title': text,
    #         'xlab': sample_name
    #     })
    #     self.bar.on('end', self.set_db)
    #     self.bar.run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        print self.busco_text_dict
        print self.busco_summary_dict
        chart.chart_busco(self.busco_summary_dict, self.busco_text_dict)
        chart.to_pdf()

        # move pdf to result dir
        pdf_files = glob.glob(self.work_dir + "/*.pdf")
        for file in pdf_files:
            os.link(file, self.output_dir + "/busco.pdf")


    def set_db(self):

        """
         保存结果标准化数据到mongo数据库中
        """
        s3_path = '{}/{}'.format(self._sheet.output, 'busco.pdf')
        busco = self.api.api('tool_lab.busco')
        self.busco_text_dict = dict()
        for i in glob.glob(self.busco.output_dir + '/*short_summary_busco_result*'):
            sample_name = os.path.basename(i).split('_short_summary_busco_result')[0]
            with open(i, 'r') as s:
                for line in s.readlines():
                    if line.startswith('#') or line.startswith('\n'):
                        continue
                    else:
                        if line.strip().startswith('C:'):
                            text = line.strip()
                            self.busco_text_dict.update({sample_name: text})
        self.busco_summary_dict = dict()
        for j in glob.glob(self.busco.output_dir + '/*summary_result*'):
            sample_name = os.path.basename(j).split('_summary_result')[0]
            self.busco_summary_dict.update({sample_name: j})




        #
        #
        # self.sample_name = os.path.basename(self.option('fa').path).split('.fa')[0]
        # summary_result = os.path.join(self.busco.output_dir, 'summary_result')
        # with open(os.path.join(self.busco.output_dir, 'short_summary_busco_result.txt'), 'r') as s:
        #     for line in s.readlines():
        #         if line.startswith('#') or line.startswith('\n'):
        #             continue
        #         else:
        #             if line.strip().startswith('C:'):
        #                 self.text = line.strip()

        busco.add_busco(self.option('main_id'), s3_path)
        self.end()

    def end(self):
        self.chart()
        for file in os.listdir(self.busco.output_dir):
            os.link(os.path.join(self.busco.output_dir, file), os.path.join(self.output_dir, file))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Busco分析文件",0],
            ["./*.pdf", "", "busco评估分析结果图",0],
            ["./*summary_result*", "xls", "busco评估分析结果表",0],
            ["./*short_summary_busco_result*", "xls", "busco评估分析详情表",0],
        ])
        super(BuscoWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'Busco{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.busco',
            'instant': False,
            'options': {
                'meta_file': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nomogram/test2medical/lung_test.txt',
                'exp_file': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nomogram/test2medical/count_test.txt',
                'factor_list': 'sex;age',
                'gene_list': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nomogram/test2medical/gene'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

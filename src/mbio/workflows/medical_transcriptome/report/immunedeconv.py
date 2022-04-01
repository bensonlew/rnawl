# -*- coding: utf-8 -*-
# __author__ = "chenyanyan, 2016.10.12"
# last_modify by khl 20170504

from biocluster.workflow import Workflow
import os, re, glob
import pandas as pd
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import unittest


class ImmunedeconvWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ImmunedeconvWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_file', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'species', 'type': 'string', 'default': 'Homo_sapiens'},
            {'name': 'method', 'type': 'string', 'default': 'xcell'},
            {'name': 'update_info', 'type': 'string'},
            {'name': "main_id", 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.immu = self.add_tool("medical_transcriptome.tool.immunedeconv")
        # self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/01 Cluster_Analysis')
        # self.inter_dirs = []

    # def send_log(self, data):
    #     # 中间目录修改
    #     m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
    #     region = m.group(1)
    #     inter_dir = m.group(2)
    #     self.logger.info("更新结果目录")
    #
    #     if "dirs" in data["data"]["sync_task_log"]:
    #         for dir_path in self.inter_dirs:
    #             dir_dict = {
    #                 "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
    #                 "size": "",
    #                 "format": dir_path[1],
    #                 "description": dir_path[2],
    #                 "region": region,
    #             }
    #             if len(dir_path) >= 5:
    #                 dir_dict.update({"code": "D" + dir_path[5]})
    #
    #             data["data"]["sync_task_log"]["dirs"].append(dir_dict)
    #     with open(self.work_dir + "/post.changed.json", "w") as f:
    #         json.dump(data, f, indent=4, cls=CJsonEncoder)
    #     super(GenesetClusterWorkflow, self).send_log(data)

    def run(self):
        self.run_immu()
        super(ImmunedeconvWorkflow, self).run()

    def run_immu(self):
        opts = {
            'exp_file': self.option('exp_file'),
            'species': self.option('species'),
            'method': self.option('method')
        }
        self.immu.set_options(opts)
        self.immu.on('end', self.set_db)
        self.immu.run()


    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        tool_immu = self.api.api("medical_transcriptome.tool_immu")
        # add result info
        tool_immu.add_immu(self.option('main_id'), self.immu.option('immu').path, self.option('method'))
        self.end()

    def end(self):
        # result_dir = self.add_upload_dir(self.tool.output_dir)
        # self.inter_dirs = [
        #     ["04 GeneSet", "", "基因集分析结果目录",0],
        #     ["04 GeneSet/01 Cluster_Analysis", "", "聚类分析", 0]
        # ]
        # result_dir.add_relpath_rules([
        #     [".", "", "聚类分析文件",0,"211530"],
        #     ["./seq.cluster_tree.txt", "txt", "基因/转录本聚类树文件",0,"211531"],
        #     ["./seq.kmeans_cluster.txt", "txt", "基因/转录本聚类树文件", 0],
        #     ["./sample.cluster_tree.txt", "txt", "样本聚类树文件",0,"211532"],
        #     ["./expression_matrix.xls", "xls", "聚类热图分析表",0,"211533"],
        #     ["./heatmap.pdf", 'pdf', "聚类热图",0],
        # ])
        # result_dir.add_regexp_rules([
        #     [r'.*subcluster_.*\.xls', 'xls', '子聚类分析表',0,"211534"],
        # ])
        super(ImmunedeconvWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.medical_transcriptome.report.immunedeconv import ImmunedeconvWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "estimate" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "medical_transcriptome.estimate",
            "output":"s3://medical_transcriptome/files/test/medical_transcriptome/medical_transcriptome/interaction_results/ASprofile_20200813_092731",
            "options": {
                'exp_file': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nomogram/test2medical/count_test.txt',
            }
        }

        wsheet = Sheet(data=data)
        wf =ImmunedeconvWorkflow(wsheet)
        wf.sheet.id = 'medical_transcriptome'
        wf.sheet.project_sn = 'medical_transcriptome'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
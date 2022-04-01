# -*- coding: utf-8 -*-
# __author__ = 'shenghe'

"""otu表的样本距离计算"""

import os
from biocluster.workflow import Workflow
import datetime


class PcaWorkflow(Workflow):
    """
    报告中调用otu计算样本距离时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PcaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", "format": "meta.otu.otu_table"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.task = self.add_tool('meta.beta_diversity.pca')
        self.output_dir = self.task
        options = {
            'otutable': self.option('otu_file')
        }
        self.task.set_options(options)
        self.task.run()
        super(PcaWorkflow, self).run()

    def set_db(self):
        """
        保存结果距离矩阵表到mongo数据库中
        """
        api_distance = self.api.distance
        matrix_path = self.output_dir + '/' + os.listdir(self.output_dir)[0]
        if not os.path.isfile(matrix_path):
            raise Exception("找不到报告文件:{}".format(matrix_path))
        params_json = {
            'otu_id': self.option('otu_id'),
            'level_id': self.option('level'),
            'distance_algorithm': self.option('method'),
            'task_type': self.option('task_type'),
            'submit_location': 'beta_sample_distance'
            }
        matrix_id = api_distance.add_dist_table(matrix_path,
                                                major=True,
                                                name='Distance_{}_{}'.format(self.option('method'),
                                                                             datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
                                                level=self.option('level'),
                                                otu_id=self.option('otu_id'),
                                                params=params_json)
        self.add_return_mongo_id('sg_beta_specimen_distance', matrix_id)
        self.logger.info(str(matrix_id))
        self.logger.info('运行self.end')
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "距离矩阵计算结果输出目录"],
            ["./%s" % os.listdir(self.output_dir)[0], "xls", "样本距离矩阵文件"],
            ])
        print self.get_upload_files()
        super(PcaWorkflow, self).end()

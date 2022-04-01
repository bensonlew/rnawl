# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import pandas as pd
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types
import time
from biocluster.file import getsize, exists
from biocluster.file import download
import unittest


class CorrNetworkWorkflow(Workflow):
    def __init__(self, wsheet):
        self._sheet = wsheet
        super(CorrNetworkWorkflow, self).__init__(wsheet)
        options = [
            {'name': 'source', 'type': 'string', 'default': 'tool_lab'},
            {'name': 'relate_id', 'type': 'string'},
            {'name': 'task_id', 'type': 'string'},
            {"name": "table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "co_type", "type": "string", "default": "pearson"},
            {"name": "co_value", "type": "float", "default": 0.08},
            {"name": "p_value", "type": "float", "default": 0.05},
            {"name": "top_abundance", "type": "int", "default": 50},
            {"name": "strategy", "type": "string", "default": "row"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {'name': 'project_task_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(wsheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.corr_network = self.add_module('tool_lab.corr_network')

    def run(self):
        self.corr_network.on('end', self.set_db)
        if self.option('source') == 'project':
            self.file_path = self.check_file_path()
            print self.file_path
            table_ = self.download_s3_file(self.file_path, 'network.txt')
            print table_
        if self.option('source') == 'tool_lab':
            table_ = self.option('table').prop['path']
        self.t_table(self.option('strategy'), table_)
        self.run_corr_network()
        super(CorrNetworkWorkflow, self).run()

    def check_file_path(self):
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        conn_upset = db[collection_name]
        status = 'start'
        count_time = 0
        while status == 'start':
            if count_time > 600:
                self.set_error('超过十分钟还没有结果文件生成，请检查是否生成文件时报错')
                break
            time.sleep(10)
            print 'sleep 10s'
            try:
                upset = conn_upset.find_one(
                    {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
                status = upset['status']
            except:
                pass
            count_time += 10
        upset = conn_upset.find_one(
            {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
        file_path = upset['file_path']
        return file_path

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        elif os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path), code = '13700502')
        return to_path

    def t_table(self, top='row', table=None):
        tb = pd.read_csv(table, sep='\t', index_col=0)
        if top == 'col':
            tb = tb.transpose()
        sorts = tb.apply(sum, axis=1).sort_values(ascending=False)
        selects = sorts.head(self.option('top_abundance')).index
        tb = tb.loc[selects, ]
        tb.to_csv(self.work_dir + '/new_table.txt', sep='\t', index_label='#id')
        self.table = self.work_dir + '/new_table.txt'
        with open(self.work_dir + '/abundance_rank.txt', 'w') as w:
            for i in zip(selects, range(1, len(selects) + 1)):
                w.write('{}\t{}\n'.format(*i))

    def run_corr_network(self):
        opts = {
            'otutable': self.table,
            'method': self.option('co_type'),
            'coefficient': self.option('co_value'),
            'significance': self.option('p_value'),
        }
        self.corr_network.set_options(opts)
        self.corr_network.run()

    def set_db(self):
        api = self.api.api('tool_lab.corr_network')
        file_list = [
            self.corr_network.output_dir + '/corr_network_calc/corr_network_centrality.txt', self.work_dir + '/abundance_rank.txt',
            self.corr_network.output_dir + '/corr_network_calc/corr_network_node_degree.txt',
            self.corr_network.output_dir + '/corr_network_calc/corr_network_clustering.txt',
            self.corr_network.output_dir + '/corr_network_calc/corr_network_by_cut.txt',
            self.corr_network.output_dir + '/otu_association/shared.0.03.{}.corr'.format(self.option('co_type'))
                ]
        api.add(self.option('main_id'), 'corr_network', file_list)
        self.end()

    def end(self):
        repaths = [
            [".", "", "物种相关性网络结果输出目录"],
            ["shared.txt", "txt", "shared文件"],
            ["corr_network_attributes.txt", "txt", "网络的单值属性表"],
            ["corr_network_by_cut.txt", "txt", "相关系数筛选后网络边文件"],
            ["corr_network_centrality.txt", "txt", "网络节点的中心系数表"],
            ["corr_network_clustering.txt", "txt", "网络节点的聚类系数表"],
            ["corr_network_degree_distribution.txt", "txt", "网络节点的度分布表"],
            ["corr_network_node_degree.txt", "txt", "网络节点的度统计表"]
        ]
        regexps = [
             [r".*\.otu\.corr", "corr", "物种相似性网络边文件"]
        ]
        sdir = self.add_upload_dir(self.corr_network.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(CorrNetworkWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.tool_lab.corr_network import CorrNetworkWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'Corr_network{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.corr_network',
            'instant': False,
            'options': {
                'source': 'project',
                'project_task_id': 'gbfl_1pln6ksbdhrmsc3qrd1u4m',
                'relate_id': '60d1a0a417b2bf1bedddc3bf',
                "co_type": "spearman",
                "co_value": 0.5,
                "p_value": 0.05,
                "strategy": 'row',
                "top_abundance": 50
                # 'strategy': 'col'
                # 'source': 'tool_lab',
                # 'table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/corr_network/otu_table.xls'
            }
        }
        wsheet = Sheet(data=data)
        wf =CorrNetworkWorkflow(wsheet)
        wf.sheet.id = 'batch_effect'
        wf.sheet.project_sn = 'batch_effect'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

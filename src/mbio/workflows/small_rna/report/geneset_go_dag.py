# -*- coding: utf-8 -*-
# __author__ = 'liubinxu, qinjincheng'

from biocluster.workflow import Workflow
from bson.objectid import ObjectId
import glob
import os
import unittest
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class GenesetGoDagWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetGoDagWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'go_enrich_detail', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'significant_diff', 'type': 'string', 'default': 'pvalue'},
            {'name': 'significant_value', 'type': 'string', 'default': '0.05'},
            {'name': 'top_num', 'type': 'string', 'default': '10'},
            {'name': 'go_list', 'type': 'string', 'default': None},
            {'name': 'go_enrich_id', 'type': 'string', "default": None},

            {"name": "task_id", "type": "string", "default": None},
            {"name": "submit_location", "type": "string", "default": None},
            {"name": "task_type", "type": "string", "default": None},

            {'name': 'update_info', 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/07 GO_Dag')
        self.inter_dirs = []
        self.tool = self.add_tool('small_rna.geneset.go_dag')
        self.api_geneset_go_dag = self.api.api('small_rna.geneset_go_dag')

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(GenesetGoDagWorkflow, self).send_log(data)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for o in self.sheet.options():
            self.logger.debug('{} - {}'.format(o, self.option(o)))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def run(self):
        self.tool.on('end', self.set_output)
        self.get_run_log()
        self.run_tool()
        super(GenesetGoDagWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("small_rna", table="sg_geneset_go_dag", main_id=self.option('go_enrich_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_tool(self):
        options = {
            'go_enrich_detail': self.option('go_enrich_detail').prop['path'],
            'significant_diff': self.option('significant_diff'),
            'significant_value': self.option('significant_value'),
            'top_num': int(self.option('top_num')),
            'go_list': self.option('go_list'),
        }
        self.tool.set_options(options)
        self.tool.run()

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        for source in glob.glob(os.path.join(self.tool.output_dir, '*')):
            basename = os.path.basename(source)
            link_name = os.path.join(self.output_dir, basename)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.set_db()

    def set_db(self):
        self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        go_enrich_id = self.option('go_enrich_id')
        output_dir = self.get_workflow_output_dir()
        self.api_geneset_go_dag.update_geneset_go_dag(go_enrich_id=go_enrich_id,
                                                      output_dir=output_dir)
        self.logger.info('finish set_db at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录",0],
            ["04 GeneSet/07 GO_Dag", "", "靶基因GO富集有向无环图", 0]
        ]
        result_dir.add_relpath_rules([
            ['.', '', '靶基因GO富集有向无环图文件', 0],
            ["./go_dag.pdf", "", "GO富集有向无环图pdf", 0, "211543"],
            ["./go_dag.png", "", "GO富集有向无环图png", 0, "211544"],
            ["./go_dag.svg", "", "GO富集有向无环图svg", 0, "211545"],
            ["./go_dag.tsv", "", "GO层级关系", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        result_dir.add_regexp_rules([
            [r'go_dag\.*', '', 'GO富集有向无环图', 0],
        ])
        super(GenesetGoDagWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test_default(self):

        from mbio.workflows.small_rna.report.geneset_go_dag import GenesetGoDagWorkflow
        from biocluster.wsheet import Sheet
        import random
        import datetime

        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'output': 's3://smallrna/files/m_188/small_rna/small_rna/interaction_results/GenesetGoDag_{}'.format(datetime.datetime.now().strftime('%Y%m%d_%H%M%S')),
            'type': 'workflow',
            'name': 'small_rna.report.geneset_go_dag',
            'options': {
                'go_enrich_detail': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/geneset_go_dag/go_enrich_geneset.xls',
                'significant_diff': 'pvalue',
                'significant_value': '0.05',
                'top_num': '10',
                'go_enrich_id': '5bf390c0a4e1af62168452ff',
            }
        }

        wsheet = Sheet(data=data)
        wf = GenesetGoDagWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

    def test_go_list(self):

        from mbio.workflows.small_rna.report.geneset_go_dag import GenesetGoDagWorkflow
        from biocluster.wsheet import Sheet
        import random
        import datetime

        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'output': 's3://smallrna/files/m_188/small_rna/small_rna/interaction_results/GenesetGoDag_{}'.format(datetime.datetime.now().strftime('%Y%m%d_%H%M%S')),
            'type': 'workflow',
            'name': 'small_rna.report.geneset_go_dag',
            'options': {
                'go_enrich_detail': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/geneset_go_dag/go_enrich_geneset.xls',
                'go_list': 'GO:0044464;GO:0044424;GO:0005488;GO:0009987;GO:0065007;GO:0043226;GO:0044699;GO:0043227;GO:0043229;GO:0044763;GO:0043231;GO:0044444;GO:0044422;GO:0044446;GO:0005515;GO:0019222;GO:0043167;GO:0031323;GO:0060255;GO:0080090',
                'go_enrich_id': '5bf390c0a4e1af62168452ff',
            }
        }

        wsheet = Sheet(data=data)
        wf = GenesetGoDagWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

if __name__ == '__main__':
    unittest.main()
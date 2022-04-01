# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,qinjincheng'

from biocluster.workflow import Workflow
from mbio.packages.ref_rna_v3.functions import workfuncdeco
import os
import shutil
import unittest
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile,linkdir
import re
import json
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class RmatsModelWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RmatsModelWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 's3_file_list', 'type': 'infile', 'format': 'ref_rna_v3.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'sample.group_table'},
            {'name': 'control_table', 'type': 'infile', 'format': 'sample.control_table'},
            {'name': 'result_dir', 'type': 'infile', 'format': 'ref_rna_v3.common_dir'},
            {'name': 'gene_id', 'type': 'string', 'default': ''},
            {'name': 'event_type', 'type': 'string', 'default': ''},
            {'name': 'as_type', 'type': 'string', 'default': ''},
            {'name': 'main_id', 'type': 'string', 'default': ''},
            {'name': 'update_info', 'type': 'string', 'default': ''}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 AS/05 AS_diff_model')
        self.inter_dirs = []
        if self._sheet.rerun:
            self.logger.info("该交互为重运行项目,先删除mongo库内容")
            interactiondelete = InteractionDelete(bind_object=self, project_type="ref_rna_v2",
                                                  main_id=self.option('main_id'))
            interactiondelete.delete_interactions_records()

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
        super(RmatsModelWorkflow, self).send_log(data)

    @workfuncdeco
    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    @workfuncdeco
    def run(self):
        self.get_run_log()
        self.run_download()
        super(RmatsModelWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_splicing_rmats_model", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    @workfuncdeco
    def run_download(self):
        self.step.add_steps('download')
        self.download = self.add_module('ref_rna_v3.download')
        self.download.set_options({'s3_file_list': self.option('s3_file_list')})
        self.download.on('start', self.set_step, {'start': self.step.download})
        self.download.on('end', self.set_step, {'end': self.step.download})
        self.download.on('end', self.run_rmats_model)
        self.download.run()

    def run_rmats_model(self):
        name2path = {
            os.path.basename(bam.strip())[:-4]: bam.strip() for bam in
            open(self.download.option('my_file_list').path)
        }
        num, vs_list = self.option('control_table').get_control_info()
        group_spname = self.option('group_table').get_group_spname()
        rmats_output = self.option('result_dir').path
        event_file = os.path.join(
            rmats_output, '{}.MATS.{}.alter_id.txt'.format(self.option('event_type'), self.option('as_type'))
        )
        ctrl, test = vs_list[0]
        l1 = ','.join(sorted([sp for sp in group_spname[test]]))
        l2 = ','.join(sorted([sp for sp in group_spname[ctrl]]))
        b1 = ','.join(sorted([name2path[sp] for sp in group_spname[test]]))
        b2 = ','.join(sorted([name2path[sp] for sp in group_spname[ctrl]]))
        self.step.add_steps('rmats_model')
        self.rmats_model = self.add_tool('ref_rna_v3.structure.rmats_model')
        self.rmats_model.set_options({
            'gene_id': self.option('gene_id'),
            'event_file': event_file,
            'event_type': self.option('event_type'),
            'l1': l1,
            'l2': l2,
            'b1': b1,
            'b2': b2
        })
        self.rmats_model.on('start', self.set_step, {'start': self.step.rmats_model})
        self.rmats_model.on('end', self.set_step, {'end': self.step.rmats_model})
        self.rmats_model.on('end', self.set_output)
        self.rmats_model.run()

    def set_output(self, event):
        for basename in os.listdir(self.rmats_model.output_dir):
            source = os.path.join(self.rmats_model.output_dir, basename)
            link_name = os.path.join(self.output_dir, basename)
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.set_db()

    @workfuncdeco
    def set_db(self):
        self.database = self.api.api('ref_rna_v3.rmats_model')
        self.database.add_rmats_model(
            output_dir=self.output_dir,
            s3_output=self._sheet.output,
            main_id=self.option('main_id')
        )
        self.end()

    @workfuncdeco
    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["04 AS", "", "可变剪切分析结果目录",0],
            ["04 AS/04 AS_diff_model", "", "差异可变剪切事件模式图文件", 0]
        ]
        relpath = [
            ['.', '', '差异可变剪切事件模式图分析',0,"211642"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ]
        result_dir.add_relpath_rules(relpath)
        super(RmatsModelWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''
    def test(self):
        from mbio.workflows.ref_rna_v3.report.rmats_model import RmatsModelWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'ref_rna_v3.report.rmats_model',
            'options': {
                's3_file_list': '/mnt/ilustre/users/sanger-dev/workspace/20190618/RmatsModel_tsg_33555_2602_4955/bam.list',
                'group_table': '/mnt/ilustre/users/sanger-dev/workspace/20190618/RmatsModel_tsg_33555_2602_4955/group.txt',
                'control_table': '/mnt/ilustre/users/sanger-dev/workspace/20190618/RmatsModel_tsg_33555_2602_4955/control.txt',
                'result_dir': 's3://refrnav2/files/m_188/188_5c820f0e8b599/tsg_33555/interaction_results/Splicing_Con_12_vs_Vit_12_20190614_142049/Con_12_vs_Vit_12/',
                'gene_id': 'ENSG00000166579',
                'event_type': 'SE',
                'as_type': 'JC',
                'main_id': '5d08497a17b2bf643e1617b3'
            }
        }
        wsheet = Sheet(data=data)
        wf = RmatsModelWorkflow(wsheet)
        wf.sheet.id = 'tsg_33555_2602_4955'
        wf.sheet.project_sn = '188_5c820f0e8b599'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,qinjincheng'

import os
import shutil
import unittest

from biocluster.workflow import Workflow

from mbio.packages.ref_rna_v3.functions import workfuncdeco
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile,linkdir
from mbio.packages.ref_rna_v2.chart import Chart
import glob


class RmatsStatWorkflow(Workflow):
    '''
    last_modify: 2019.06.14
    '''

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RmatsStatWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'root', 'type': 'infile', 'format': 'ref_rna_v3.common_dir'},
            {'name': 'pvalue_fdr', 'type': 'string', 'default': None},
            {'name': 'fdr', 'type': 'float', 'default': None},
            {'name': 'psi', 'type': 'float', 'default': None},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self._sheet.rerun:
            self.logger.info("该交互为重运行项目,先删除mongo库内容")
            interactiondelete = InteractionDelete(bind_object=self, project_type="ref_rna_v2",
                                                  main_id=self.option('main_id'))
            interactiondelete.delete_interactions_records()


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
        self.run_rmats_stat()
        super(RmatsStatWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_splicing_rmats_stats", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    @workfuncdeco
    def run_rmats_stat(self):
        self.step.add_steps('rmats_stat')
        self.rmats_stat = self.add_tool('ref_rna_v3.structure.rmats_stat')
        self.rmats_stat.set_options({
            'root': self.option('root'),
            'pvalue_fdr': self.option('pvalue_fdr'),
            'fdr': self.option('fdr'),
            'psi': self.option('psi'),
            'main_id': self.option('main_id')
        })
        self.rmats_stat.on('start', self.set_step, {'start': self.step.rmats_stat})
        self.rmats_stat.on('end', self.set_step, {'end': self.step.rmats_stat})
        self.rmats_stat.on('end', self.set_output)
        self.rmats_stat.run()

    def set_output(self, event):
        for src in [self.rmats_stat.option('event_stats').path, self.rmats_stat.option('psi_stats').path]:
            dst = os.path.join(self.output_dir, os.path.basename(src))
            shutil.copy(src, dst)
            self.logger.info('succeed in copying {} to {}'.format(src, dst))
        else:
            self.set_db()

    @workfuncdeco
    def set_db(self):
        api = self.api.api('ref_rna_v3.rmats')
        api.add_sg_splicing_rmats_stats_for_controller(stat_id=self.option('main_id'), outpath=self.output_dir)
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        splice_diff = self.output_dir + "/event_stats.file.txt"
        splice_psi = self.output_dir + "/psi_stats.file.txt"
        chart.chart_splice_diff_stat(splice_diff, splice_psi, cmp_name=self.option("root").prop["path"].split("/")[-2])
        chart.to_pdf()

        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            linkfile(p, self.output_dir + "/" + os.path.basename(p))


    @workfuncdeco
    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        relpath = [
            ['.', '', '差异可变剪切事件统计结果目录', 0, "211643"],
            ['event_stats.file.txt', '', '差异可变剪切事件统计表', 0, "211644"],
            ['psi_stats.file.txt', '', '差异可变剪切模式变化统计表', 0, "211645"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ['*.pdf', 'txt', "差异可变剪切统计图", 0],
        ]
        result_dir.add_relpath_rules(relpath)
        super(RmatsStatWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.ref_rna_v3.report.rmats_stat import RmatsStatWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'ref_rna_v3.report.rmats_stat',
            'options': {
                'root': '/mnt/ilustre/users/sanger-dev/workspace/20190617/RmatsStat_tsg_33555_4445_6310/remote_input/root/Con_34_vs_Vit_34',
                'pvalue_fdr': 'fdr',
                'fdr': '0.05',
                'psi': '0.1',
                'main_id': '5d07476117b2bf4d142e33b9'
            }
        }
        wsheet = Sheet(data=data)
        wf = RmatsStatWorkflow(wsheet)
        wf.sheet.id = 'tsg_33555_4445_6310'
        wf.sheet.project_sn = '188_5c820f0e8b599'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

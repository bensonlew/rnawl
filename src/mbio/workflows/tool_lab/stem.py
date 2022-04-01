# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.workflow import Workflow
from mbio.packages.ref_rna_v2.functions import workfuncdeco
import shutil
import subprocess
import os
import unittest
from mbio.packages.ref_rna_v2.chart_advance import ChartAdvance
import pandas as pd
import glob



class StemWorkflow(Workflow):
    '''
    last_modify: 2019.06.18
    '''

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(StemWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'matrix', 'type': 'infile', 'format': 'ref_rna_v2.matrix'},
            {'name': 'geneset', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'pheno', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'method', 'type': 'string', 'default': 'SCM'},  # ['SCM', 'K']
            {'name': 'log', 'type': 'int', 'default': 10},
            #标准化方法["log_normal"，"normal","no_normal"]
            {"name" :"normalized_method",'type': 'string', 'default': 'log_normal'},
            {'name': 'number', 'type': 'int', 'default': 50},
            {'name': 'unit', 'type': 'int', 'default': 2},
            {'name': 'significance', 'type': 'float', 'default': 0.05},
            {'name': 'correction', 'type': 'string', 'default': 'Bonferroni'},  # ['Bonferroni', 'FDR', 'None']
            {'name': 'clusters', 'type': 'int', 'default': 10},
            {'name': 'starts', 'type': 'int', 'default': 20},
            {'name': 'detail', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'cluster', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

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
        self.run_prepare()
        super(StemWorkflow, self).run()

    @workfuncdeco
    def run_prepare(self):
        self.step.add_steps('prepare')
        self.prepare = self.add_tool('tool_lab.timeseries.prepare')
        options = {
            'matrix': self.option('matrix'),
            'geneset': self.option('geneset'),
            'pheno': self.option('pheno'),
            'method': self.option('method'),
            # 'log': self.option('log'),
            'normalized_method':self.option('normalized_method'),
            'number': self.option('number'),
            'unit': self.option('unit'),
            'significance': self.option('significance'),
            'correction': self.option('correction'),
            'clusters': self.option('clusters'),
            'starts': self.option('starts'),
        }
        self.prepare.set_options(options)
        self.prepare.on('start', self.set_step, {'start': self.step.prepare})
        self.prepare.on('end', self.set_step, {'end': self.step.prepare})
        self.prepare.on('end', self.run_stem)
        self.prepare.run()

    def run_stem(self):
        self.java = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/miniconda2/bin/java')
        self.script = {
            'stem': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v2/stem/stem.jar'),
            'bash': os.path.join(self.work_dir, 'stem.sh')
        }
        self.file = {
            'identity_file': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v2/stem/identity_file'),
            'setting': self.prepare.option('setting').path,
            'genetable': os.path.join(self.work_dir, 'setting_genetable.txt'),
            'profiletable': os.path.join(self.work_dir, 'setting_profiletable.txt'),
            'kmeansclustertable': os.path.join(self.work_dir, 'setting_kmeansclustertable.txt')
        }
        os.chdir(self.work_dir)
        shutil.copy(self.file['setting'], 'setting.txt')
        cmd = '{} -mx4096M -jar {} -b {} {}'.format(
            self.java, self.script['stem'], 'setting.txt', self.work_dir
        )
        self.logger.debug('command of stem:\n{}'.format(cmd))
        hostname = {'sanger-dev': '192.168.12.101', 'sanger': '10.2.0.110', 'isanger': '10.2.0.115'}[
            self.config.wpm_user]
        self.logger.info('#' * 64)
        self.logger.info('start iteration of running STEM before it is finished')
        flag = True
        n = {'sanger-dev': 70, 'sanger': 60, 'isanger': 50}[self.config.wpm_user]
        while flag and n:
            n -= 1
            s = '''ssh -i {} {}@{} "cd {}; export DISPLAY=localhost:{}.0; {}; exit"\n'''.format(
                self.file['identity_file'], self.config.wpm_user, hostname, self.work_dir, n, cmd
            )
            open(self.script['bash'], 'w').write(s)
            spc = subprocess.Popen(
                'sh {}'.format(self.script['bash']), shell=True,
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            ret = spc.wait()
            self.logger.debug('stdout of stem:\n{}'.format(spc.stdout.read()))
            self.logger.debug('stderr of stem:\n{}'.format(spc.stderr.read()))
            if ret:
                self.logger.debug('fail to run stem at localhost:{}'.format(n))
            else:
                self.logger.info('succeed in running at localhost:{}'.format(n))
                flag = False
                self.logger.info('stop iteration after succeed in running STEM')
                self.logger.info('#' * 64)
        else:
            self.run_statistics()

    @workfuncdeco
    def run_statistics(self):
        self.step.add_steps('statistics')
        self.statistics = self.add_tool('tool_lab.timeseries.statistics')
        options = {
            'method': self.option('method'),
            'genetable': self.file['genetable'],
            'matrix': self.option('matrix')
        }
        if self.option('method') == 'SCM':
            options.update({'profiletable': self.file['profiletable']})
        elif self.option('method') == 'K':
            options.update({'kmeansclustertable': self.file['kmeansclustertable']})
        self.statistics.set_options(options)
        self.statistics.on('start', self.set_step, {'start': self.step.statistics})
        self.statistics.on('end', self.set_step, {'end': self.step.statistics})
        self.statistics.on('end', self.set_output)
        self.statistics.run()

    def set_output(self):
        self.logger.info('begin of the function (set_output) at ({})'.format(self.__class__.__name__))
        shutil.copy(self.statistics.option('detailtable').path, os.path.join(self.output_dir, 'detail.xls'))
        self.option('detail').set_path(os.path.join(self.output_dir, 'detail.xls'))
        shutil.copy(self.statistics.option('clustertable').path, os.path.join(self.output_dir, 'cluster.xls'))
        self.option('cluster').set_path(os.path.join(self.output_dir, 'cluster.xls'))
        self.logger.info('final of the function (set_output) at ({})'.format(self.__class__.__name__))
        self.set_db()

    @workfuncdeco
    def set_db(self):
        self.database = self.api.api('tool_lab.stem')
        pdf_s3_path = os.path.join(self._sheet.output, "stem.trend.pdf")
        self.database.add_stem(
            detail_table=self.option('detail').path,
            cluster_table=self.option('cluster').path,
            main_id=self.option('main_id'),
            s3_out_put=pdf_s3_path
        )
        self.modify_result()

    def modify_result(self):
        if self.option('method') == 'SCM':
            stat_pd = pd.read_table(self.option('cluster').path)
            stat_pd = stat_pd[["profile", "assigned", "cluster"]]
            stat_pd.columns = ["Profile_id", "members", "Cluster_ID"]
            stat_pd["members"] = stat_pd["members"].astype("int")
            stat_pd.to_csv(os.path.join(self.output_dir,"trend_stat.xls"),index=False,sep="\t")
            detail_df = pd.read_table(self.option('detail').path)
            list_pro = sorted(set(list(detail_df["profile"].tolist())))
            for i in list_pro:
                profile = i
                select_df = detail_df[detail_df["profile"] == i]
                select_df.to_csv(os.path.join(self.output_dir, "profile_{}_detail.xls".format(str(profile))), index=False, sep="\t")
        else:
            stat_pd = pd.read_table(self.option('cluster').path)
            stat_pd = stat_pd[["cluster", "number"]]
            stat_pd.columns = ["Cluster_ID", "members"]
            stat_pd["members"] = stat_pd["members"].astype("int")
            stat_pd.to_csv(os.path.join(self.output_dir, "trend_stat.xls"), index=False, sep="\t")
            detail_df = pd.read_table(self.option('detail').path)
            list_pro = sorted(set(list(detail_df["cluster"].tolist())))
            for i in list_pro:
                cluster = i
                select_df = detail_df[detail_df["cluster"] == i]
                select_df.to_csv(os.path.join(self.output_dir, "cluster_{}_detail.xls".format(str(cluster))),
                                 index=False, sep="\t")
        self.end()

    @workfuncdeco
    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r'trend_stat.xls', 'xls', '趋势分析统计表', 0],
            [r'profile.*.xls', 'xls', 'profile详情表', 0],
            [r'cluster.*.xls', 'xls', 'cluster详情表', 0],
            [r'run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        result_dir.add_regexp_rules([
            [r'.', 'xls', '时序差异分析详情信息', 0],
            [r'.*/trend_stat\.xls', 'xls', '趋势分析统计表', 0],
            [r'.*/cluster.*\.xls', 'xls', 'cluster详情表', 0],
            [r'.*/detail\.xls', 'xls', '趋势分析统计表', 0],
            [r'.*/profile.*\.xls', 'xls', 'profile详情表', 0],
            [r'.*/.*\.pdf', 'pdf', '时序趋势分析图', 0],

        ])
        super(StemWorkflow, self).end()

    def chart(self):
        chart = ChartAdvance()
        chart.work_dir = self.work_dir + "/"
        if self.option('method') == 'SCM':
            profile_list = glob.glob(self.output_dir + '/profile*detail.xls')
        elif self.option('method') == 'K':
            profile_list = glob.glob(self.output_dir + '/cluster*detail.xls')
        cluster_table = self.option('cluster').path
        chart.chart_stem_cluster(profile_list, cluster_table, self.option('method'))
        chart.to_pdf()

        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, os.path.join(self.output_dir, os.path.basename(p)))



class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test_scm(self):
        from mbio.workflows.tool_lab.stem import StemWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'tool_lab.stem',
            'options': {
                'matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/matrix.tsv',
                'geneset': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/geneset.list',
                'pheno': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/pheno.txt',
                'method': 'SCM',
                'number': 50,
                'unit': 2,
                'significance': 0.05,
                'correction': 'Bonferroni'
            }
        }
        wsheet = Sheet(data=data)
        wf = StemWorkflow(wsheet)
        wf.sheet.id = 'ref_rna_v2_upgrade'
        wf.sheet.project_sn = 'ref_rna_v2_upgrade'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

    def test_k(self):
        from mbio.workflows.tool_lab.stem import StemWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'tool_lab.stem',
            'options': {
                'matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/matrix.tsv',
                'geneset': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/geneset.list',
                'pheno': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/pheno.txt',
                'method': 'K',
                'clusters': 10,
                # 'starts': 20
            }
        }
        wsheet = Sheet(data=data)
        wf = StemWorkflow(wsheet)
        wf.sheet.id = 'tool_lab_test'
        wf.sheet.project_sn = 'tool_lab_test'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

    # def test_sanger(self):
    #     from mbio.workflows.ref_rna_v2.report.stem import StemWorkflow
    #     from biocluster.wsheet import Sheet
    #     import random
    #     data = {
    #         'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
    #         'type': 'workflow',
    #         'name': 'ref_rna_v2.report.stem',
    #         'options': {
    #             'matrix': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/stem/matrix.tsv',
    #             'pheno': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/stem/s3.pheno.txt',
    #             'method': 'SCM',
    #             'number': 50,
    #             'unit': 2,
    #             'significance': 0.05,
    #             'correction': 'Bonferroni'
    #         }
    #     }
    #     wsheet = Sheet(data=data)
    #     wf = StemWorkflow(wsheet)
    #     wf.sheet.id = 'i-sanger_155039'
    #     wf.sheet.project_sn = 'i-sanger_155039'
    #     wf.IMPORT_REPORT_DATA = False
    #     wf.IMPORT_REPORT_AFTER_DATA = False
    #     wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_scm')])
    unittest.TextTestRunner(verbosity=2).run(suite)

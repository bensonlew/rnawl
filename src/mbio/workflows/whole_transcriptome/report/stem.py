# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.workflow import Workflow
from mbio.packages.ref_rna_v2.functions import workfuncdeco
import shutil
import subprocess
import os
import unittest
import pandas as pd
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json
from mbio.packages.ref_rna_v2.chart_advance import ChartAdvance
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
            {'name': 'number', 'type': 'int', 'default': 50},
            {'name': 'unit', 'type': 'int', 'default': 2},
            {'name': 'significance', 'type': 'float', 'default': 0.05},
            {'name': 'correction', 'type': 'string', 'default': 'Bonferroni'},  # ['Bonferroni', 'FDR', 'None']
            {'name': 'clusters', 'type': 'int', 'default': 10},
            {'name': 'starts', 'type': 'int', 'default': 20},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'detail', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'cluster', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'group_dict', 'type': 'string'},
            {'name': 'update_info', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/07 Advanced_Analysis/05 Time_series_Analysis')
        self.inter_dirs = []

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
        super(StemWorkflow, self).send_log(data)

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
        super(StemWorkflow, self).send_log(data)

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
        self.run_prepare()
        super(StemWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="stem", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    @workfuncdeco
    def run_prepare(self):
        self.step.add_steps('prepare')
        self.option("correction", 'None')
        self.prepare = self.add_tool('ref_rna_v2.timeseries.prepare')
        options = {
            'matrix': self.option('matrix'),
            'geneset': self.option('geneset'),
            'pheno': self.option('pheno'),
            'method': self.option('method'),
            'log': self.option('log'),
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
        self.statistics = self.add_tool('ref_rna_v2.timeseries.statistics')
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
        '''
        shutil.copy(self.statistics.option('detailtable').path, os.path.join(self.output_dir, 'detail.tsv'))
        self.option('detail').set_path(os.path.join(self.output_dir, 'detail.tsv'))
        shutil.copy(self.statistics.option('clustertable').path, os.path.join(self.output_dir, 'cluster.tsv'))
        self.option('cluster').set_path(os.path.join(self.output_dir, 'cluster.tsv'))
        '''

        self.option('detail').set_path(self.statistics.option('detailtable').path)
        self.option('cluster').set_path(self.statistics.option('clustertable').path)
        self.logger.info('final of the function (set_output) at ({})'.format(self.__class__.__name__))
        self.set_db()

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
    def set_db(self):
        self.database = self.api.api('whole_transcriptome.stem')
        self.database.add_stem(
            detail_table=self.option('detail').path,
            cluster_table=self.option('cluster').path,
            main_id=self.option('main_id')
        )
        self.modify_result()

    @workfuncdeco
    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["07 Advanced_Analysis", "", "高级分析结果目录", 0],
            ["07 Advanced_Analysis/05 Time_series_Analysis", "", "时序分析结果目录", 0]
        ]
        result_dir.add_relpath_rules([
             [r'.', '', '时序表达趋势分析文件', 0]
         ])
        result_dir.add_regexp_rules([
            [r'trend_stat.xls', 'xls', '趋势分析统计表', 0],
            [r'profile.*.xls', 'xls', 'profile详情表', 0],
            [r'cluster.*.xls', 'xls', 'cluster详情表', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
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
        from mbio.workflows.ref_rna_v2.report.stem import StemWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'ref_rna_v2.report.stem',
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
        from mbio.workflows.ref_rna_v2.report.stem import StemWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'ref_rna_v2.report.stem',
            'options': {
                'matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/matrix.tsv',
                'geneset': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/geneset.list',
                'pheno': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/pheno.txt',
                'method': 'K',
                'clusters': 10,
                'starts': 20
            }
        }
        wsheet = Sheet(data=data)
        wf = StemWorkflow(wsheet)
        wf.sheet.id = 'ref_rna_v2_upgrade'
        wf.sheet.project_sn = 'ref_rna_v2_upgrade'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

    def test_sanger(self):
        from mbio.workflows.ref_rna_v2.report.stem import StemWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'ref_rna_v2.report.stem',
            'options': {
                'matrix': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/stem/matrix.tsv',
                'pheno': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/stem/s3.pheno.txt',
                'method': 'SCM',
                'number': 50,
                'unit': 2,
                'significance': 0.05,
                'correction': 'Bonferroni'
            }
        }
        wsheet = Sheet(data=data)
        wf = StemWorkflow(wsheet)
        wf.sheet.id = 'i-sanger_155039'
        wf.sheet.project_sn = 'i-sanger_155039'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_sanger')])
    unittest.TextTestRunner(verbosity=2).run(suite)

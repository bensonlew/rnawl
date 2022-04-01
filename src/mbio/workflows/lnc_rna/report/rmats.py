# -*- coding: utf-8 -*-
# __author__ = 'shicaiping, qinjincheng'

from biocluster.workflow import Workflow
import os
import shutil
import unittest
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json
import re
from biocluster.config import Config


class RmatsWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RmatsWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 's3_file_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'sample.group_table'},
            {'name': 'control_table', 'type': 'infile', 'format': 'sample.control_table'},
            {'name': 'ref_gtf', 'type': 'string', 'default': ''},
            # ['paired', 'single']
            {'name': 'read_type', 'type': 'string', 'default': 'paired'},
            # ['fr-unstranded', 'fr-firststrand', 'fr-secondstrand']
            {'name': 'lib_type', 'type': 'string', 'default': 'fr-unstranded'},
            {'name': 'splicing_id', 'type': 'string', 'default': ''},
            {'name': 'type_file', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'update_info', 'type': 'string', 'default': ''}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/05 AS/AS_diff')
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
        super(RmatsWorkflow, self).send_log(data)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        if not os.path.exists(self.option('ref_gtf')):
            self.logger.info("before: {}".format(self.option('ref_gtf')))
            ref_gtf = os.path.join(Config().SOFTWARE_DIR, self.option('ref_gtf').split("/app/")[1])
            self.option('ref_gtf', ref_gtf)
            self.logger.info("after: {}".format(self.option('ref_gtf')))
        self.get_run_log()
        self.run_download()
        super(RmatsWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_splicing_rmats", main_id=self.option('splicing_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_download(self):
        self.step.add_steps('download')
        self.download = self.add_module('lnc_rna.download')
        samples = list()
        target_list = os.path.join(self.work_dir, "target_list")
        with open(self.option("group_table").prop["path"], "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                sample = line.strip().split("\t")[0]
                if sample not in samples:
                    samples.append(sample)
        with open(self.option("s3_file_list").prop["path"], "r") as f, open(target_list, "w") as w:
            for line in f:
                for sample in samples:
                    if sample in line:
                        w.write(line)
        options = {
            's3_file_list': target_list
        }
        self.download.set_options(options)
        self.download.on('start', self.set_step, {'start': self.step.download})
        self.download.on('end', self.set_step, {'end': self.step.download})
        self.download.on('end', self.run_rmats)
        self.download.run()

    def run_rmats(self):
        bam_list = self.download.option('local_file_list').path
        bam_loc = os.path.join(self.work_dir, 'bam_loc.txt')
        open(bam_loc, 'w').writelines(
            ['{}\t{}\n'.format(line.strip(), os.path.basename(line.strip())[:-4]) for line in open(bam_list)]
        )
        self.step.add_steps('rmats')
        self.rmats = self.add_module('lnc_rna.rmats')
        options = {
            'group_table': self.option('group_table'),
            'control_table': self.option('control_table'),
            'bam_loc': bam_loc,
            'ref_gtf': self.option('ref_gtf'),
            'read_type': self.option('read_type'),
            'lib_type': self.option('lib_type')
        }
        self.rmats.set_options(options)
        self.rmats.on('start', self.set_step, {'start': self.step.rmats})
        self.rmats.on('end', self.set_step, {'end': self.step.rmats})
        self.rmats.on('end', self.set_output)
        self.rmats.run()

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        for dirname in os.listdir(self.rmats.output_dir):
            src = os.path.join(self.rmats.output_dir, dirname)
            dst = os.path.join(self.output_dir, dirname)
            if os.path.isdir(dst):
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
            self.logger.info('succeed in copying {} to {}'.format(src, dst))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.set_db()

    def set_db(self):
        self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        self.database = self.api.api('lnc_rna.rmats')
        num, vs_list = self.option('control_table').get_control_info()
        group_spname = self.option('group_table').get_group_spname()
        for ctrl, test in vs_list:
            self.database.add_sg_splicing_rmats(
                splicing_id=self.option('splicing_id'),
                outpath=os.path.join(self.output_dir, '{}_vs_{}'.format(test, ctrl)),
                s3_output=self._sheet.output,
                gene_type_tsv=self.option('type_file').path
            )
        self.logger.info('finish set_db at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["05 AS", "", "可变剪切分析结果目录", 0],
            ["05 AS/AS_diff", "", "差异可变剪切分析文件", 0],
        ]
        result_dir.add_relpath_rules([
            [r'.', '', '差异可变剪切事件分析', 0, '211251'],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        result_dir.add_regexp_rules([
            [r'.*_vs_.*', '', '差异比较组别',0],
            [r'.*_vs_.*/all_events_detail_big_table.txt', '', '结果详情表', 0],
            [r'.*_vs_.*/event_stats.file.txt', '', '差异可变剪切事件统计表', 0],
            [r'.*_vs_.*/event_type.file.txt', '', '可变剪切事件类型统计表', 0],
            [r'.*_vs_.*/psi_stats.file.txt', '', '差异可变剪切模式变化统计表', 0],
            [r'.*_vs_.*/splicing_stats.xls', 'xls', '组内差异可变剪切事件统计表', 0],
            [r'.*_vs_.*/fromGTF.(RI|A3SS|A5SS|SE|MXE).alter_id.txt', '', '全部可变剪接事件表', 0],
            [r'.*_vs_.*/fromGTF.novelEvents.(RI|A3SS|A5SS|SE|MXE).alter_id.txt', '', '新发现可变剪接事件表', 0],
            [r'.*_vs_.*/(RI|A3SS|A5SS|SE|MXE).MATS.JCEC.alter_id.psi_info.txt', '', '差异事件详情表（JCEC）', 0],
            [r'.*_vs_.*/(RI|A3SS|A5SS|SE|MXE).MATS.JC.alter_id.psi_info.txt', '', '差异事件详情表（JC）', 0],
            [r'.*_vs_.*/(RI|A3SS|A5SS|SE|MXE).MATS.JCEC.alter_id.txt', '', '事件详情表（JCEC）', 0],
            [r'.*_vs_.*/(RI|A3SS|A5SS|SE|MXE).MATS.JC.alter_id.txt', '', '事件详情表（JC）', 0],
            [r'.*_vs_.*/.*_splicing_JC_stats\.xls', 'xls', '组内差异可变剪切模式变化统计表(JC)', 0],
            [r'.*_vs_.*/.*_splicing_JCEC_stats\.xls', 'xls', '组内差异可变剪切模式变化统计表(JCEC)', 0],
        ])
        super(RmatsWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''
    def test(self):
        from mbio.workflows.lnc_rna.report.rmats import RmatsWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'lnc_rna.report.rmats',
            'options': {
                'group_table': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/rmats/group.txt',
                'control_table': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/rmats/control.txt',
                's3_file_list': '/mnt/ilustre/users/sanger-dev/workspace/20190313/SnpIndel_lnc_rna_6850_8188/s3_bam.list',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'read_type': 'paired',
                'lib_type': 'fr-unstranded',
                'control_id': '5c6dfb5517b2bf6465297fca',
                'group_id': '5c6dff2d17b2bf72d19b6db4',
                'submit_location': 'splicingrmats',
                'task_id': 'lnc_rna',
                'task_type': 2,
            }
        }
        wsheet = Sheet(data=data)
        wf = RmatsWorkflow(wsheet)
        wf.sheet.id = 'lnc_rna'
        wf.sheet.project_sn = 'lnc_rna'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

if __name__ == '__main__':
    unittest.main()
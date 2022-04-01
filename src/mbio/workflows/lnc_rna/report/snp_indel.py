# -*- coding: utf-8 -*-
# __author__ = 'chenyanyan, qindanhua, liubinxu, qinjincheng'

from biocluster.workflow import Workflow
import json
import time
import glob
import pandas as pd
import os
import shutil
import re
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.config import Config
from biocluster.api.file.lib.transfer import MultiFileTransfer
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder

class SnpIndelWorkflow(Workflow):
    '''
    last_modify: 2019.04.23
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SnpIndelWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'method_type', 'type': 'string', 'default': ''},
            {'name': 's3_file_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'ref_fa', 'type': 'string', 'default': ''},
            {'name': 'ref_gtf', 'type': 'string', 'default': ''},
            {'name': 'ref_genome', 'type': 'string', 'default': ''},
            {'name': 'des', 'type': 'string', 'default': ''},
            {'name': 'des_type', 'type': 'string', 'default': ''},
            {'name': 'main_id', 'type': 'string', 'default': ''},
            {'name': 'type_file', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'update_info', 'type': 'string', 'default': ''}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/06 SNP_InDel_Analysis')
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
        super(SnpIndelWorkflow, self).send_log(data)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} = {}'.format('method_type', self.option('method_type')))
        if self.option('s3_file_list').is_set:
            self.logger.debug('{} = {}'.format('s3_file_list', self.option('s3_file_list').path))
        if self.option('ref_fa'):
            self.logger.debug('{} = {}'.format('ref_fa', self.option('ref_fa')))
        if self.option('ref_gtf'):
            self.logger.debug('{} = {}'.format('ref_gtf', self.option('ref_gtf')))
        self.logger.debug('{} = {}'.format('ref_genome', self.option('ref_genome')))
        if self.option('des'):
            self.logger.debug('{} = {}'.format('des', self.option('des')))
        self.logger.debug('{} = {}'.format('des_type', self.option('des_type')))
        self.logger.debug('{} = {}'.format('main_id', self.option('main_id')))
        self.logger.debug('{} = {}'.format('type_file', self.option('type_file').path))
        self.logger.debug('{} = {}'.format('update_info', self.option('update_info')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        if not os.path.exists(self.option('ref_fa')):
            self.logger.info("before: {}".format(self.option('ref_fa')))
            ref_fa = os.path.join(Config().SOFTWARE_DIR, self.option('ref_fa').split("/app/")[1])
            self.option('ref_fa', ref_fa)
            self.logger.info("after: {}".format(self.option('ref_fa')))
        if not os.path.exists(self.option('ref_gtf')):
            self.logger.info("before: {}".format(self.option('ref_gtf')))
            ref_gtf = os.path.join(Config().SOFTWARE_DIR, self.option('ref_gtf').split("/app/")[1])
            self.option('ref_gtf', ref_gtf)
            self.logger.info("after: {}".format(self.option('ref_gtf')))
        if not os.path.exists(self.option('des')):
            self.logger.info("before: {}".format(self.option('des')))
            des = os.path.join(Config().SOFTWARE_DIR, self.option('des').split("/app/")[1])
            self.option('des', des)
            self.logger.info("after: {}".format(self.option('des')))
        self.get_run_log()
        self.run_download()
        super(SnpIndelWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_snp", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_download(self):
        self.step.add_steps('download')
        self.download = self.add_module('lnc_rna.download')
        options = {
            's3_file_list': self.option('s3_file_list')
        }
        self.download.set_options(options)
        self.download.on('start', self.set_step, {'start': self.step.download})
        self.download.on('end', self.set_step, {'end': self.step.download})
        self.download.on('end', self.run_snp_indel)
        self.download.run()

    def run_snp_indel(self):
        self.step.add_steps('snp_indel')
        self.snp_indel = self.add_module('lnc_rna.snp_indel')
        options = {
            'method_type': self.option('method_type'),
            'ref_fa': self.option('ref_fa'),
            'bam_list': self.download.option('local_file_list'),
            'ref_gtf': self.option('ref_gtf'),
            'ref_genome': self.option('ref_genome'),
            'des': self.option('des'),
            'des_type': self.option('des_type'),
        }
        self.snp_indel.set_options(options)
        self.snp_indel.on('start', self.set_step, {'start': self.step.snp_indel})
        self.snp_indel.on('end', self.set_step, {'end': self.step.snp_indel})
        self.snp_indel.on('end', self.set_output)
        self.snp_indel.run()

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        for basename in os.listdir(self.snp_indel.output_dir):
            source = os.path.join(self.snp_indel.output_dir, basename)
            link_name = os.path.join(self.output_dir, basename)
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.set_db()

    def set_db(self):
        self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        self.database = self.api.api('lnc_rna.snp_indel')
        self.my_upload_dir = os.path.join(self.work_dir, 'upload')
        if os.path.isdir(self.my_upload_dir):
            shutil.rmtree(self.my_upload_dir)
        os.mkdir(self.my_upload_dir)
        self.database.add_snp(
            module_output=self.output_dir,
            method_type=self.option('method_type'),
            upload_dir=self.my_upload_dir,
            main_id=self.option('main_id'),
            gene_type_tsv=self.option('type_file').path
        )
        self.logger.info('finish set_db at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.my_upload_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.my_upload_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.my_upload_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.my_upload_dir)
        self.inter_dirs = [
            ["06 SNP_InDel_Analysis", "", "SNP/InDel分析结果目录", 0],
        ]
        result_dir.add_relpath_rules([
            ['.', '', 'SNP/InDel分析文件', 0],
            ['snp_annotation.xls', '', 'SNP分析注释详情结果表格', 0],
            ['snp_transition_tranversion_statistics.xls', '', 'SNP类型统计结果表格', 0],
            ['snp_freq_statistics.xls', '', 'SNP频率统计结果表格', 0],
            ['snp_depth_statistics.xls', '', 'SNP深度统计结果表格', 0],
            ['snp_position_distribution.xls', '', 'SNP不同区域分布结果表格', 0],
            ['indel_position_distribution.xls', '', 'InDel不同区域分布结果表格', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        super(SnpIndelWorkflow, self).end()

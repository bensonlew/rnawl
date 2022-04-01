# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.workflow import Workflow
import os
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import json


class LncrnaSeqConsWorkflow(Workflow):
    '''
    last_modify: 2019.04.24
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(LncrnaSeqConsWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'lncrna_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'keep_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'old_genome', 'type': 'string', 'default': None},
            {'name': 'new_genome', 'type': 'string', 'default': None},
            {'name': 'type_tsv', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'main_id', 'type': 'string', 'default': ''},
            {'name': 'update_info', 'type': 'string', 'default': ''}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/07 Advanced_Analysis/07 Lncrna_seq_cons')
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
        super(LncrnaSeqConsWorkflow, self).send_log(data)

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
        self.get_run_log()
        self.run_lncrna_seq_cons()
        super(LncrnaSeqConsWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_lncrna_seq_cons", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_lncrna_seq_cons(self):
        self.step.add_steps('lncrna_seq_cons')
        self.lncrna_seq_cons = self.add_tool('lnc_rna.structure.lncrna_seq_cons')
        options = {
            'lncrna_gtf': self.option('lncrna_gtf'),
            'keep_list': self.option('keep_list'),
            'old_genome': self.option('old_genome'),
            'new_genome': self.option('new_genome'),
            'type_tsv': self.option('type_tsv')
        }
        self.lncrna_seq_cons.set_options(options)
        self.lncrna_seq_cons.on('start', self.set_step, {'start': self.step.lncrna_seq_cons})
        self.lncrna_seq_cons.on('end', self.set_step, {'end': self.step.lncrna_seq_cons})
        self.lncrna_seq_cons.on('end', self.set_output)
        self.lncrna_seq_cons.run()

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        for basename in os.listdir(self.lncrna_seq_cons.output_dir):
            source = os.path.join(self.lncrna_seq_cons.output_dir, basename)
            link_name = os.path.join(self.output_dir, basename)
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.set_db()

    def set_db(self):
        self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        self.database = self.api.api('lnc_rna.lncrna_seq_cons')
        self.database.add_lncrna_seq_cons(
            tabular=os.path.join(self.output_dir, 'lncRNA_sequence_conservation.tabular'),
            main_id=self.option('main_id')
        )
        self.logger.info('finish set_db at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))

        os.rename(self.output_dir + '/lncRNA_sequence_conservation.tabular', self.output_dir + '/lncRNA_sequence_conservation.xls')
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["07 Advanced_Analysis", "", "高级分析结果目录", 0],
            ["07 Advanced_Analysis/07 Lncrna_seq_cons", "", "lncRNA序列保守性结果目录", 0],
        ]
        result_dir.add_relpath_rules([
            ['.', '', 'lncRNA序列保守性文件', 0],
            ['./lncRNA_sequence_conservation.xls', '', 'lncRNA序列保守性结果表', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        super(LncrnaSeqConsWorkflow, self).end()

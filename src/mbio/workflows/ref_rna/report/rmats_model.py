# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/5/17 09:37

import re, os, Bio, argparse, sys, fileinput, urllib2

'''
跑事件模式图用到的工作流
'''

import json
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from bson import ObjectId
import datetime
import pandas as pd
import shutil, subprocess
import re
from biocluster.workflow import Workflow



class RmatsModelWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        self.new_options = {}
        super(RmatsModelWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "event_id", "type": "string"},
            {"name": "rmats_out_root_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {"name": "label_a", "type": "string"},
            {"name": "label_b", "type": "string"},
            {"name": "event_type", 'type': "string"},
            {"name": "intron_s", 'type': "int", "default": 1},
            {"name": "exon_s", 'type': "int", "default": 1},
            {"name": "splicing_id", 'type': "string"},
            {"name": "rmats_model_id", 'type': "string"},
            {"name": "update_info", 'type': "string"}
        ]
        self.logger.info(options)
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run_rmats_model(self):
        a_bam_str = ','.join([line.strip() for line in subprocess.check_output(
            os.path.join('ls {}'.format(os.path.join(self.option('rmats_out_root_dir').prop['path'], 'SAMPLE_1/*/*.bam'))),
            shell=True).strip().split(
            '\n')])
        b_bam_str = ','.join([line.strip() for line in subprocess.check_output(
            os.path.join('ls {}'.format(os.path.join(self.option('rmats_out_root_dir').prop['path'], 'SAMPLE_2/*/*.bam'))),
            shell=True).strip().split(
            '\n')])
        rmats_event_file = os.path.join(self.option("rmats_out_root_dir").prop['path'],
                                        "MATS_output/" + self.option(
                                            'event_type') + '.MATS.ReadsOnTargetAndJunctionCounts.alter_id.txt')
        tmp_event_file = os.path.join(self.option("rmats_out_root_dir").prop['path'],
                                      "MATS_output/" + self.option('event_id') + "." + self.option(
                                          'event_type') + '.MATS.ReadsOnTargetAndJunctionCounts.alter_id.txt')
        self.logger.info(rmats_event_file)
        self.logger.info(tmp_event_file)

        tmp_str = subprocess.check_output(
            "grep '%s' %s  " % (self.option('event_id'), rmats_event_file), shell=True).strip()

        if tmp_str:
            open(tmp_event_file, 'w').write(tmp_str)
        else:
            raise Exception('这个事件：{}没有模式信息记录！'.format(self.option('event_id')))

        self.new_options = {
            "a_bam_str": a_bam_str,
            "b_bam_str": b_bam_str,
            "label_a": self.option("label_a"),
            "label_b": self.option("label_b"),
            "event_type": self.option("event_type"),
            "event_file": tmp_event_file,
            "intron_s": self.option("intron_s"),
            "exon_s": self.option("exon_s")
        }
        self.rmats_model = self.add_tool("gene_structure.rmats_model")
        self.rmats_model.set_options(self.new_options)
        self.rmats_model.on('end', self.set_db)
        self.rmats_model.run()
    
    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        relpath = [
            [".", "", "结果输出目录"],
        
        ]
        result_dir.add_regexp_rules([
            [r"\S+\.pdf$", "xls", "rmats事件模式图"]
        ])
        result_dir.add_relpath_rules(relpath)
        super(RmatsModelWorkflow, self).end()
    
    def set_db(self):
        """
        保存结果表保存到mongo数据库中
        """
        api_rmats_model = self.api.refrna_rmats_model
        result_files = [os.path.join(self.rmats_model.output_dir + '/Sashimi_plot', f.strip()) for f in
                        os.listdir(self.rmats_model.output_dir + '/Sashimi_plot') if re.match('^\S+\.(pdf|png)$', f.strip())]
        self.logger.info("要向mongo数据库中导入的文件是： %s" % result_files)
        splicing_id = ObjectId(self.option('splicing_id'))
        self.logger.info("准备开始向mongo数据库中导入rmats 事件模式图的rmats_model和graph信息！")
        # rmats_model_id = api_rmats_model.add_rmats_model(params=self.new_options, splicing_id=splicing_id)
        rmats_model_id = self.option('rmats_model_id')
        
        # for pdf in result_files:
        # rmats_model_id = ObjectId(self.option('rmats_model_id'))
        # self.logger.info('得到的rmats_model主表id是：%s' % rmats_model_id)
        for pic in result_files:
            file_type = re.match('^\S+\.(png|pdf)$', os.path.basename(pic)).group(1)
            api_rmats_model.add_sg_fs(file_path=pic, rmats_model_id=ObjectId(rmats_model_id), file_type=file_type)
        self.logger.info("向mongo数据库中导入rmats 事件模式图的rmats_model和graph信息成功！")
        self.end()
        
        # self.update_rmats_model(table_id=self.option('diff_express_id'), compare_column=compare_column,
        #                          group_detail=self.group_spname, samples=self.samples)
    
    def run(self):
        self.run_rmats_model()
        super(RmatsModelWorkflow, self).run()

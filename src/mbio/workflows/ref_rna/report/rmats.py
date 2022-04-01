# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/5/17 08:57

import re, os, Bio, argparse, sys, fileinput

import json
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from bson import ObjectId
import datetime
import pandas as pd
import shutil, subprocess
import re
from biocluster.workflow import Workflow
import importlib
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.gene_structure.rmats_process_func import *
from mbio.packages.gene_structure.rmats_process_func import process_single_rmats_output_dir
import re
from mbio.files.sequence.file_sample import FileSampleFile
from mainapp.models.mongo.ref_rna import RefRna
from biocluster.file import getsize, exists
from biocluster.file import download
from mbio.packages.ref_rna_v2.copy_file import CopyFile
from biocluster.api.file.lib.transfer import MultiFileTransfer


class RmatsWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        # self.logger.info(wsheet_object.options)
        self.rpc = False
        super(RmatsWorkflow, self).__init__(wsheet_object)
        self.rmats = self.add_tool("gene_structure.rmats_bam")
        self.step.add_steps("rmats_bam")
        self.rmats.on("end", self.set_output)
    
    def check_options(self):
        pass

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
            self.set_error('file can not find')
        return to_path

    def download_bam(self):
        case_bam_download = list()
        control_bam_download = list()
        transfer = MultiFileTransfer()
        if os.path.exists(self.work_dir + "/bam_path"):
            shutil.rmtree(self.work_dir + "/bam_path")
        os.mkdir(self.work_dir + "/bam_path")
        for bam in self.option("case_group_bam_str").split(","):
            bam_basename = os.path.basename(bam.strip())
            transfer.add_download(bam.strip(), self.work_dir + "/bam_path/")
            transfer.perform()
            bam_new = self.work_dir + "/bam_path/" + bam_basename
            #bam_new = self.download_s3_file(bam.strip(), os.path.basename(bam.strip()))
            case_bam_download.append(bam_new)
        self.option("case_group_bam_str", ",".join(case_bam_download))
        for bam in self.option("control_group_bam_str").split(","):
            bam_basename = os.path.basename(bam.strip())
            transfer.add_download(bam.strip(), self.work_dir + "/bam_path/")
            transfer.perform()
            bam_new = self.work_dir + "/bam_path/" + bam_basename
            #bam_new = self.download_s3_file(bam.strip(), os.path.basename(bam.strip()))
            control_bam_download.append(bam_new)
        self.option("control_group_bam_str", ",".join(control_bam_download))
    
    def run_rmats(self):
        opts = {
            "A_group_bam": self.option("case_group_bam_str"),
            "B_group_bam": self.option("control_group_bam_str"),
            "lib_type": self.option("lib_type"),
            "ref_gtf": self.option("ref_gtf"),
            "seq_type": self.option("seq_type"),
            "read_length": self.option("read_length"),
            "novel_as": self.option("novel_as"),
            "cut_off": self.option("cut_off"),
            "keep_temp": self.option("keep_temp"),
            "analysis_mode": self.option("analysis_mode")
        }
        self.rmats.set_options(opts)
        self.rmats.run()
    
    def run(self):
        options = [
            {"name": "seq_type", "type": "string", "default": "paired"},  # 两个选项：'paired'  or ’single‘
            {"name": "analysis_mode", "type": "string", "default": "U"},
            {"name": "read_length", "type": "int", "default": 160},
            {"name": "ref_gtf", "type": "string"},  # 一定要设置
            {"name": "novel_as", "type": "int", "default": 1},  # 是否发现新的AS事件，默认为是
            {"name": "lib_type", "type": "string", "default": "fr-unstranded"},  # 建库类型
            {"name": "cut_off", "type": "float", "default": 0.05},
            {"name": "keep_temp", "type": "int", "default": 1},
            {"name": "update_info", "type": "string"},
            {"name": "case_group_bam_str", "type": "string"},
            {"name": "control_group_bam_str", "type": "string"},
            {"name": "case_group_name", "type": "string"},
            {"name": "control_group_name", "type": "string"},
            {"name": "splicing_id", "type": "string"},
            {"name": "control_id", "type": "string"},
            {"name": "control_file", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "group_detail", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.logger.info(self._sheet.options())
        self.download_bam()
        self.run_rmats()
        super(RmatsWorkflow, self).run()
    
    def test_run(self):
        self._work_dir = '/mnt/ilustre/users/sanger-dev/workspace/20170522/Rmats_tsg_1000_6149_6691'
        self._output_path = os.path.join(self._work_dir, 'output')
        options = [
            # {"name": "sample_bam_dir", "type": "infile", "format": "align.bwa.bam_dir"},
            # {"name": "rmats_control", "type": "infile", "format": "sample.control_table"},
            # {"name": "group_table", "type": "string"},
            # {"name": "group_detail", "type": "string"},
            # {"name": "gname", "type": "string", "default": "group"},  # 分组方案名称
            {"name": "seq_type", "type": "string", "default": "paired"},  # 两个选项：'paired'  or ’single‘
            {"name": "analysis_mode", "type": "string", "default": "U"},
            {"name": "read_length", "type": "int", "default": 160},
            {"name": "ref_gtf", "type": "string"},  # 一定要设置
            {"name": "novel_as", "type": "int", "default": 1},  # 是否发现新的AS事件，默认为是
            {"name": "lib_type", "type": "string", "default": "fr-unstranded"},  # 建库类型
            {"name": "cut_off", "type": "float", "default": 0.05},
            {"name": "keep_temp", "type": "int", "default": 1},
            {"name": "update_info", "type": "string"},
            # {"name": "chr_set", "type": "string"},
            {"name": "case_group_bam_str", "type": "string"},
            {"name": "control_group_bam_str", "type": "string"},
            {"name": "case_group_name", "type": "string"},
            {"name": "control_group_name", "type": "string"},
            {"name": "splicing_id", "type": "string"},
            {"name": "control_id", "type": "string"},
            {"name": "control_file", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "group_detail", "type": "string"}
            # {"name":"ref_gtf",}
            # {"name": "case_name", "type": "string"},
            # {"name": "control_name", "type": "string"}
        ]
        # self.logger.info(self.init_options)
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.logger.info(self._sheet.options())
        self.rmats = self.add_tool("gene_structure.rmats_bam")
        self.step.add_steps("rmats_bam")
        self.logger.info('set db')
        self.logger.info('结果目录为：%s' % os.path.join(self.output_dir, self.rmats.name))
        # process_single_rmats_output_dir(root=os.path.join(self.output_dir, self.rmats.name))
        self.set_db()
    
    def set_output(self):
        self.logger.info('set output')
        output_dir = os.path.join(self.output_dir, self.rmats.name)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        outfiles = os.listdir(self.rmats.output_dir)
        CopyFile().linkdir(self.rmats.output_dir, output_dir)
        '''
        for f in outfiles:
            if f.startswith("log.RNASeq-MATS"):
                continue
            f_path = os.path.join(self.rmats.output_dir, f)
            target = os.path.join(output_dir, f)

            os.symlink(f_path, target)
        '''
        process_single_rmats_output_dir(root=output_dir)
        self.set_db()
        self.end()
    
    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r'.', '', '可变剪切分析结果文件']
        ])
        result_dir.add_regexp_rules([
            ["fromGTF\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt", 'txt', '可变剪接事件基本表'],
            ["fromGTF\.novelEvents\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt", 'txt', '新发现可变剪接事件基本表'],
            ["(RI|A3SS|A5SS|SE|MXE)\.MATS\.ReadsOnTargetAndJunctionCounts\.alter_id\.psi_info\.txt", 'txt',
             '差异事件详情表（ReadsOnTargetAndJunctionCounts证据）'],
            ["(RI|A3SS|A5SS|SE|MXE)\.MATS\.JunctionCountOnly\.alter_id\.psi_info\.txt", 'txt',
             '差异事件详情表（JunctionCountOnly证据）'],
            ['all_events_detail_big_table.txt', 'txt', '全部结果整合大表'],
            ['config.txt', 'txt', '运行配置详情文件'],
            ['all_events_detail_big_table.txt', 'txt', '结果综合详情表']
        ])
        super(RmatsWorkflow, self).end()
    
    def set_db(self):
        """
        保存结果表保存到mongo数据库中
        """
        
        api_rmats = self.api.refrna_splicing_rmats
        self.logger.info("准备开始向mongo数据库中导入rmats分析的信息！")
        rmats_out_root_dir = os.path.join(self.output_dir, self.rmats.name)
        self.logger.info('结果根目录为： %s' % rmats_out_root_dir)
        api_rmats.db['sg_splicing_rmats'].update({'_id': ObjectId(self.option('splicing_id'))},
                                                 {'$set': {'rmats_out_root_dir': rmats_out_root_dir}})
        self.logger.info('结果根目录为： %s' % rmats_out_root_dir)
        api_rmats.add_sg_splicing_rmats_for_controller(splicing_id=ObjectId(self.option('splicing_id')),
                                                       outpath=rmats_out_root_dir)
        self.logger.info("向mongo数据库中导入rmats的信息成功！")
        # self.end()

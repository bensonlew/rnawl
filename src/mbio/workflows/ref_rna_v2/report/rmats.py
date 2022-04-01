# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/5/17 08:57
# last modified by shicaiping at 20180511

from bson import ObjectId
from biocluster.workflow import Workflow
import os
from mbio.packages.ref_rna_v2.rmats_process_func import *
from mbio.packages.ref_rna_v2.rmats_process_func import process_single_rmats_output_dir
import re,shutil
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.ref_rna_v2.chart import Chart


class RmatsWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(RmatsWorkflow, self).__init__(wsheet_object)
        self.rmats = self.add_tool("ref_rna_v2.rmats_bam")
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
            if not os.path.exists(bam):
                bam_basename = os.path.basename(bam.strip())
                transfer.add_download(bam.strip(), self.work_dir + "/bam_path/")
                transfer.perform()
                bam_new = self.work_dir + "/bam_path/" + bam_basename
            else:
                bam_new = bam
            #bam_new = self.download_s3_file(bam.strip(), os.path.basename(bam.strip()))
            case_bam_download.append(bam_new)
        self.option("case_group_bam_str", ",".join(case_bam_download))
        for bam in self.option("control_group_bam_str").split(","):
            if not os.path.exists(bam):
                bam_basename = os.path.basename(bam.strip())
                transfer.add_download(bam.strip(), self.work_dir + "/bam_path/")
                transfer.perform()
                bam_new = self.work_dir + "/bam_path/" + bam_basename
            else:
                bam_new = bam
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
        }
        self.rmats.set_options(opts)
        self.rmats.run()
    
    def run(self):
        options = [
            {"name": "seq_type", "type": "string", "default": "paired"},  # 两个选项：'paired'  or ’single‘
            {"name": "read_length", "type": "int", "default": 150},
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 一定要设置
            {"name": "lib_type", "type": "string", "default": "fr-unstranded"},  # 建库类型
            {"name": "update_info", "type": "string"},
            {"name": "case_group_bam_str", "type": "string"},
            {"name": "control_group_bam_str", "type": "string"},
            {"name": "case_group_name", "type": "string"},
            {"name": "control_group_name", "type": "string"},
            {"name": "control_id", "type": "string"},
            {"name": "control_file", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "group_dict", "type": "string"},
            {"name": "splicing_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.logger.info(self._sheet.options())
        self.get_run_log()
        self.download_bam()
        self.run_rmats()
        super(RmatsWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_splicing_rmats", main_id=self.option('splicing_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()
    
    def set_output(self):
        self.logger.info('set output')
        outfiles = os.listdir(self.rmats.output_dir)
        for f in outfiles:
            f_path = os.path.join(self.rmats.output_dir, f)
            target = os.path.join(self.output_dir, f)
            if os.path.exists(target):
                os.remove(target)
            os.link(f_path, target)
        self.set_db()
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"

        chart.chart_geneset_enrich_circ(sys.argv[1], sys.argv[2])
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")[0]
        os.link(pdf_file, self.tool.output_dir + "/sample_correlation.pdf")
    
    def end(self):
        # self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r'.', '', '可变剪切分析结果文件', 0,  "211251"]
        ])
        result_dir.add_regexp_rules([
            ["fromGTF\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt", 'txt', '全部可变剪接事件表', 0, "211227"],
            ["fromGTF\.novelEvents\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt", 'txt', '新发现可变剪接事件表', 0, "211228"],
            ["(RI|A3SS|A5SS|SE|MXE)\.MATS\.JCEC\.alter_id\.psi_info\.txt", 'txt',
             '差异事件详情表（JCEC）', 0, "211229"],
            ["(RI|A3SS|A5SS|SE|MXE)\.MATS\.JC\.alter_id\.psi_info\.txt", 'txt',
             '差异事件详情表（JC）', 0, "211230"],
            ["(RI|A3SS|A5SS|SE|MXE)\.MATS\.JCEC\.alter_id\.txt", 'txt',
             '事件详情表（JCEC）', 0, "211231"],
            ["(RI|A3SS|A5SS|SE|MXE)\.MATS\.JC\.alter_id\.txt", 'txt',
             '事件详情表（JC）', 0, "211232"],
            ['all_events_detail_big_table.txt', 'txt', '结果详情表', 0, "211233"],
            ['psi_stats.file.txt', 'txt', '差异可变剪切模式变化统计表', 0, "211234"],
            # ['*.pdf', 'pdf', '差异可变剪切模式统计图', 0],
            ['event_type.file.txt', 'txt', '可变剪切事件类型统计表', 0, "211235"],
            ['event_stats.file.txt', 'txt', '差异可变剪切事件统计表', 0, "211236"],
            ['A_group_bam.txt', 'txt', '', 0, "211252"],
            ['B_group_bam.txt', 'txt', '', 0, "211237"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(RmatsWorkflow, self).end()
    
    def set_db(self):
        """
        保存结果表保存到mongo数据库中
        """
        api_rmats = self.api.api('ref_rna_v2.splicing_rmats')
        self.logger.info("准备开始向mongo数据库中导入rmats分析的信息！")
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        rmats_out_root_dir = os.path.join(self.workflow_output)
        self.logger.info('结果根目录为： %s' % rmats_out_root_dir)
        api_rmats.db['sg_splicing_rmats'].update({'_id': ObjectId(self.option('splicing_id'))},
                                                 {'$set': {'result_dir': rmats_out_root_dir, 'main_id': ObjectId(self.option('splicing_id'))}}, upsert=True)
        api_rmats.add_sg_splicing_rmats_for_controller(splicing_id=ObjectId(self.option('splicing_id')),
                                                       outpath=self.rmats.output_dir)
        self.logger.info("向mongo数据库中导入rmats的信息成功！")

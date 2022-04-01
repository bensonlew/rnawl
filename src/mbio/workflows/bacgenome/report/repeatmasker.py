# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'
# version 1.0
# last_modify: 2020.10.22

from biocluster.workflow import Workflow
import os
import gevent
import shutil
from biocluster.file import download,exists
from pymongo import MongoClient
from mbio.packages.bacgenome.common import link_dir

class RepeatmaskerWorkflow(Workflow):
    """
        细菌基因组散在重复序列预测工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RepeatmaskerWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "specimen_id", "type": "string"},  # 样品名称
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.common_api = self.api.api('bacgenome.common_api')
        self.faa_dir = os.path.join(self.work_dir, "faa_dir")
        self.modules = []

    def down_load_files(self):

        self.logger.info("开始下载组装序列文件！")
        (self.analysis_type,self.sample_id) = self.common_api.get_seq_files(self.faa_dir, self.option("task_id"))
        self.logger.info("下载完成！")
        if self.analysis_type in ['complete']:
            self.cat_reads()

    def run(self):
        self.down_load_files()
        self.run_repeatmasker()
        super(RepeatmaskerWorkflow, self).run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        self.logger.info("开始导表！")
        self.up_repeatmasker = self.api.api('bacgenome.repeatmasker')
        repeat_id = self.up_repeatmasker.get_id(task_id = self.option('task_id')) # test
        for sample in self.sample_id:
            sample_dir = os.path.join(self.output_dir, sample)
            sample_path = os.path.join(sample_dir, sample + ".gff")
            if os.path.exists(sample_path):
                if self.analysis_type in ['uncomplete']:
                    self.up_repeatmasker.add_repeat_predict_detail(repeat_id, sample, sample_path,
                                                                   self.faa_dir + "/" + sample +"_scaf.fna")
                else:
                    self.up_repeatmasker.add_repeat_predict_detail(repeat_id, sample, sample_path,
                                                                  self.work_dir + "/" + "faa_dir2" + '/' + sample  + "/" + sample + "_scaf.fna")
        self.end()

    def cat_reads(self):
        """
        合并reads 完成图合并fasta文件，形成单个样本文件一个
        :return:
        """
        for sample in self.sample_id:
            faa_dir2 = os.path.join(self.work_dir, "faa_dir2/"+sample)
            if os.path.exists(faa_dir2):
                shutil.rmtree(faa_dir2)
            os.makedirs(faa_dir2)
            sample_path = os.path.join(faa_dir2, sample + "_scaf.fna")
            assemble_dir = os.path.join(self.faa_dir, sample)
            dir_list = os.listdir(assemble_dir)
            for file2 in dir_list:
                os.system("cat {} >> {}".format(os.path.join(assemble_dir, file2), sample_path))

    def end(self):
        repaths = [
            [".", "", "基因组散在重复序列预测结果目录",0,"130205"],
        ]
        regexps = [
            [r".*\.*.gff$", 'gff', '散在重复序列预测gff文件',0,"130206"],
            [r'.*\.*.out$', 'out', '散在重复序列预测结果详情表',0,"130207"],
            [r'.*\.*.tbl$', 'tbl', '散在重复序列预测结果统计表',0,"130208"],
        ]
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules(repaths)
        result_dir.add_regexp_rules(regexps)
        super(RepeatmaskerWorkflow, self).end()

    def set_output(self):
        """
        设置结果目录文件
        :return:
        """
        self.logger.info("开始设置结果文件目录！")
        for module in self.modules:
            for sample_dir in os.listdir(module.output_dir):
                print sample_dir
                sample_dir_path = module.output_dir + "/" + sample_dir
                output_sample = self.output_dir + "/" + sample_dir
                if os.path.exists(output_sample):
                    shutil.rmtree(output_sample)
                if os.path.exists(sample_dir_path):
                    if len(os.listdir(sample_dir_path)) != 0:
                        link_dir(sample_dir_path, output_sample)
        self.logger.info("设置结果文件目录完成!")
        self.set_db()

    def run_repeatmasker(self):
        for sample in self.sample_id:
            self.repeatmasker = self.add_module('bacgenome.repeatmasker')
            if self.analysis_type in ['uncomplete']:
                file = self.faa_dir + '/' +sample + "_scaf.fna"
                opts = {
                    "genome_fa": file,
                    "sample": sample,
                    "analysis": self.analysis_type,
                }
            else:
                faa_dir2 = os.path.join(self.work_dir, "faa_dir2/"+sample)
                file = os.path.join(faa_dir2, sample + "_scaf.fna")
                opts = {
                    "genome_fa": file,
                    "sample": sample,
                    "analysis": self.analysis_type,
                }
            self.repeatmasker.set_options(opts)
            self.modules.append(self.repeatmasker)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.set_output)
        else:
            self.modules[0].on('end', self.set_output)
        for module in self.modules:
            module.run()
            gevent.sleep(0)
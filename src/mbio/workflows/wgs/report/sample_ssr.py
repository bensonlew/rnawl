# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 20180425

import os
import re
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from biocluster.api.file.lib.transfer import MultiFileTransfer


class SampleSsrWorkflow(Workflow):
    """
    交互分析：样本基因组SSR分析及引物设计
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SampleSsrWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "specimen_names", "type": "string"},
            {"name": "tm1", "type": "float", "default": 57.0},  # float, Tm1 (℃)
            {"name": "tm2", "type": "float", "default": 63.0},  # float, Tm2 (℃),要大于tm1
            {"name": "product_size", "type": "string", "default": "300-500"},  # Product Size(bp),数值要为int,范围间用-分隔，多个用,分隔,如：100-300
            {"name": "primer_num", "type": "int", "default": 3},  # Max Pairs Primer Number,范围:[1,5]
            {"name": "update_info", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "project_type", "type": "string"},
            {'name': "clean_fastq_path", "type": "string"},
            {'name': 'is_bucket', "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.target_dir = ''
        self.fastq_list = ''

    def check_options(self):
        if not self.option("specimen_names"):
            raise OptionError("请设置specimen_names", code="14500901")
        if not self.option("task_id"):
            raise OptionError("请设置task_id", code="14500902")

    def get_fastq_list(self):
        """
        得到样本对应的fastq_l和fastq_r
        """
        sample_ssr_api = self.api.api("wgs.sample_ssr_analysis")
        if self.option("project_type"):
            sample_ssr_api._project_type = self.option("project_type")
        self.fastq_list = os.path.join(self.work_dir, "clean_fastq.list")
        self.logger.info(self.option("task_id"))
        # if self._sheet.client == 'client01':
        #     self.target_dir = "/mnt/ilustre/data"
        # else:
        #     self.target_dir = "/mnt/ilustre/tsanger-data"
        self.target_dir = self.option('clean_fastq_path')
        sample_ssr_api.get_specimen_fastq_list(task_id=self.option("task_id"),
                                               specimen_names=self.option("specimen_names"),
                                               fastq_list=self.fastq_list, target_dir=self.target_dir)
        if self.option('is_bucket') == 'true':
            self.download_from_s3_()

    def download_from_s3_(self):
        """
        从对象存储中下载文件到指定路径
        :return:
        """
        if not os.path.exists(self.work_dir + "/temp"):
            os.mkdir(self.work_dir + "/temp")
        self.logger.info("开始下载对象存储中的文件！")
        transfer = MultiFileTransfer()
        with open(self.fastq_list, 'r') as r, open(self.work_dir + "/new_clean_fastq.list", 'w') as w:
            data = r.readlines()
            for line in data:
                tmp = line.strip().split('\t')
                transfer.add_download(tmp[1], "{}/temp/".format(self.work_dir))
                transfer.add_download(tmp[2], "{}/temp/".format(self.work_dir))
                w.write("{}\t{}\t{}\n".format(tmp[0], self.work_dir + "/temp/{}".format(os.path.basename(tmp[1])),
                                              self.work_dir + "/temp/{}".format(os.path.basename(tmp[2]))))
        transfer.perform()
        self.logger.info("下载对象存储中的文件成功！")

    def run_ssr_primer(self):
        options = {
            "fastq_list": self.work_dir + "/new_clean_fastq.list" if self.option('is_bucket') == 'true' else
            self.fastq_list,
            "tm1": self.option("tm1"),
            "tm2": self.option("tm2"),
            "product_size": self.option("product_size"),
            "primer_num": self.option("primer_num")
        }
        self.ssr_primer = self.add_module("wgs.ssr_primer")
        self.ssr_primer.set_options(options)
        self.ssr_primer.on("end", self.set_output)
        self.ssr_primer.run()

    def set_output(self):
        for f in os.listdir(self.ssr_primer.output_dir):
            f1 = os.path.join(self.ssr_primer.output_dir, f)
            f2 = os.path.join(self.output_dir, f)
            if os.path.isfile(f1):
                if os.path.exists(f2):
                    os.remove(f2)
                os.link(f1, f2)
            else:
                if not os.path.exists(f2):
                    os.mkdir(f2)
                for f_ in os.listdir(f1):
                    f3 = os.path.join(f1, f_)
                    f4 = os.path.join(f2, f_)
                    if os.path.exists(f4):
                        os.remove(f4)
                    os.link(f3, f4)
        self.set_db()
        self.end()

    def set_db(self):
        """
        将结果保存到mongo
        """
        self.logger.info("将结果保存到mongo")
        sample_ssr_api = self.api.api("wgs.sample_ssr_analysis")
        if self.option("project_type"):
            sample_ssr_api._project_type = self.option("project_type")
        ssr_specimen_id = self.option("main_id")
        task_id = self.option("task_id")
        for f in os.listdir(self.output_dir):
            m = re.match(r"(.+).ssr.stat", f)
            # n = re.match(r"(.+).result", f)
            if m:
                specimen_id = m.group(1)
                ssr_stat = os.path.join(self.output_dir, f)
                sample_ssr_api.add_sg_ssr_specimen_stat(ssr_specimen_id, specimen_id, ssr_stat)
            # if n:
            #     specimen_id = n.group(1)
            #     ssr_detail = os.path.join(self.output_dir, f)
            #     sample_ssr_api.add_sg_ssr_specimen_detail(ssr_specimen_id, specimen_id, ssr_detail)
        ssr_compare = os.path.join(self.output_dir, "ssr_compare")
        for f in os.listdir(ssr_compare):
            n = re.match(r"(.+).result", f)
            if n:
                specimen_id = n.group(1)
                ssr_detail = os.path.join(ssr_compare, f)
                sample_ssr_api.add_sg_ssr_specimen_detail(ssr_specimen_id, specimen_id, ssr_detail)
        download_file = self._sheet.output.strip().split(':')[1] + "/ssr_compare"
        sample_ssr_api.update_ssr_output_dir(ssr_specimen_id, ssr_compare, download_file)

    def run(self):
        self.get_fastq_list()
        self.run_ssr_primer()
        super(SampleSsrWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SampleSsrWorkflow, self).end()

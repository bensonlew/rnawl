# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""从序列文件或者序列文件夹当中获取样本信息"""
import os
import gevent
import json
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.meta.common_function import copy_task


class SampleExtractWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SampleExtractWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_fasta", "type": "infile", "format": "sequence.fasta, sequence.fasta_dir"},
            {"name": "in_fastq", "type": "infile", "format": "sequence.fastq, sequence.fastq_dir"},
            {"name": "table_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "from_task", "type": "string", "default": ""},##来源于哪个task_id
            {"name": "query_id", "type": "string"}, #根据任务id和query_id确定样本检测的唯一性
            {"name": "task_id", "type": "string", 'default': ""}, ## 任务的id
            {"name": "file_name", "type": "string", 'default': ""}, ## 文件名称
            {"name": "main_id", "type": "string", 'default': ""}, ## 兼容之前版本的样本检测
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.seq_extract = None
        self.setted_option = ""


    def check_options(self):
        # if self.option("in_fasta").is_set and self.option("in_fastq").is_set
        # 由于多文件format的问题， 当两个文件format检测都没有通过的时候，
        # 这个option的值是bool， is_set属性不能使用
        # 所以在这里，当option的值是bool的时候， 就说明这个变量没有设置
        if not isinstance(self.option("in_fastq"), bool) and not isinstance(self.option("in_fasta"), bool):
            raise OptionError("不能同时设置参数in_fasta和in_fastq")
        if self.option("from_task") == "":
            if isinstance(self.option("in_fastq"), bool):
                try:
                    self.option("in_fasta").is_set
                    self.setted_option = "fasta"
                except Exception as e:
                    raise OptionError("参数in_fasta错误:{}".format(e))
            if isinstance(self.option("in_fasta"), bool):
                try:
                    self.option("in_fastq").is_set
                    self.setted_option = "fastq"
                except Exception as e:
                    raise OptionError("参数in_fastq错误:{}".format(e))

    def run_seq_extract(self):
        if self.setted_option == "fasta":
            opts = {
                "in_fasta": self.option("in_fasta")
            }
        elif self.setted_option == "fastq":
            opts = {
                "in_fastq": self.option("in_fastq"),
                "table_id": self.option("table_id")
            }
        self.seq_extract.set_options(opts)
        self.seq_extract.on("end", self.set_db)
        self.seq_extract.run()

    def copy_task(self):
        """
        根据from_task直接从MongoDB中copy结果，不需要样本检测
        不需要导表，直接copy完事就行了
        :return:
        """
        self.logger.info("开始根据from_task生成结果表：{}".format(self.option("from_task")))
        if self.option("task_id") != "":
            task_id = self.option("task_id")
        else:
            task_id = "_".join(self._sheet.id.split('_')[0:2])
        from_task = self.option("from_task")
        insert_id = self.option("table_id")
        copy_task(task_id, from_task, insert_id, "sg_sample_check_detail", "sg_sample_check", "check_id", query=self.option("query_id"))
        gevent.spawn_later(2, self.end)

    def set_db(self):
        """
        更改MongoDB的导入形式，增加detail表用于前端记录
        :return:
        """
        api_sample = self.api.sample_extract
        api_sample.update_sg_seq_sample(self.seq_extract.option("file_sample_list").prop["path"], self.option("main_id"))
        if self.option("file_name") != "":
            s3_origin = self.get_s3_name("fastq")
            self.logger.info("s3_origin:{}".format(s3_origin))
            api_sample.add_seq_sample_detail(self.seq_extract.option("file_sample_list").prop["path"], self.option("table_id"), real_name=self.option("file_name"), s3_name=s3_origin)
        else:
            s3_origin = self.get_s3_name("fastq_dir")
            api_sample.add_seq_sample_detail(self.seq_extract.option("file_sample_list").prop["path"], self.option("table_id"), s3_name=s3_origin)
        self.end()

    def run(self):
        if self.option("from_task") == "":
            self.check_options()  # 通常不应该直接在调用， 但这里有些特殊，需要确定in_fasta和in_fastq哪个进行了设定
        if self.setted_option == "fasta":
            self.seq_extract = self.add_tool("sequence.fasta_sample_extract")
        elif self.setted_option == "fastq":
            self.seq_extract = self.add_module("meta.sample_extract.sample_extract")
        if self.option("from_task") != "":
            self.copy_task()
        else:
            self.run_seq_extract()
        super(SampleExtractWorkflow, self).run()

    def get_s3_name(self,type):
        """
        根据输入的s3路径得到s3的文件名称
        :return:
        """
        s3_origin = {}
        if type in ['fastq']:
            data_json = os.path.join(self.work_dir, "data.json")
            with open(data_json, 'r') as f:
                js = f.read()
                json_dict = json.loads(js)
            in_fastq = str(json_dict['options']['in_fastq']).split("||")[-1]
            true_file_name = in_fastq
            if self.option("file_name") not in s3_origin:
                s3_origin[self.option("file_name")] = true_file_name

        elif type in ['fastq_dir']:
            api_sample = self.api.sample_extract
            task_id = "_".join(self._sheet.id.split('_')[0:2])
            fastq_list = api_sample.get_file_info(task_id, self.option("query_id"))
            for fastq_dict in fastq_list:
                new_file_name = fastq_dict['alias']
                true_file_name = fastq_dict['file_path']
                if new_file_name not in s3_origin:
                    s3_origin[new_file_name] = true_file_name
        return s3_origin

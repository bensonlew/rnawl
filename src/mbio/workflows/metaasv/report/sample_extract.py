# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import os,re
import shutil
import json
import gevent
from biocluster.file import download
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir,copy_task,check_file


class SampleExtractWorkflow(Workflow):
    """
    metaasv 样本检测  交互分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SampleExtractWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fastq_file", "type": "infile", "format": "sequence.fastq, sequence.fastq_dir"},
            {"name": "fastq_dir", "type": "infile", "format": "metaasv.json_table"},
            {"name": "fastq_type", "type": "string"},##类型fastq和fastq_dir
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "from_task", "type": "string", "default": ""},##来源于哪个task_id
            {"name": "query_id", "type": "string"}, #根据任务id和query_id确定样本检测的唯一性
            {"name": "file_name", "type": "string", "default": ""},## 选择fastq时的真实的文件名称
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option("main_id"):
            raise OptionError("没有设置参数main_id")
        # if self.option("fastq_type") not in ["fastq", "fastq_dir"]:
        #     raise OptionError("fastq_type参数不符合参数类型fastq和fastq_dir")

    def download_file(self):
        """
        根据前端传过来的fastq_dir参数下载文件夹
        :return:
        'fastq_dir': u'{"file_list":[{"alias":"NASA5_M5_Fec1.fq","file_path":"s3:\\/\\/commonbucket\\/files\\/m_188\\/188_5edf26082cf44\\/Gut_data\\/d6c0a7e50699ed3874a0c2891d444959.fq"}]]}
        """
        self.fastq_dir = os.path.join(self.work_dir, "input_dir")
        if os.path.exists(self.fastq_dir):
            shutil.rmtree(self.fastq_dir)
        os.mkdir(self.fastq_dir)
        with open(self.option("fastq_dir").prop["path"], "r") as f:
            params_fastq_dir = json.load(f)
        if "file_list" in params_fastq_dir:
            file_list = params_fastq_dir["file_list"]
            for file_dict in file_list:
                file_name = file_dict["alias"]
                file_path = file_dict["file_path"]
                old_file_path = "/".join(file_path.split("/")[0:-1] + [file_name])
                new_file_path = os.path.join(self.fastq_dir, file_name)
                if os.path.exists(new_file_path):
                    os.remove(new_file_path)
                if old_file_path.startswith("/mnt/"):
                    os.system("cp {} {}".format(old_file_path, new_file_path))
                else:
                    old_file_path = file_path
                    download(old_file_path, new_file_path)

    def run_sequence(self,dir):
        """
        对上传的文件夹进行解压，上传序列为压缩格式
        :return:
        """
        opts = {
            'fastq_dir': dir,
        }
        self.sequence.set_options(opts)
        self.sequence.on("end", self.run_seq_extract)
        self.sequence.run()

    def run_file_sequence(self,file):
        """
        对上传的文件夹进行解压，上传序列为压缩格式
        :return:
        """
        file_name = os.path.basename(file)
        file_name = file_name.strip(".gz") if re.search(r'\.gz', file_name) else  file_name.strip(".tar.gz")
        result_path = os.path.join(self.work_dir, "results")
        if os.path.exists(result_path):
            shutil.rmtree(result_path)
        os.mkdir(result_path)
        opts = {
            'fastq': file,
            'sample_name': file_name,
            'direction': 's',
            'result_path': result_path,
            'pipeline': 'metaasv',
        }
        self.sequence.set_options(opts)
        self.sequence.on("end", self.run_seq_extract)
        self.sequence.run()

    def run_seq_extract(self):
        """
        进行样本检测fastq和fastq_dir
        :return:
        """
        if self.option("fastq_type") in ["fastq_dir"]:
            self.seq_extract = self.add_module("meta.sample_extract.sample_extract")
            if self.is_zip:
                in_fastq = os.path.join(self.sequence.output_dir, 'data')
            else:
                in_fastq = self.fastq_dir
            opts = {
                "in_fastq": in_fastq
                # "file_list": self.option("fastq_dir")
            }
        elif self.option("fastq_type") in ["fastq"]:
            self.seq_extract = self.add_module("meta.sample_extract.sample_extract")
            if self.is_zip:
                dir = os.path.join(self.work_dir, 'results')
                for file in os.listdir(dir):
                    file_path = os.path.join(dir, file)
                    break
                in_fastq = file_path
            else:
                in_fastq = self.option("fastq_file").prop['path']
            opts = {
                "in_fastq": in_fastq
            }
        self.seq_extract.set_options(opts)
        self.seq_extract.on("end", self.set_db)
        self.seq_extract.run()

    def copy_task(self):
        """
        根据from_task直接从MongoDB中copy结果，不需要样本检测
        :return:
        """
        self.logger.info("开始根据from_task生成结果表：{}".format(self.option("from_task")))
        # try:
        task_id = "_".join(self._sheet.id.split('_')[0:2])
        from_task = self.option("from_task")
        insert_id = self.option("main_id")
        copy_task(task_id, from_task, insert_id, "sample_check_detail", "sample_check", "check_id", query=self.option("query_id"))
        gevent.spawn_later(1, self.end)
        # except:
        #     self.set_error("copy数据失败！")

    def set_db(self):
        """
        导入MongoDB, 不上传结果文件，只是导入数据
        :return:
        """
        api_sample = self.api.api("metaasv.sample_check")
        task_id = "_".join(self._sheet.id.split('_')[0:2])
        s3_origin = self.get_s3_name(self.option("fastq_type"))
        if self.option("fastq_type") in ['fastq']:
            api_sample.add_seq_sample_detail(self.seq_extract.option("file_sample_list").prop["path"], self.option("main_id"), real_name=self.option("file_name"), s3_name=s3_origin)
        else:
            api_sample.add_seq_sample_detail(self.seq_extract.option("file_sample_list").prop["path"], self.option("main_id"), s3_name=s3_origin)
        api_sample.delete_sg_status(task_id, self.option("query_id"))
        self.end()

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
            in_fastq = str(json_dict['options']['fastq_file']).split("||")[-1]
            true_file_name = in_fastq
            if self.option("file_name") not in s3_origin:
                s3_origin[self.option("file_name")] = true_file_name

        elif type in ['fastq_dir']:
            data_json = os.path.join(self.work_dir, "fastq.json")
            with open(data_json, 'r') as f:
                js = f.read()
                json_dict = json.loads(js)
            fastq_list = json_dict['file_list']
            for fastq_dict in fastq_list:
                new_file_name = fastq_dict['alias']
                true_file_name = fastq_dict['file_path']
                if new_file_name not in s3_origin:
                    s3_origin[new_file_name] = true_file_name
        return s3_origin

    def run(self):
        """
        运行
        :return:
        """
        if self.option("from_task") != "":
            self.copy_task()
        else:
            if self.option("fastq_dir").is_set:
                self.download_file()
            if self.option("fastq_type") in ["fastq_dir"]:
                self.is_zip = check_file('dir', self.fastq_dir)
                self.sequence = self.add_module('metaasv.reads_unzip')
                if self.is_zip:
                    self.run_sequence(self.fastq_dir)
                else:
                    self.run_seq_extract()
            elif self.option("fastq_type") in ["fastq"]:
                self.is_zip = check_file('file', self.option("fastq_file").prop["path"])
                self.sequence = self.add_tool('sequence.fastq_ungz')
                if self.is_zip:
                    self.run_file_sequence(self.option("fastq_file").prop['path'])
                else:
                    self.run_seq_extract()
        super(SampleExtractWorkflow, self).run()

    def end(self):
        """
        结束
        :return:
        """
        super(SampleExtractWorkflow, self).end()
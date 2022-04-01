# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
## @20200303
from biocluster.core.exceptions import OptionError
from mbio.packages.toolapps.common import link_dir
from biocluster.module import Module
import os,re
import shutil
import gevent


class TransgenescanModule(Module):
    """
    小工具：宏转录组TransGeneScan软件进行基因预测
    """
    def __init__(self, work_id):
        super(TransgenescanModule, self).__init__(work_id)
        options = [
            {'name': 'in_fasta', 'type': 'infile', 'format': 'toolapps.fasta_dir'},  # 输入的fasta序列文件夹
            {'name': 'fasta', 'type': 'infile', 'format': 'toolapps.fasta'},  # 输入的fasta序列
        ]
        self.add_option(options)
        self.step.add_steps('transgenescan')
        self.tool_list = []

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if (not self.option('in_fasta').is_set) and (not self.option("fasta").is_set):
            raise OptionError('请输入fasta序列文件或者fasta_dir序列文件夹', code="24400301")

    def run_transgenescan(self):
        """
        并行运行TransGeneScan基因预测
        :return:
        """
        self.logger.info("start transgenescan tool")
        if self.option("fasta").is_set:
            transgenescan = self.add_tool("predict.transgenescan")
            fasta_path = self.option("fasta").prop['fasta']
            basename = os.path.basename(fasta_path)
            sample_name = basename.strip(".fa")
            opts = {
                'input_genome': fasta_path,
                'sample': sample_name
                }
            transgenescan.set_options(opts)
            self.tool_list.append(transgenescan)
        else:
            fasta_path = self.option("in_fasta").prop['fasta_dir']
            sample_dict = self.get_list()
            self.logger.info("sample_dict: {}".format(sample_dict))
            for sample in sample_dict.keys():
                transgenescan = self.add_tool("predict.transgenescan")
                file_name = sample_dict[sample]
                file_path = os.path.join(fasta_path, file_name)
                opts = {
                    'input_genome': file_path,
                    'sample': sample
                    }
                transgenescan.set_options(opts)
                self.tool_list.append(transgenescan)
        if len(self.tool_list) > 1:
            self.on_rely(self.tool_list, self.set_output)
        else:
            self.tool_list[0].on("end", self.set_output)
        for tool in self.tool_list:
            tool.run()
            gevent.sleep(0)

    def run(self):
        """
        运行
        :return:
        """
        super(TransgenescanModule, self).run()
        self.run_transgenescan()

    def get_list(self):
        """
        根据list文件或者上传的文件名称获取样本名称
        list格式 file_name sample_name
        :return:
        """
        sample_dict = {}
        input_dir_path = self.option("in_fasta").prop['fasta_dir']
        self.logger.info("in_fasta: {}".format(input_dir_path))
        file_list = os.listdir(input_dir_path)
        if 'list.txt' in file_list: ##如果有list.txt文件，则从list.txt文件中直接获取样本与文件的对应关系
            list_path = os.path.join(input_dir_path, 'list.txt')
            with open(list_path, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if "file_name" in line:
                        pass
                    else:
                        line = line.strip().split("\t")
                        file_name = line[0]
                        sample_name = line[1]
                        if sample_name not in sample_dict.keys():
                            sample_dict[sample_name] = file_name
        else:##如果不存在list.txt文件，则从根据文件名称获取样本与文件的对应关系
            for file in file_list:
                if re.search(r'.fa', file):
                    file_name = file.strip(".fa")
                elif re.search(r'.fasta', file):
                    file_name = file.strip(".fasta")
                elif re.search(r'.ffn', file):
                    file_name = file.strip(".ffn")
                if file_name not in sample_dict.keys():
                    sample_dict[file_name] = file
        return sample_dict

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("开始将结果文件链接到结果文件下")
        predict = os.path.join(self.output_dir, "predict")
        if os.path.exists(predict):
            shutil.rmtree(predict)
        else:
            os.mkdir(predict)
        for tool in self.tool_list:
            link_dir(tool.output_dir, predict)
        self.logger.info("链接结果文件完成")
        self.end()

    def end(self):
        """
        运行结束
        :return:
        """
        result_dir = self.add_upload_dir(os.path.join(self.output_dir, 'predict'))
        self.logger.info("开始上传结果文件")
        result_dir.add_relpath_rules([
            [".", "", "TransGenScan结果文件目录"],
        ])
        super(TransgenescanModule, self).end()
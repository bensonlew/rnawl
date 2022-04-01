# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20190311
from biocluster.module import Module
import os, re
import shutil
from biocluster.core.exceptions import OptionError


class BaseInfoModule(Module):
    """
    用于统计多个fastq文件的碱基质量信息
    """
    def __init__(self, work_id):
        super(BaseInfoModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type":"infile", "format": "sequence.fastq_dir"},
        ]
        self.add_option(options)
        self.tools =[]

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("fastq_dir").is_set:
            raise OptionError("必须输入样本文件夹！")

    def run_base_info(self):
        """
        循环同时每10个样本投递
        :return:
        """
        fastq_path = self.option("fastq_dir").prop['path']
        """
        if not os.path.exists(self.work_dir + '/base_fastq'):
            os.mkdir(self.work_dir + '/base_fastq')
        else:
            pass
        base_path = os.path.join(self.work_dir, 'base_fastq')
        n = 1
        for file in os.listdir(fastq_path):
            if re.search(r'1.fq', file):
                name = str(file).split('.')[0]
                fq1_path = os.path.join(fastq_path, file)
                fq2_path = os.path.join(fastq_path, name + '.2.fq')
                if os.path.exists(os.path.join(base_path, "fastq_"+ str(n))):
                    pass
                else:
                    os.mkdir(os.path.join(base_path, "fastq_" + str(n)))
                new_file_path = os.path.join(base_path, "fastq_" + str(n))
                if len(os.listdir(new_file_path)) >= 10:
                    n += 1
                else:
                    os.link(fq1_path, os.path.join(new_file_path, file))
                    os.link(fq2_path, os.path.join(new_file_path, name + '.2.fq'))
        """
        for file in os.listdir(fastq_path):
            fastq_base_path = os.path.join(fastq_path, file)
            base_info = self.add_tool("meta.qc.mg_base_info")
            opts = ({
                "fastq_path": fastq_base_path,
            })
            base_info.set_options(opts)
            self.tools.append(base_info)
        if len(self.tools) >1:
            self.on_rely(self.tools, self.set_output)
        else:
            self.tools[0].on('end', self.set_output)
        for tool in self.tools:
            tool.run()

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        for i in self.tools:
            for f in os.listdir(i.output_dir):
                file_path = os.path.join(i.output_dir, f)
                new_path = os.path.join(self.output_dir, f)
                if os.path.exists(new_path):
                    os.remove(new_path)
                os.link(file_path, new_path)
        self.end()

    def run(self):
        super(BaseInfoModule, self).run()
        self.run_base_info()

    def end(self):
        super(BaseInfoModule, self).end()
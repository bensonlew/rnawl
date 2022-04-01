# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import shutil
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.toolapps.common import get_info,link_dir


class HumannModule(Module):
    """
    对每个文件质控后的reads进行分析
    """
    def __init__(self, work_id):
        super(HumannModule, self).__init__(work_id)
        options = [
            {"name": "fa_dir", "type": "infile", "format": "sequence.fasta_dir"},#指控后的reads序列的目录
            {"name": "search_mode", "type": "string", "default": "uniref50"},  # uniref50 or uniref90
            {"name": "prescreen_threshold", "type": 'float', "default": 0.01},
            {"name": "identity_threshold", "type": 'int', "default": 50},
            {"name": "subject_coverage", "type": 'int', "default": 50},
            {"name": "query_coverage", "type": 'int', "default": 90},
            {"name": "pathways", "type": "string", "default": "metacyc"},  # metacyc,unipathway
            {"name": "translated_alignment", "type": "string", "default": "diamond"},  # usearch,rapsearch,diamond
        ]
        self.modules = []
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fa_dir').is_set:
            raise OptionError('必须输入fa_dir文件夹')

    def run_humanns(self):
        self.sample_path = get_info(self.option('fa_dir').prop['path'])
        n=0
        for sample in sorted(self.sample_path.keys()):
            self.humann = self.add_tool('metagenomic.humann')
            opts = {
                "read1": self.sample_path[sample][0],
                "read2": self.sample_path[sample][1],
                "sample": sample,
                "search_mode": self.option("search_mode"),
                "prescreen_threshold": self.option("prescreen_threshold"),
                "identity_threshold": self.option("identity_threshold"),
                "subject_coverage": self.option("subject_coverage"),
                "query_coverage": self.option("query_coverage"),
                "pathways": self.option("pathways"),
                "translated_alignment": self.option("translated_alignment"),
            }
            self.humann.set_options(opts)
            self.modules.append(self.humann)
            self.humann.on('end', 'humann{}'.format(n))
            n += 1
        self.logger.info(self.modules)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.set_output)
        elif len(self.modules) == 1:
            self.modules[0].on("end", self.set_output)
        for module in self.modules:
            module.run()

    def run(self):
        """
        运行
        :return:
        """
        super(HumannModule, self).run()
        self.run_humanns()


    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if len(os.listdir(self.output_dir)) >= 1:
            for i in os.listdir(self.output_dir):
                os.remove(self.output_dir + "/" +i)
        if len(self.modules) > 1:
            for module in self.modules:
                file =os.listdir(module.output_dir)[0]
                link_dir(module.output_dir + "/" + file, self.output_dir + "/" + file)
        elif len(self.modules) == 1:
            file = os.listdir(self.modules[0].output_dir)[0]
            link_dir(self.modules[0].output_dir + "/" + file, self.output_dir + "/" + file)
        self.end()

    def end(self):
        super(HumannModule, self).end()
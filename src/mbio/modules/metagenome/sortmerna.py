# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import shutil
from collections import defaultdict
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.toolapps.common import link_dir


class SortmernaModule(Module):
    """
    对每个文件质控后的reads进行分析
    """
    def __init__(self, work_id):
        super(SortmernaModule, self).__init__(work_id)
        options = [
            {"name": "fa_dir", "type": "infile", "format": "sequence.fasta_dir"},#指控后的reads序列的目录
            {"name": "database", "type": 'string', "default": "rfam_5.8s,rfam_5s,arc_16s,arc_23s,bac_16s,bac_23s,euk_18s,euk_28s"},
        ]
        self.modules = []
        self.add_option(options)
        self.sample_path = defaultdict(list)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fa_dir').is_set:
            raise OptionError('必须输入fa_dir文件夹')

    def run_sortmerna(self):
        self.get_info()
        n=0
        for sample in self.sample_path:
            self.sortmerna = self.add_module("metagenome.single_sortmerna")
            opts = {
                "read1": self.sample_path[sample][0],
                "read2": self.sample_path[sample][1],
                "sample": sample,
                "database": self.option("database"),
                }
            self.sortmerna.set_options(opts)
            self.modules.append(self.sortmerna)
            self.sortmerna.on('end', 'sortmerna{}'.format(n))
            n += 1
        self.logger.info(self.modules)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.set_output)
        elif len(self.modules) == 1:
            self.modules[0].on("end", self.set_output)
        for module in self.modules:
            module.run()

    def get_info(self):
        with open(self.option('fa_dir').prop['path'] + "/list.txt") as fr:
            for line in fr:
                tmp = line.strip().split('\t')
                if tmp[1] in self.sample_path.keys():
                    if tmp[2] == 'l':
                        self.sample_path[tmp[1]].insert(0, self.option('fa_dir').prop['path']+ '/' + tmp[0])
                    else:
                        self.sample_path[tmp[1]].append(self.option('fa_dir').prop['path'] + '/' + tmp[0])
                else:
                    self.sample_path[tmp[1]].append(self.option('fa_dir').prop['path'] + '/' + tmp[0])

    def run(self):
        """
        运行
        :return:
        """
        super(SortmernaModule, self).run()
        self.run_sortmerna()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if len(os.listdir(self.output_dir)) >= 1:
            for i in os.listdir(self.output_dir):
                os.remove(self.output_dir + '/' + i)
        if len(self.modules) > 1:
            for module in self.modules:
                for i in os.listdir(module.output_dir):
                    os.link(module.output_dir + "/" + i, self.output_dir + "/" + i)
        else:
            link_dir(self.modules[0].output_dir, self.output_dir)
        self.end()

    def end(self):
        super(SortmernaModule, self).end()
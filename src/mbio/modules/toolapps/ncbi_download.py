# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.01.14

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir


class NcbiDownloadModule(Module):
    """
    从NCBI下载构建看家基因
    """
    def __init__(self, work_id):
        super(NcbiDownloadModule, self).__init__(work_id)
        options = [
            {"name": "sample_list", "type": "string"},  # 样品名称，sample1;sample2;...
            {"name": "s16", "type": "outfile", "format": "sequence.fasta"},
            {"name": "database_stat", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "cus_table", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "genomes", "type": "outfile", "format": "sequence.fasta_dir"},
        ]
        self.add_option(options)
        self.ncbi_download = self.add_tool('toolapps.ncbi_download')
        self.ncbi_stat = self.add_tool('toolapps.ncbi_stat')
        self.modules =[]


    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("sample_list"):
            raise OptionError("必须设置参数sample_list")

    def run_genome(self):
        opts = {
            'sample_list': self.option('sample_list'),
        }
        self.ncbi_download.set_options(opts)
        self.ncbi_download.on('end', self.run_housekeeping)
        self.ncbi_download.run()

    def run_housekeeping(self):
        for file in os.listdir(self.ncbi_download.option("genomes").prop['path']):
            sample = file.split(".fasta")[0]
            self.housekeeping = self.add_tool('toolapps.genome_housekeeping')
            opts = {
                'input_genome': self.ncbi_download.option("genomes").prop['path']+"/"+file,
                'sample': sample,
            }
            self.housekeeping.set_options(opts)
            self.modules.append(self.housekeeping)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.run_stat)
        elif len(self.modules) == 1:
            self.modules[0].on("end", self.run_stat)
        for module in self.modules:
            module.run()

    def run_stat(self):
        """
        合并16s和housekeeping结果
       :return:
        """
        if os.path.exists(self.work_dir+"/house"):
            shutil.rmtree(self.work_dir+"/house")
        os.mkdir(self.work_dir+"/house")
        for module in self.modules:
            link_dir(module.output_dir, self.work_dir+"/house")
        opts = {
            'sample_list': self.option('sample_list'),
            's16': self.ncbi_download.option("table").prop['path'],
            'house_dir': self.work_dir+"/house",
        }
        self.ncbi_stat.set_options(opts)
        self.ncbi_stat.on('end', self.set_output)
        self.ncbi_stat.run()

    def run(self):
        """
        运行
        :return:
        """
        super(NcbiDownloadModule, self).run()
        self.run_genome()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """

        self.option("genomes", self.ncbi_download.option("genomes").prop['path'])
        self.option("s16", self.ncbi_download.option("s16").prop['path'])
        self.option("database_stat", self.ncbi_stat.output_dir+"/all_database.stat.xls")
        self.option("cus_table", self.ncbi_download.option("table").prop['path'])
        self.end()

    def end(self):
        super(NcbiDownloadModule, self).end()
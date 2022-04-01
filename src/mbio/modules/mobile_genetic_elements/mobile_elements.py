# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import shutil
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir,link_file
from Bio import SeqIO

class MobileElementsModule(Module):
    """
    单个基因组可移动元件的预测.分五部分，IS、island、integron、prompter、transposon
    author: gaohao
    last_modify: 2020.07.23
    """
    def __init__(self, work_id):
        super(MobileElementsModule, self).__init__(work_id)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "gene_fna", "type": "infile",  "format": "sequence.fasta"},  # 预测的基因核酸文件
            {"name": "faa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "ptt", "type": "string"},  # 类似gff
            {"name": "rnt", "type": "string"},  # rNA的统计文件
            {"name": "sample", "type": "string"},
        ]
        self.island = self.add_module('mobile_genetic_elements.island')
        self.integron = self.add_module('mobile_genetic_elements.integron_predict')
        #self.prompter = self.add_tool('mobile_genetic_elements.prompter_predict')
        self.transposon = self.add_module('mobile_genetic_elements.transposon_predict')
        self.enzyme = self.add_tool('mobile_genetic_elements.diamond_enzyme')
        self.mobile_stat = self.add_tool('mobile_genetic_elements.mobile_stat')
        self.list = [self.transposon, self.island, self.integron, self.enzyme]
        self.step.add_steps('island', 'integron', 'prompter','transposon','enzyme', 'mobile_stat')
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("genome_fa").is_set:
            raise OptionError('必须输入genome_fa文件!')
        if not self.option("gene_fna").is_set:
            raise OptionError('必须输入gene_fna文件!')
        if not self.option("faa").is_set:
            raise OptionError('必须输入faa文件!')
        if not self.option("ptt"):
            raise OptionError('必须输入ptt文件!')
        if not self.option("rnt"):
            raise OptionError('必须输入rnt文件!')


    def run_transposon(self):
        opts = {
            "genome_fa": self.option("genome_fa"),
            "gene_faa": self.option("faa"),
            "gene_fna": self.option("gene_fna"),
            "gene_gff": self.option("ptt"),
            "sample": self.option("sample")
        }
        self.transposon.set_options(opts)
        self.transposon.on("end", self.set_output, "transposon")
        self.transposon.run()

    def run_island(self):
        opts = {
            "genome_fa": self.option("genome_fa"),
            "gene_fna": self.option("gene_fna"),
            "faa": self.option("faa"),
            "ptt": self.option("ptt"),
            "rnt": self.option("rnt"),
            "sample": self.option("sample")
        }
        self.island.set_options(opts)
        self.island.on("end", self.set_output, "island")
        self.island.run()

    def run_integron(self):
        opts = {
            "genome_fa": self.option("genome_fa"),
            "gene_faa": self.option("faa"),
            "gene_gff": self.option("ptt"),
        }
        self.integron.set_options(opts)
        self.integron.on("end", self.set_output, "integron")
        self.integron.run()

    """
    def run_prompter(self):
        opts = {
            "assemble": self.option("genome_fa"),
            "gff": self.option("ptt"),
            "sample": self.option("sample")
        }
        self.prompter.set_options(opts)
        self.prompter.on("end", self.set_output, "prompter")
        self.prompter.run()
    """

    def run_enzyme(self):
        opts = {
            "query":self.option("faa"),
            "ptt": self.option("ptt"),
        }
        self.enzyme.set_options(opts)
        self.enzyme.on("end",self.set_output, 'enzyme')
        self.enzyme.run()

    def run_mobile_stat(self):
        if os.path.exists(self.work_dir + "/files"):
            shutil.rmtree(self.work_dir + "/files")
        os.mkdir(self.work_dir + "/files")
        if os.path.exists(self.enzyme.output_dir + "/" + self.option("sample") + ".enzyme.xls"):
            link_file(self.enzyme.output_dir + "/" + self.option("sample") + ".enzyme.xls", self.work_dir + "/files/" + self.option("sample") + ".enzyme.xls")
        """
        if os.path.exists(self.prompter.output_dir + "/" + self.option("sample") + ".prompter.xls"):
            link_file(self.prompter.output_dir + "/" + self.option("sample") + ".prompter.xls", self.work_dir + "/files/" + self.option("sample") + ".prompter.xls")
        """
        if os.path.exists(self.integron.output_dir + "/" + self.option("sample") + ".integrons"):
            link_file(self.integron.output_dir + "/" + self.option("sample") + ".integrons", self.work_dir + "/files/" + self.option("sample") + ".integrons")
        if os.path.exists(self.island.output_dir + "/" + self.option("sample") + ".island.xls"):
            link_file(self.island.output_dir + "/" + self.option("sample") + ".island.xls", self.work_dir + "/files/" + self.option("sample") + ".island.xls")
        if os.path.exists(self.transposon.output_dir + "/" + self.option("sample") + ".transposon.xls"):
            link_file(self.transposon.output_dir + "/" + self.option("sample") + ".transposon.xls", self.work_dir + "/files/" + self.option("sample") + ".transposon.xls")
        if len(os.listdir(self.work_dir + "/files")) >=1:
            opts = {
                "genome": self.option("genome_fa"),
                "dir": self.work_dir + "/files",
                "sample_name": self.option("sample")
            }
            self.mobile_stat.set_options(opts)
            self.mobile_stat.on("end", self.set_output, 'mobile_stat')
            self.mobile_stat.run()
        else:
            self.end()


    def run(self):
        """
        运行
        :return:
        """
        super(MobileElementsModule, self).run()
        self.on_rely(self.list, self.run_mobile_stat)
        self.run_island()
        self.run_integron()
        #self.run_prompter()
        self.run_transposon()
        self.run_enzyme()

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        obj = event['bind_object']
        if event['data'] == 'island':
            if len(os.listdir(obj.output_dir)) >0:
                link_dir(obj.output_dir, self.output_dir + "/island")
        elif event['data'] == 'transposon':
            if len(os.listdir(obj.output_dir)) > 0:
                link_dir(obj.output_dir, self.output_dir + "/transposon")
        elif event['data'] == 'integron':
            if len(os.listdir(obj.output_dir)) > 0:
                link_dir(obj.output_dir, self.output_dir + "/integron")
        elif event['data'] == 'enzyme':
            if len(os.listdir(obj.output_dir)) > 0:
                link_dir(obj.output_dir, self.output_dir + "/enzyme")
        elif event['data'] == 'mobile_stat':
            if len(os.listdir(obj.output_dir)) > 0:
                link_dir(obj.output_dir, self.output_dir)
            self.end()

        """
        elif event['data'] == 'prompter':
            link_dir(obj.output_dir, self.output_dir + "/prompter")
        """



    def end(self):
        super(MobileElementsModule, self).end()
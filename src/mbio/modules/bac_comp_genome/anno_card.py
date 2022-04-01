# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify:20180916

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.bacgenome.common import sum_stat
from mbio.packages.bac_comp_genome.common_function import link_file

class AnnoCardModule(Module):
    def __init__(self, work_id):
        super(AnnoCardModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 100000},  # 将fasta序列拆分此行数的多个文件
            {"name": "sample", "type": "string"},  # 样品名称
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "card", "type": "outfile","format": "sequence.profile_table"}
            # 结果类型，aro_category or drug_class 细菌升级使用
        ]
        self.add_option(options)  #####检查option是否list格式，其中每个opt是否字典格式
        self.split_fasta = self.add_tool("sequence.split_fasta")
        self.card_align_tools = []
        self.card_anno_tool = self.add_tool("bac_comp_genome.card_anno")
        self.add_info_tools = []
        self.align_result_path = ''

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="21200301")
        if not self.option("sample"):
            raise OptionError("必须设置样品名称", code="21200302")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="21200303")
        return True

    def run_split_fasta(self):
        self.split_fasta.set_options({
            "fasta": self.option("query"),
            "lines": self.option("lines"),
        })
        self.split_fasta.on('end', self.card_align)
        self.split_fasta.run()

    def card_align(self):
        self.align_result_path = os.path.join(self.work_dir, "card_algin")
        if os.path.exists(self.align_result_path):
            pass
        else:
            os.mkdir(self.align_result_path)
        for f in os.listdir(self.split_fasta.output_dir):
            file_path = os.path.join(self.split_fasta.output_dir, f)
            align_diamond = self.add_tool('align.meta_diamond')
            align_diamond.set_options({
                "query": file_path,
                "database": "card",
                "evalue": self.option("evalue"),
                "sensitive": 0,
                "target_num": 1
            })
            self.card_align_tools.append(align_diamond)
        if len(self.card_align_tools) > 1:
            self.on_rely(self.card_align_tools, self.card_anno)
        else:
            self.card_align_tools[0].on('end', self.card_anno)
        for tool in self.card_align_tools:
            tool.run()

    def card_anno(self):
        xml_dir = self.align_result_path
        for i in self.card_align_tools:
            # print os.listdir(i.output_dir)
            for f in os.listdir(i.output_dir):
                if os.path.splitext(f)[1] == '.xml_new':
                    file_path = os.path.join(i.output_dir, f)
                    new_path = os.path.join(self.align_result_path, os.path.basename(file_path))
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(file_path, new_path)
        options = {
            "card_xml_dir": xml_dir,
        }
        self.card_anno_tool.set_options(options)
        self.card_anno_tool.on('end', self.set_output)
        self.card_anno_tool.run()

    def set_output(self):
        link_file(self.card_anno_tool.output_dir + "/gene_card_anno.xls", self.output_dir + "/" + self.option("sample") + ".card_anno.xls")
        self.option('card', self.output_dir + "/" + self.option("sample") + ".card_anno.xls")
        self.end()

    def run(self):
        super(AnnoCardModule, self).run()
        self.run_split_fasta()

    def end(self):
        super(AnnoCardModule, self).end()

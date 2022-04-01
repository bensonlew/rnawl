# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:2017.9.13

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class CardAnnotationModule(Module):
    def __init__(self, work_id):
        super(CardAnnotationModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 100000},  # 将fasta序列拆分此行数的多个文件
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},  # 基因丰度表
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "card_result_dir", "type": "outfile", "format": "annotation.mg_anno_dir"}  # 设置结果文件后面要用
        ]
        self.add_option(options)  #####检查option是否list格式，其中每个opt是否字典格式
        self.split_fasta = self.add_tool("sequence.split_fasta")
        # self.step.add_steps('split_fasta', 'diamond','link','anno','anno_stat')
        self.card_align_tools = []
        self.card_anno_tool = self.add_tool("annotation.card_anno")
        self.card_anno_stat_tool = self.add_tool("annotation.card_anno_stat")
        self.align_result_path = ''

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="21200201")
        if not self.option("reads_profile_table").is_set:
            raise OptionError("必须设置参数reads_profile_table", code="21200202")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="21200203")
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
                "database": "card_v3.0.9",## fix_by qingchen.zhang @20200811 更新数据库
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
        self.card_anno_tool.set_options({
            "card_xml_dir": xml_dir
        })
        self.card_anno_tool.on('end', self.card_anno_stat)
        self.card_anno_tool.run()

    def card_anno_stat(self):
        if self.card_anno_tool.option('card_anno_result').is_set:
            self.card_anno_stat_tool.set_options({
                'card_anno_table': self.card_anno_tool.option('card_anno_result'),
                'reads_profile_table': self.option('reads_profile_table'),
            })
            self.card_anno_stat_tool.on('end', self.set_output)
            self.card_anno_stat_tool.run()
        else:
            self.set_output()

    def set_output(self):
        self.option('card_result_dir', self.output_dir)
        anno_allfiles = os.listdir(self.card_anno_tool.output_dir)
        stat_allfiles = os.listdir(self.card_anno_stat_tool.output_dir)
        out_oldfiles = [os.path.join(self.card_anno_tool.output_dir, i) for i in anno_allfiles]
        stat_oldfiles = [os.path.join(self.card_anno_stat_tool.output_dir, i) for i in stat_allfiles]
        out_oldfiles.extend(stat_oldfiles)
        out_name = anno_allfiles + stat_allfiles
        output_newfiles = [os.path.join(self.output_dir, i) for i in out_name]
        for i in range(len(output_newfiles)):
            if os.path.exists(output_newfiles[i]):
                os.remove(output_newfiles[i])
            # print "old:",out_oldfiles[i]
            # print "new:",output_newfiles[i]
            os.link(out_oldfiles[i], output_newfiles[i])
        self.end()

    def run(self):
        super(CardAnnotationModule, self).run()
        self.run_split_fasta()

    def end(self):
        super(CardAnnotationModule, self).end()

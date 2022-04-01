# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'
# last_modify:2018.10.16

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class MvirdbAnnotationModule(Module):
    def __init__(self, work_id):
        super(MvirdbAnnotationModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 100000},  # 将fasta序列拆分此行数的多个文件
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},  # 基因丰度表
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "mvirdb_result_dir", "type": "outfile", "format": "annotation.mg_anno_dir"}  # 设置结果文件后面要用
        ]
        self.add_option(options)  ##### 检查option是否list格式，其中每个opt是否字典格式
        self.split_fasta = self.add_tool("sequence.split_fasta")
        self.mvirdb_align_tools = []
        self.mvirdb_anno_tool = self.add_tool("annotation.mvirdb_anno")

        self.mvirdb_anno_stat_tool = self.add_tool("annotation.mvirdb_anno_stat")
        self.align_result_path = ''

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="21202601")
        if not self.option("reads_profile_table").is_set:
            raise OptionError("必须设置参数reads_profile_table", code="21202602")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="21202603")
        return True

    def run_split_fasta(self):
        self.split_fasta.set_options({
            "fasta": self.option("query"),
            "lines": self.option("lines")
        })
        self.split_fasta.on('end', self.mvirdb_align)
        self.split_fasta.run()

    def mvirdb_align(self):
        self.align_result_path = os.path.join(self.work_dir, "mvirdb_algin")
        if os.path.exists(self.align_result_path):
            pass
        else:
            os.mkdir(self.align_result_path)
        for f in os.listdir(self.split_fasta.output_dir):
            file_path = os.path.join(self.split_fasta.output_dir, f)
            align_diamond = self.add_tool('align.meta_diamond')
            align_diamond.set_options({
                "query": file_path,
                "database": "mvirdb",
                "evalue": self.option("evalue"),
                "sensitive": 0,
                "target_num": 1
            })
            self.mvirdb_align_tools.append(align_diamond)
        if len(self.mvirdb_align_tools) > 1:
            self.on_rely(self.mvirdb_align_tools, self.mvirdb_anno)
        else:
            self.mvirdb_align_tools[0].on('end', self.mvirdb_anno)
        for tool in self.mvirdb_align_tools:
            tool.run()

    def mvirdb_anno(self):
        xml_dir = self.align_result_path
        for i in self.mvirdb_align_tools:
            # print os.listdir(i.output_dir)
            for f in os.listdir(i.output_dir):
                if os.path.splitext(f)[1] == '.xml_new':
                    file_path = os.path.join(i.output_dir, f)
                    new_path = os.path.join(self.align_result_path, os.path.basename(file_path))
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(file_path, new_path)  #####创建硬链
        # xml_dir = self.mvirdb_align.output_dir.prop['path']
        self.mvirdb_anno_tool.set_options({
            "mvirdb_xml_dir": xml_dir
        })
        self.mvirdb_anno_tool.on('end', self.mvirdb_anno_stat)
        self.mvirdb_anno_tool.run()

    def mvirdb_anno_stat(self):
        self.mvirdb_anno_stat_tool.set_options({
            'mvirdb_anno_table': self.mvirdb_anno_tool.option('mvirdb_anno_result'),
            'reads_profile_table': self.option('reads_profile_table'),
            'is_modify' : 'T'
        })
        self.mvirdb_anno_stat_tool.on('end', self.set_output)
        self.mvirdb_anno_stat_tool.run()

    def set_output(self):
        self.option('mvirdb_result_dir', self.output_dir)
        anno_allfiles = os.listdir(self.mvirdb_anno_tool.output_dir)
        stat_allfiles = os.listdir(self.mvirdb_anno_stat_tool.output_dir)
        out_oldfiles = [os.path.join(self.mvirdb_anno_tool.output_dir, i) for i in anno_allfiles]
        stat_oldfiles = [os.path.join(self.mvirdb_anno_stat_tool.output_dir, i) for i in stat_allfiles]
        out_oldfiles.extend(stat_oldfiles)
        out_name = anno_allfiles + stat_allfiles
        output_newfiles = [os.path.join(self.output_dir, i) for i in out_name]
        for i in range(len(output_newfiles)):
            if os.path.exists(output_newfiles[i]):
                os.remove(output_newfiles[i])
            # print "old:",out_oldfiles[i]
            # print "new:",output_newfiles[i]
            os.link(out_oldfiles[i], output_newfiles[i])
        if os.path.exists(self.output_dir + '/gene_mvirdb_anno.xls'):
            os.remove(self.output_dir + '/gene_mvirdb_anno.xls')
        os.link(self.mvirdb_anno_stat_tool.work_dir + '/mvirdb_anno_rm_1.xls',self.output_dir + '/gene_mvirdb_anno.xls')

        self.end()

    def run(self):
        super(MvirdbAnnotationModule, self).run()
        self.run_split_fasta()

    def end(self):
        super(MvirdbAnnotationModule, self).end()

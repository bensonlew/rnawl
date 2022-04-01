# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:2017.9.29

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class VfdbDnaAnnoModule(Module):
    def __init__(self, work_id):
        super(VfdbDnaAnnoModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 100000},  # 将fasta序列拆分此行数的多个文件
            {"name": "sample", "type": "string"},  # 样品名称
            {"name": "evalue", "type": "float", "default": 1e-5}  # evalue值
        ]
        self.add_option(options)
        self.add_info_tools = []
        self.split_fasta = self.add_tool("sequence.split_fasta")
        self.vfdb_core_anno = self.add_tool("annotation.vfdb_anno")
        self.vfdb_align_core_tools = []
        self.vfdb_align_predict_tools = []
        self.vfdb_anno_tools = []
        self.vfdb_anno_tool_core = ''
        self.vfdb_anno_tool_pre = ''
        self.align_core_path = ''
        self.align_predict_path = ''

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="21201801")
        if not self.option("sample"):
            raise OptionError("必须设置样品名称", code="21201802")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="21201803")
        return True

    ######  分割query fasta序列
    def run_split_fasta(self):
        self.split_fasta.set_options({
            "fasta": self.option("query"),
            "lines": self.option("lines"),
        })
        self.split_fasta.on('end', self.vfdb_align_core)
        self.split_fasta.run()

    #### 序列先比对核心库，输出没比对上的fasta序列进行predict 库比对
    def vfdb_align_core(self):
        self.align_core_path = os.path.join(self.work_dir, "vfdb_algin_core")
        for file in os.listdir(self.split_fasta.output_dir):
            file_path = os.path.join(self.split_fasta.output_dir, file)
            core_align_tool = self.add_tool('align.meta_diamond')
            core_align_tool.set_options({
                "query": file_path,
                "database": "vfdb_v20200703_core",## fix_by qingchen.zhang @20200811  数据库更新版本为vfdb_20200703
                "evalue": self.option("evalue"),
                "sensitive": 0,
                "target_num": 1,
                "unalign": 1
            })
            self.vfdb_align_core_tools.append(core_align_tool)
            self.logger.info("run core align tool")
        if len(self.vfdb_align_core_tools) > 1:
            self.on_rely(self.vfdb_align_core_tools, self.vfdb_anno)
        else:
            self.vfdb_align_core_tools[0].on('end', self.vfdb_anno)
        for tool in self.vfdb_align_core_tools:
            tool.run()

    #### 比对全部完成后，调用vfdb_anno工具根据输入参数不同对核心库和预测库
    def vfdb_anno(self):
        core_dir = self.align_core_path
        self.link(self.vfdb_align_core_tools, core_dir)
        self.vfdb_core_anno.set_options({
            "vfdb_xml_dir": core_dir,
            "database": "core",
            "result_format":"dna"
        })
        self.vfdb_core_anno.on('end', self.change_table)
        self.vfdb_core_anno.run()

    def change_table(self):
        for file in ["gene_vfdb_core_anno.xls", "vfdb_level.xls"]:
            add_info = self.add_tool("bacgenome.add_gene_info")
            add_info.set_options({
                'sample': self.option('sample'),
                'sequence': self.option('query'),
                'table': self.vfdb_core_anno.output_dir + '/' + file
            })
            self.add_info_tools.append(add_info)
        if len(self.add_info_tools) > 1:
            self.on_rely(self.add_info_tools, self.set_output)
        else:
            self.add_info_tools[0].on('end', self.set_output)
        for tool in self.add_info_tools:
            tool.run()

    def link(self, align_result, xml_dir):
        if os.path.exists(xml_dir):
            pass
        else:
            os.mkdir(xml_dir)
        for i in align_result:
            # print os.listdir(i.output_dir)
            for f in os.listdir(i.output_dir):
                if os.path.splitext(f)[1] == '.xml_new':
                    file_path = os.path.join(i.output_dir, f)
                    #new_path = os.path.join(xml_dir, os.path.basename(file_path))
                    new_path = os.path.join(xml_dir, f)
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(file_path, new_path)  #####创建硬链

    def set_output(self):
        self.link_file()
        self.end()

    def link_file(self):
        """
        link文件到本module的output目录
        """
        allfiles = ["vfdb_core_align_table.xls","gene_vfdb_core_anno.xls","vfdb_level.xls"]
        newfile = [self.option('sample')+"_vfdb_align.xls",self.option('sample')+"_vfdb_anno.xls",self.option('sample')+"_vfdb_level.xls"]
        oldfiles = [os.path.join(self.vfdb_core_anno.output_dir, i) for i in allfiles]
        newfiles = [os.path.join(self.output_dir, i) for i in newfile]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def run(self):
        super(VfdbDnaAnnoModule, self).run()
        self.run_split_fasta()

    def end(self):
        super(VfdbDnaAnnoModule, self).end()

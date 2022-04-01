# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify:2019.0911

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_file,vfdb_anno


class AnnoVfdbModule(Module):
    def __init__(self, work_id):
        super(AnnoVfdbModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 100000},  # 将fasta序列拆分此行数的多个文件
            {"name": "sample", "type": "string"},  # 样品名称
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "vfdb", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)
        self.vfdb_align_core_tools = []
        self.vfdb_align_predict_tools = []
        self.vfdb_anno_tools = []
        self.vfdb_anno_tool_core = ''
        self.vfdb_anno_tool_pre = ''
        self.split_fasta = self.add_tool("sequence.split_fasta")
        self.align_core_path = ''
        self.align_predict_path = ''

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数")
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
                "database": "vfdb_core",
                "evalue": self.option("evalue"),
                "sensitive": 0,
                "target_num": 1,
                "unalign": 1
            })
            self.vfdb_align_core_tools.append(core_align_tool)
            self.logger.info("run core align tool")
        if len(self.vfdb_align_core_tools) > 1:
            self.on_rely(self.vfdb_align_core_tools, self.predict_align)
        else:
            self.vfdb_align_core_tools[0].on('end', self.predict_align)
        for tool in self.vfdb_align_core_tools:
            tool.run()

    def predict_align(self):
        self.align_predict_path = os.path.join(self.work_dir, "vfdb_algin_predict")
        for tool in self.vfdb_align_core_tools:
            files = os.listdir(tool.output_dir)
            for file in files:
                self.logger.info(file)
                if "_unalign.fasta" in file:
                    infile = os.path.join(tool.output_dir, file)
                    self.logger.info(infile)
                    self.logger.info("add " + file + " to align")
                else:
                    self.logger.info("the file is not unalign file!")
            predict_tool = self.add_tool('align.meta_diamond')
            predict_tool.set_options({
                "query": infile,
                "database": "vfdb_predict",
                "evalue": self.option("evalue"),
                "sensitive": 0,
                "target_num": 1,
            })
            self.vfdb_align_predict_tools.append(predict_tool)
        if len(self.vfdb_align_predict_tools) > 1:
            self.on_rely(self.vfdb_align_predict_tools, self.vfdb_anno)
        else:
            self.vfdb_align_predict_tools[0].on('end', self.vfdb_anno)
        for tool in self.vfdb_align_predict_tools:
            tool.run()

    #### 比对全部完成后，调用vfdb_anno工具根据输入参数不同对核心库和预测库
    def vfdb_anno(self):
        core_dir = self.align_core_path
        self.link(self.vfdb_align_core_tools, core_dir)
        vfdb_core_anno = self.add_tool("bac_comp_genome.anno_vfdb")
        vfdb_core_anno.set_options({
            "vfdb_xml_dir": core_dir,
            "database": "core"
        })
        self.vfdb_anno_tools.append(vfdb_core_anno)
        self.vfdb_anno_tool_core = vfdb_core_anno
        pre_dir = self.align_predict_path
        self.link(self.vfdb_align_predict_tools, pre_dir)
        vfdb_pre_anno = self.add_tool("bac_comp_genome.anno_vfdb")
        vfdb_pre_anno.set_options({
            "vfdb_xml_dir": pre_dir,
            "database": "predict"
        })
        self.vfdb_anno_tools.append(vfdb_pre_anno)
        self.vfdb_anno_tool_pre = vfdb_pre_anno
        self.on_rely(self.vfdb_anno_tools, self.set_output)
        for tool in self.vfdb_anno_tools:
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
                    new_path = os.path.join(xml_dir, f)
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(file_path, new_path)  #####创建硬链

    def set_output(self):
        for tool in self.vfdb_anno_tools:
            if os.path.exists(tool.output_dir + "/gene_vfdb_core_anno.xls"):
                link_file(tool.output_dir + "/gene_vfdb_core_anno.xls", self.output_dir + "/" +self.option("sample") + ".vfdb_core_anno.xls")
            elif os.path.exists(tool.output_dir + "/gene_vfdb_predict_anno.xls"):
                link_file(tool.output_dir + "/gene_vfdb_predict_anno.xls", self.output_dir + "/" +self.option("sample") + ".vfdb_predict_anno.xls")
        vfdb_anno(self.output_dir + "/" +self.option("sample") + ".vfdb_core_anno.xls", self.output_dir + "/" +self.option("sample") + ".vfdb_predict_anno.xls", self.output_dir + "/" +self.option("sample") + ".vfdb_anno.xls")
        self.option("vfdb", self.output_dir + "/" +self.option("sample") + ".vfdb_anno.xls")
        self.end()

    def run(self):
        super(AnnoVfdbModule, self).run()
        self.run_split_fasta()

    def end(self):
        super(AnnoVfdbModule, self).end()
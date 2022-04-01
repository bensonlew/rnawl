# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
# last_modify:20181026

from biocluster.module import Module
#from mbio.modules.annotation.phi_dna_anno import PhiDnaAnnoModule
import os
import shutil
from biocluster.core.exceptions import OptionError

class PhiAnnotationModule(Module):
    def __init__(self, work_id):
        Module.__init__(self, work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 100000},  # 将fasta序列拆分此行数的多个文件
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},  # 基因丰度表
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "phi_result_dir", "type": "outfile", "format": "annotation.mg_anno_dir"}  # 设置结果文件后面要用
        ]
        self.add_option(options)
        self.logger.info(options)
        self.split_fasta = self.add_tool("sequence.split_fasta")
        self.phi_align_tools = []
        self.phi_anno_tool = self.add_tool("annotation.phi_anno")
        self.phi_anno_stat_tool = self.add_tool("annotation.phi_anno_stat")
        self.align_result_path = ''

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="21202501")
        if not self.option("reads_profile_table").is_set:
            raise OptionError("必须设置参数reads_profile_table", code="21202502")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="21202503")
        return True

    def run_split_fasta(self):
        self.split_fasta.set_options({
            "fasta": self.option("query"),
            "lines": self.option("lines"),
        })
        self.split_fasta.on('end', self.phi_align)
        self.split_fasta.run()

    def phi_align(self):
        self.align_result_path = os.path.join(self.work_dir, "phi_algin")
        if os.path.exists(self.align_result_path):
            pass
        else:
            os.mkdir(self.align_result_path)
        for f in os.listdir(self.split_fasta.output_dir):
            file_path = os.path.join(self.split_fasta.output_dir, f)
            align_diamond = self.add_tool('align.meta_diamond')
            align_diamond.set_options({
                "query": file_path,
                "database": "phi_v49",
                "evalue": self.option("evalue"),
                "sensitive": 0,
                "target_num": 1
            })
            self.phi_align_tools.append(align_diamond)
        if len(self.phi_align_tools) > 1:
            self.on_rely(self.phi_align_tools, self.phi_anno)
        else:
            self.phi_align_tools[0].on('end', self.phi_anno)
        for tool in self.phi_align_tools:
            tool.run()

    def phi_anno(self):
        xml_dir = self.align_result_path
        for i in self.phi_align_tools:
            # print os.listdir(i.output_dir)
            for f in os.listdir(i.output_dir):
                if os.path.splitext(f)[1] == '.xml_new':
                    file_path = os.path.join(i.output_dir, f)
                    new_path = os.path.join(self.align_result_path, os.path.basename(file_path))
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(file_path, new_path)
        self.phi_anno_tool.set_options({
            "phi_xml_dir": xml_dir,
            "project": "meta"
            # "result_format":"dna"
        })
        self.phi_anno_tool.on('end', self.phi_anno_stat)
        self.phi_anno_tool.run()

    def phi_anno_stat(self):
        self.phi_anno_stat_tool.set_options({
            'phi_anno_table': self.phi_anno_tool.option('phi_anno_result'),
            'reads_profile_table': self.option('reads_profile_table'),
        })
        self.phi_anno_stat_tool.on('end', self.set_output)
        self.phi_anno_stat_tool.run()

    def link_file(self):
        self.option('phi_result_dir', self.output_dir)
        # anno_allfiles = os.listdir(self.phi_anno_tool.output_dir)
        anno_allfiles = ["gene_phi_anno.xls"]
        stat_allfiles = os.listdir(self.phi_anno_stat_tool.output_dir)
        out_oldfiles = [os.path.join(self.phi_anno_tool.output_dir, i) for i in anno_allfiles]
        stat_oldfiles = [os.path.join(self.phi_anno_stat_tool.output_dir, i) for i in stat_allfiles]
        out_oldfiles.extend(stat_oldfiles)
        out_name = anno_allfiles + stat_allfiles
        output_newfiles = [os.path.join(self.output_dir, i) for i in out_name]
        for i in range(len(output_newfiles)):
            if os.path.exists(output_newfiles[i]):
                os.remove(output_newfiles[i])
            # print "old:",out_oldfiles[i]
            # print "new:",output_newfiles[i]
            os.link(out_oldfiles[i], output_newfiles[i])

    def set_output(self):
        self.link_file()
        self.end()

    def run(self):
        super(PhiAnnotationModule, self).run()
        self.run_split_fasta()

    def end(self):
        super(PhiAnnotationModule, self).end()


# 所继承的父类代码
class PhiDnaAnnoModule(Module):
    def __init__(self, work_id):
        super(PhiDnaAnnoModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 100000},  # 将fasta序列拆分此行数的多个文件
            {"name": "sample", "type": "string"},  # 样品名称
            {"name": "evalue", "type": "float", "default": 1e-5}  # evalue值
        ]
        self.add_option(options)  #####检查option是否list格式，其中每个opt是否字典格式
        self.split_fasta = self.add_tool("sequence.split_fasta")
        self.add_info_tools = []
        # self.step.add_steps('split_fasta', 'diamond','link','anno','anno_stat')
        self.phi_align_tools = []
        self.phi_anno_tool = self.add_tool("annotation.phi_anno")
        # self.phi_anno_stat_tool = self.add_tool("annotation.phi_anno_stat")
        self.align_result_path = ''

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="21202504")
        if not self.option("sample"):
            raise OptionError("必须设置样品名称", code="21202505")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="21202506")
        return True

    def run_split_fasta(self):
        self.split_fasta.set_options({
            "fasta": self.option("query"),
            "lines": self.option("lines"),
        })
        self.split_fasta.on('end', self.phi_align)
        self.split_fasta.run()

    def phi_align(self):
        self.align_result_path = os.path.join(self.work_dir, "phi_algin")
        if os.path.exists(self.align_result_path):
            pass
        else:
            os.mkdir(self.align_result_path)
        for f in os.listdir(self.split_fasta.output_dir):
            file_path = os.path.join(self.split_fasta.output_dir, f)
            align_diamond = self.add_tool('align.meta_diamond')
            align_diamond.set_options({
                "query": file_path,
                "database": "phi",
                "evalue": self.option("evalue"),
                "sensitive": 0,
                "target_num": 1
            })
            self.phi_align_tools.append(align_diamond)
        if len(self.phi_align_tools) > 1:
            self.on_rely(self.phi_align_tools, self.phi_anno)
        else:
            self.phi_align_tools[0].on('end', self.phi_anno)
        for tool in self.phi_align_tools:
            tool.run()

    def phi_anno(self):
        xml_dir = self.align_result_path
        for i in self.phi_align_tools:
            # print os.listdir(i.output_dir)
            for f in os.listdir(i.output_dir):
                if os.path.splitext(f)[1] == '.xml_new':
                    file_path = os.path.join(i.output_dir, f)
                    new_path = os.path.join(self.align_result_path, os.path.basename(file_path))
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(file_path, new_path)
        self.phi_anno_tool.set_options({
            "phi_xml_dir": xml_dir
            # "result_format":"dna"
        })
        self.phi_anno_tool.on('end', self.change_table)
        self.phi_anno_tool.run()

    def change_table(self):
        for file in ["gene_phi_anno.xls", "phi_stat.xls", "phi_phenotype.xls"]:
            add_info = self.add_tool("bacgenome.add_gene_info")
            add_info.set_options({
                'sample': self.option('sample'),
                'sequence': self.option('query'),
                'table': self.phi_anno_tool.output_dir + '/' + file
            })
            self.add_info_tools.append(add_info)
        if len(self.add_info_tools) > 1:
            self.on_rely(self.add_info_tools, self.set_output)
        else:
            self.add_info_tools[0].on('end', self.set_output)
        for tool in self.add_info_tools:
            tool.run()

    def set_output(self):
        self.run_link_file()
        self.end()

    def run_link_file(self):
        """
        link文件到本module的output目录
        """
        allfiles = ["phi_align_table.xls", "gene_phi_anno.xls", "phi_stat.xls","phi_phenotype.xls"]
        newfile = [self.option('sample') + "_phi_align.xls", self.option('sample') + "_phi_anno.xls",
                   self.option('sample') + "_phi_stat.xls",self.option('sample') + "_phi_phenotype.xls"]
        oldfiles = [os.path.join(self.phi_anno_tool.output_dir, i) for i in allfiles]
        newfiles = [os.path.join(self.output_dir, i) for i in newfile]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def run(self):
        super(PhiDnaAnnoModule, self).run()
        self.run_split_fasta()

    def end(self):
        super(PhiDnaAnnoModule, self).end()

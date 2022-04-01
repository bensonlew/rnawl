# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20180213

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
import glob


class TcdbDnaAnnoModule(Module):
    def __init__(self, work_id):
        super(TcdbDnaAnnoModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 100000},  # 将fasta序列拆分此行数的多个文件
            {"name": "sample", "type": "string"},  # 样品名称
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "analysis", "type": "string", "default":""}, #complete or  uncomplete
            {"name": "has_coverage", "type": "bool", 'default': False}
        ]
        self.add_option(options)  #####检查option是否list格式，其中每个opt是否字典格式
        self.split_fasta = self.add_tool("sequence.split_fasta")
        self.add_info_tools = []
        # self.step.add_steps('split_fasta', 'diamond','link','anno','anno_stat')
        self.tcdb_align_tools = []
        self.tcdb_anno_tool = self.add_tool("annotation.tcdb_anno")
        # self.tcdb_anno_stat_tool = self.add_tool("annotation.tcdb_anno_stat")
        self.align_result_path = ''

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="21201601")
        if not self.option("sample"):
            raise OptionError("必须设置样品名称", code="21201602")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="21201603")
        return True

    def run_split_fasta(self):
        self.split_fasta.set_options({
            "fasta": self.option("query"),
            "lines": self.option("lines"),
        })
        self.split_fasta.on('end', self.tcdb_align)
        self.split_fasta.run()

    def tcdb_align(self):
        self.align_result_path = os.path.join(self.work_dir, "tcdb_algin")
        if os.path.exists(self.align_result_path):
            pass
        else:
            os.mkdir(self.align_result_path)
        for f in os.listdir(self.split_fasta.output_dir):
            file_path = os.path.join(self.split_fasta.output_dir, f)
            align_diamond = self.add_tool('align.meta_diamond')
            align_diamond.set_options({
                "query": file_path,
                "database": "tcdb_v20200917",
                "evalue": self.option("evalue"),
                "sensitive": 0,
                "target_num": 5
            })
            self.tcdb_align_tools.append(align_diamond)
        if len(self.tcdb_align_tools) > 1:
            self.on_rely(self.tcdb_align_tools, self.tcdb_anno)
        else:
            self.tcdb_align_tools[0].on('end', self.tcdb_anno)
        for tool in self.tcdb_align_tools:
            tool.run()

    def tcdb_anno(self):
        self.logger.info(">>>in func tcdb_anno")
        xml_dir = self.align_result_path
        self.logger.info(">>>check tcdbaligntools")
        for i in self.tcdb_align_tools:
            # print os.listdir(i.output_dir)
            for f in os.listdir(i.output_dir):
                if os.path.splitext(f)[1] == '.xml_new':
                    file_path = os.path.join(i.output_dir, f)
                    new_path = os.path.join(self.align_result_path, os.path.basename(file_path))
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(file_path, new_path)
        self.logger.info(">>>start set option")
        self.tcdb_anno_tool.set_options({
            "tcdb_xml_dir": xml_dir,
            "has_coverage" : self.option('has_coverage')
            # "result_format":"dna"
        })
        self.tcdb_anno_tool.on('end', self.change_table)
        self.tcdb_anno_tool.run()
        self.logger.info(">>>tcdb_anno is start running")

    def change_table(self):
        self.logger.info(">>>change_table is start running")
        for file in ["gene_tcdb_anno.xls"]:
            table_file = self.tcdb_anno_tool.output_dir + '/' + file
            add_info = self.add_tool("bacgenome.add_gene_info")
            opts = {
                'sample': self.option('sample'),
                'sequence': self.option('query'),
                'table': table_file,
                #'split': True
            }
            if self.option("analysis") in ["complete"]:
                opts["split"] = True
            else:
                opts["split"] = False
            add_info.set_options(opts)
            self.add_info_tools.append(add_info)
        if len(self.add_info_tools) > 1:
            self.on_rely(self.add_info_tools, self.set_output)
        else:
            self.add_info_tools[0].on('end', self.set_output)
        for tool in self.add_info_tools:
            tool.run()

    def set_output(self):
        self.linkdir(self.tcdb_anno_tool.output_dir, "")
        self.rename("gene_tcdb_anno.xls", self.option('sample') + "_tcdb_anno.xls")
        self.rename("tcdb_align_table.xls", self.option('sample') + "_tcdb_align.xls")
        self.rename("top1_tcdb_align_table.xls", self.option('sample') + "_tcdb_align_top1.xls")

        ###避免重运行时因为gene_tcdb_anno.xls不存在而报错
        table_files = glob.glob(self.tcdb_anno_tool.output_dir+'/*_whole_genome_tcdb_anno.xls')
        if len(table_files) != 0:
            table_file = table_files[0]
            tcdb_tmp = self.tcdb_anno_tool.output_dir+'/gene_tcdb_anno.xls'
            if not os.path.exists(tcdb_tmp):
                os.link(table_file,tcdb_tmp)

        self.end()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def rename(self, old, new):
        if os.path.exists(self.output_dir + '/' + old):
            os.rename(self.output_dir + '/' + old, self.output_dir + '/' + new)

    def run(self):
        super(TcdbDnaAnnoModule, self).run()
        self.run_split_fasta()

    def end(self):
        super(TcdbDnaAnnoModule, self).end()

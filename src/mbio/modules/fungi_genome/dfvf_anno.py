# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modify:2018.5.26

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class DfvfAnnoModule(Module):
    def __init__(self, work_id):
        super(DfvfAnnoModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 100000},  # 将fasta序列拆分此行数的多个文件
            {"name": "sample", "type": "string"},  # 样品名称
            {"name": "evalue", "type": "float", "default": 1e-5} , # evalue值
            {"name": "database","type":"string", "default": "dfvf"}
        ]
        self.add_option(options)
        self.split_fasta = self.add_tool("sequence.split_fasta")
        self.dfvf_add_info = self.add_tool("fungi_genome.dfvf_add_info")
        self.dfvf_align_tools = []


    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="22100401")
        if not self.option("sample"):
            raise OptionError("必须设置样品名称", code="22100402")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="22100403")
        return True

    ######  分割query fasta序列
    def run_split_fasta(self):
        self.split_fasta.set_options({
            "fasta": self.option("query"),
            "lines": self.option("lines"),
        })
        self.split_fasta.on('end', self.run_dfvf_align)
        self.split_fasta.run()

    #### 序列先比对核心库，输出没比对上的fasta序列进行predict 库比对
    # def vfdb_align_core(self):
    def run_dfvf_align(self):
        self.align_path = os.path.join(self.work_dir, "dfvf_algin_result")
        for file in os.listdir(self.split_fasta.output_dir):
            file_path = os.path.join(self.split_fasta.output_dir, file)
            dfvf_align_tool = self.add_tool('align.meta_diamond')
            dfvf_align_tool.set_options({
                "query": file_path,
                "database": self.option("database"),
                "evalue": self.option("evalue"),
                "sensitive": 0,
                "target_num": 1,
                "unalign": 1,
                "outfmt": 6
            })
            self.dfvf_align_tools.append(dfvf_align_tool)
            self.logger.info("run diamond align tool")
        if len(self.dfvf_align_tools) > 1:
            self.on_rely(self.dfvf_align_tools, self.run_dfvf_add_info)
        else:
            self.dfvf_align_tools[0].on('end', self.run_dfvf_add_info)
        for tool in self.dfvf_align_tools:
            tool.run()

    def run_dfvf_add_info(self):
        self.link(self.dfvf_align_tools, self.align_path)
        splitfiles = os.path.join(self.align_path,"*.xls")
        catname = "{}_{}.xls".format(self.option("sample"), self.option("database"))
        catfile = os.path.join(self.align_path,catname)
        os.system("cat {} > {} ".format(splitfiles, catfile))
        database_path = self.path.join(self.config.SOFTWARE_DIR,"database/dfvf/dfvf.data")
        self.dfvf_add_info.set_options({
            "bsnout": catfile,
            "database": database_path
        })
        self.dfvf_add_info.on("end",self.set_output)
        self.dfvf_add_info.run()



    def link(self, align_result, new_dir):
        if os.path.exists(new_dir):
            pass
        else:
            os.mkdir(new_dir)
        for i in align_result:
            # print os.listdir(i.output_dir)
            for f in os.listdir(i.output_dir):
                # if os.path.splitext(f)[1] == '.xml_new':
                if os.path.splitext(f)[1] == 'xls' :
                    file_path = os.path.join(i.output_dir, f)
                    #new_path = os.path.join(xml_dir, os.path.basename(file_path))
                    new_path = os.path.join(new_dir, f)
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
        #allfiles = ["vfdb_core_align_table.xls","gene_vfdb_core_anno.xls","vfdb_level.xls"]
        allfiles = ["result.out.xls"]
        #newfile = [self.option('sample')+"_vfdb_align.xls",self.option('sample')+"_vfdb_anno.xls",self.option('sample')+"_vfdb_level.xls"]
        newfiles_all = ["{}_{}.result.xls".format(self.option('sample'), self.option("database"))]
        oldfiles = [os.path.join(self.dfvf_add_info.output_dir, i) for i in allfiles]
        newfiles = [os.path.join(self.output_dir, i) for i in newfiles_all]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def run(self):
        super(DfvfAnnoModule, self).run()
        self.run_split_fasta()

    def end(self):
        super(DfvfAnnoModule, self).end()

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
            {"name": "database","type":"string", "default": "DFVF"}
        ]
        self.add_option(options)
        self.split_fasta = self.add_tool("sequence.split_fasta")
        self.dfvf_add_info = self.add_tool("annotation.dfvf_add_info")
        self.dfvf_align_tools = []
        self.catfile = ""

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="21200801")
        if not self.option("sample"):
            raise OptionError("必须设置样品名称", code="21200802")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="21200803")
        return True

    ######  分割query fasta序列
    def run_split_fasta(self,faa):
        self.split_fasta.set_options({
            #"fasta": self.option("query"),
            "fasta": faa,
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
        self.logger.info("开始运行run_dfvf_add_info")
        new_path_str = self.link(self.dfvf_align_tools, self.align_path)
        self.logger.info(new_path_str)
        #catname = "align_table.xls".format(self.option("sample"), self.option("database").lower())
        catname = "align_table.xls"
        self.catfile = os.path.join(self.align_path,catname)
        os.system("cat {} > {} ".format(new_path_str, self.catfile))
        self.dfvf_add_info.set_options({
            "bsnout": self.catfile
        })
        self.dfvf_add_info.on("end",self.set_output)
        self.logger.info("开始运行dfvf_add_info")
        self.dfvf_add_info.run()



    def link(self, align_result, new_dir):
        if os.path.exists(new_dir):
            pass
        else:
            os.mkdir(new_dir)
        new_path_str = ""
        for i in align_result:
            for f in os.listdir(i.output_dir):
                if f.split('.')[1] == "xls":  # "outfmt": 6 时生成的比对结果文件后缀为 .xls
                    file_path = os.path.join(i.output_dir, f)
                    #new_path = os.path.join(xml_dir, os.path.basename(file_path))
                    new_path = os.path.join(new_dir, f)
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(file_path, new_path)  #####创建硬链
                    new_path_str += new_path + " "
        return new_path_str

    def set_output(self):
        self.link_file()
        self.end()

    def link_file(self):
        """
        link文件到本module的output目录
        """
        allfiles = ["result.out.xls"]
        #newfiles_all = ["{}_{}.result.xls".format(self.option('sample'), self.option("database"))]
        newfiles_all = ["{}_dfvf.anno.xls".format(self.option('sample'))]
        oldfiles = [os.path.join(self.dfvf_add_info.output_dir, i) for i in allfiles]
        oldfiles.append(self.work_dir + '/dfvf_algin_result/align_table.xls')
        self.logger.info(self.work_dir + '/dfvf_algin_result/align_table.xls')
        # if not os.path.exists(self.output_dir+"/DFVF"):
        #     os.mkdir(self.output_dir+"/DFVF")
        newfiles = [os.path.join(self.output_dir, i) for i in newfiles_all]
        oldfiles.append(self.catfile)
        newfiles.append(os.path.join(self.output_dir,"{}_{}_align_table.xls".format(self.option("sample"), self.option("database").lower())))
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(newfiles)):
            self.logger.info(oldfiles[i])
            os.link(oldfiles[i], newfiles[i])

    def run(self):
        super(DfvfAnnoModule, self).run()
        os.system("sed 's/ /__/g' {} > {}/new.faa".format(self.option('query').prop['path'],self.work_dir))
        self.run_split_fasta('{}/new.faa'.format(self.work_dir))

    def end(self):
        super(DfvfAnnoModule, self).end()


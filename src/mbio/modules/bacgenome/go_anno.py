# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.module import Module
import os,re
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table


class GoAnnoModule(Module):
    """
    细菌、真菌针对go数据库升级改用，主要是nr注释的结果变动导致后面注释存在问题
    与宏基因组区别，不需要统计丰度
    """
    def __init__(self, work_id):
        super(GoAnnoModule, self).__init__(work_id)
        options = [
            {"name":"blastout","type":"infile","format":"align.blast.blast_xml"},#blast table表
            {"name": "go_result_dir", "type": "outfile", "format": "annotation.mg_anno_dir"}  # 设置结果文件后面要用
        ]
        self.add_option(options)
        self.go_align_tools = []
        self.blast2go_tool = self.add_tool("annotation.go.nr_go")
        self.go_anno_tool = self.add_tool("annotation.mg_go_annotation")
        self.go_anno_stat_tool = self.add_tool("bacgenome.go_annot")
        self.align_result_path = ''

    def check_options(self):
        if not self.option("blastout").is_set:
            raise OptionError("Must provide blastout file!")
        return True

    def run_blastgo(self):
        """
        运行blast2go得到注释结果
        :return:
        """
        blastout=self.option("blastout").prop['path']
        table = xml2table(blastout, self.work_dir + '/temp_blastable.xls')
        self.blast2go_tool.set_options({
            'blastout':table,
        })
        self.blast2go_tool.on('end',self.run_go_anno)
        self.blast2go_tool.run()

    def run_go_anno(self):
        """
        根据注释的id结果得到的一张level层级对应关系的大表生成各分类水平的表
        :return:
        """
        self.go_anno_tool.set_options({
            'blast2go_annot': self.blast2go_tool.option('blast2go_annot'),
        })
        self.go_anno_tool.on('end', self.run_go_stat)
        self.go_anno_tool.run()

    def run_go_stat(self):
        file = self.blast2go_tool.option('blast2go_annot')
        self.go_anno_stat_tool.set_options({
            'go2level_infile':file
        })
        self.go_anno_stat_tool.on('end',self.set_output)
        self.go_anno_stat_tool.run()

    def set_output(self):
        self.option('go_result_dir', self.output_dir)
        anno_allfiles = os.listdir(self.go_anno_tool.output_dir)
        out_oldfiles = [os.path.join(self.go_anno_tool.output_dir, i) for i in anno_allfiles]
        output_newfiles = [os.path.join(self.output_dir, i) for i in anno_allfiles]
        for i in range(len(output_newfiles)):
            if os.path.exists(output_newfiles[i]):
                os.remove(output_newfiles[i])
            os.link(out_oldfiles[i], output_newfiles[i])
        if os.path.exists(self.output_dir + "/go12level_statistics.xls"):
            os.remove(self.output_dir + "/go12level_statistics.xls")
        os.link(self.go_anno_tool.work_dir + "/go12level_statistics.xls", self.output_dir + "/go12level_statistics.xls")
        if os.path.exists(self.output_dir + "/go123level_statistics.xls"):
            os.remove(self.output_dir + "/go123level_statistics.xls")
        os.link(self.go_anno_tool.work_dir + "/go123level_statistics.xls", self.output_dir + "/go123level_statistics.xls")
        if os.path.exists(self.output_dir + "/query_gos.list"):
            os.remove(self.output_dir + "/query_gos.list")
        os.link(self.go_anno_tool.work_dir + "/GO.list", self.output_dir + "/query_gos.list")
        if os.path.exists(self.output_dir + "/blast2go.annot"):
            os.remove(self.output_dir + "/blast2go.annot")
        os.link(self.blast2go_tool.output_dir + "/go_annot.xls", self.output_dir + "/blast2go.annot")
        if os.path.exists(self.output_dir + "/go_statistics.xls"):
            os.remove(self.output_dir + "/go_statistics.xls")
        os.link(self.go_anno_stat_tool.output_dir + "/go_statistics.xls", self.output_dir + "/go_statistics.xls")
        self.end()

    def run(self):
        super(GoAnnoModule, self).run()
        self.run_blastgo()

    def end(self):
        super(GoAnnoModule, self).end()
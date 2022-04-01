# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify:2018.12.3

from biocluster.module import Module
import os,re
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import table_to_table


class AnnoGoStatModule(Module):
    def __init__(self, work_id):
        super(AnnoGoStatModule, self).__init__(work_id)
        options = [
            {"name": "go1234level_out", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_anno", "type": "infile", "format": "sequence.profile_table"},  # 基因注释file
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},  # 基因丰度表
            {"name": "go_result_dir", "type": "outfile", "format": "annotation.mg_anno_dir"},  # 设置结果文件后面要用
            {"name": "level", "type": "string","default":"all"},  # 各个分析计算单个水平丰度用
            {"name": "out_table", "type": "outfile", "format": "meta.otu.otu_table"} # 单水平输出丰度文件
        ]
        self.add_option(options)
        self.remove_gene = self.add_tool("annotation.go.remove_gene_go")
        self.go_anno_stat_tool = self.add_tool("annotation.mg_go_anno_stat")
        self.align_result_path = ''

    def check_options(self):
        if not self.option("go1234level_out").is_set:
            raise OptionError("Must provide blastout file!", code="21201901")
        if not self.option("reads_profile_table").is_set:
            raise OptionError("Must provide reads_profile_table file", code="21201902")
        return True

    def run_remove_gene(self):
        self.remove_gene.set_options({
            'go1234level_out':self.option("go1234level_out"),
            'gene_anno':self.option("gene_anno"),
        })
        self.remove_gene.on('end',self.run_go_stat)
        self.remove_gene.run()

    def run_go_stat(self):
        self.go_anno_stat_tool.set_options({
            'go_anno':self.remove_gene.option("gogene_out"),
            'reads_profile_table':self.option('reads_profile_table'),
            'level': self.option("level")
        })
        self.go_anno_stat_tool.on('end',self.set_output)
        self.go_anno_stat_tool.run()

    def set_output(self):
        if self.option("level") != "all":
            self.option("out_table",self.go_anno_stat_tool.option("level_out"))
        else:
            self.option('go_result_dir', self.output_dir)
            anno_allfiles = os.listdir(self.remove_gene.output_dir)
            stat_allfiles = os.listdir(self.go_anno_stat_tool.output_dir)
            out_oldfiles = [os.path.join(self.remove_gene.output_dir, i) for i in anno_allfiles]
            stat_oldfiles = [os.path.join(self.go_anno_stat_tool.output_dir, i) for i in stat_allfiles]
            out_oldfiles.extend(stat_oldfiles)
            out_name = anno_allfiles + stat_allfiles
            output_newfiles = [os.path.join(self.output_dir, i) for i in out_name]
            for i in range(len(output_newfiles)):
                if os.path.exists(output_newfiles[i]):
                    os.remove(output_newfiles[i])
                os.link(out_oldfiles[i], output_newfiles[i])
        self.end()

    def run(self):
        super(AnnoGoStatModule, self).run()
        self.run_remove_gene()

    def end(self):
        super(AnnoGoStatModule, self).end()

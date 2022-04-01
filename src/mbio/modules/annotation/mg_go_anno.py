# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify:2019.9.30

from biocluster.module import Module
import os,re
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import table_to_table


class MgGoAnnoModule(Module):
    def __init__(self, work_id):
        super(MgGoAnnoModule, self).__init__(work_id)
        options = [
            {"name":"blastout","type":"infile","format":"sequence.profile_table"},#blast table表
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},  # 基因丰度表
            {"name": "go_result_dir", "type": "outfile", "format": "annotation.mg_anno_dir"}  # 设置结果文件后面要用
        ]
        self.add_option(options)
        self.go_align_tools = []
        self.blast2go_tool = self.add_tool("annotation.go.nr_go")
        self.go_anno_tool = self.add_tool("annotation.mg_go_annotation")
        self.go_anno_stat_tool = self.add_tool("annotation.mg_go_anno_stat")
        self.align_result_path = ''

    def check_options(self):
        if not self.option("blastout").is_set:
            raise OptionError("Must provide blastout file!", code="21202001")
        if not self.option("reads_profile_table").is_set:
            raise OptionError("Must provide reads_profile_table file", code="21202002")
        return True

    def run_blastgo(self):
        with open(self.option('blastout').prop['path'], 'r') as f:
            files = f.readlines()
            des = files[2].split('\t')[10]
        if re.search(r'gi\|', des):
            blastout = self.blast2go_tool.work_dir + '/blast_nr_table.xls'
            table_to_table(self.option("blastout").prop['path'], blastout)
        else:
            blastout=self.option("blastout").prop['path']
        self.blast2go_tool.set_options({
            'blastout':blastout,
        })
        self.blast2go_tool.on('end',self.run_go_anno)
        self.blast2go_tool.run()

    def run_go_anno(self):
        self.go_anno_tool.set_options({
            'blast2go_annot': self.blast2go_tool.option('blast2go_annot'),
        })
        self.go_anno_tool.on('end', self.run_go_stat)
        self.go_anno_tool.run()

    def run_go_stat(self):
        file =self.go_anno_tool.output_dir + '/go1234level_statistics.xls'
        self.go_anno_stat_tool.set_options({
            'go_anno':file,
            'reads_profile_table':self.option('reads_profile_table')
        })
        self.go_anno_stat_tool.on('end',self.set_output)
        self.go_anno_stat_tool.run()

    def set_output(self):
        self.option('go_result_dir', self.output_dir)
        anno_allfiles = os.listdir(self.go_anno_tool.output_dir)
        stat_allfiles = os.listdir(self.go_anno_stat_tool.output_dir)
        out_oldfiles = [os.path.join(self.go_anno_tool.output_dir, i) for i in anno_allfiles]
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
        super(MgGoAnnoModule, self).run()
        self.run_blastgo()

    def end(self):
        super(MgGoAnnoModule, self).end()
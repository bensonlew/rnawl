#-*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'20180930

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class PfamAnnotationModule(Module):
    def __init__(self, work_id):
        super(PfamAnnotationModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},    #非冗余基因集的输入
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},  #基因丰度表
            {"name": "evalue", "type": "float", "default": 1e-5},   #默认的evalue值
            {"name": "lines", "type": "int", "default": 200000},
            {"name": "database", "type": "string"},  #pfam注释的数据库必须要加此参数--pfam
            {"name": "align_result", "type": "outfile", "format": "meta_genomic.hmmscan_table"},  # 设置结果文件后面要用
            {"name": "pfam_result_dir", "type": "outfile", "format": "annotation.mg_anno_dir"}, #设置pfam注释的结果文件目录
            {"name": "best", "type": "bool", "default": True}
        ]
        self.add_option(options)    #检查option是否是list格式，体重每个option是否是字典类型
        self.split_fasta = self.add_tool("sequence.split_fasta") #不使用统一的module，因为不试用于本分析
        self.pfam_align_tools = []
        self.cat_out = self.add_tool("align.cat_hmmscanout")
        self.anno = self.add_tool("annotation.pfam_anno")   #质控前后序列信息的统计
        self.pfam_anno_stat_tool = self.add_tool("annotation.pfam_anno_stat")
        self.align_result_path = ''

    def check_options(self):
        if not self.option('query').is_set:
            raise OptionError("必须设置参数query", code="21202801")
        if not self.option("reads_profile_table").is_set:
            raise OptionError("必须设置参数reads_profile_table", code="21202802")
        return True

    def run_split_fasta(self):
        self.split_fasta.set_options({  #######批量设置参数值，并传给tool的option
                                        "fasta": self.option("query"),
                                        "lines": self.option("lines"),
                                        })
        self.split_fasta.on('end', self.run_align)  #########split发生end事件触发run.align
        self.split_fasta.run()

    def run_align(self):
        self.align_result_path = os.path.join(self.work_dir, "pfam_align")
        if os.path.exists(self.align_result_path):
            pass
        else:
            os.mkdir(self.align_result_path)
        for f in os.listdir(self.split_fasta.output_dir):
            file_path = os.path.join(self.split_fasta.output_dir, f)
            align_hmmscan = self.add_tool('align.mg_hmmscan')
            align_hmmscan.set_options({
                "faa_file": file_path,
                "database": "pfam_v33.1", ## fix by qingchen.zhang
            })
            self.pfam_align_tools.append(align_hmmscan)
        if len(self.pfam_align_tools) > 1:
            self.on_rely(self.pfam_align_tools, self.run_cat_result)
        else:
            #self.align_result_path = self.pfam_align_tools[0].output_dir
            self.pfam_align_tools[0].on('end', self.run_cat_result)
        for tool in self.pfam_align_tools:
            tool.run()

    def run_cat_result(self):
        for i in self.pfam_align_tools:
            for f in os.listdir(i.output_dir):
                if os.path.splitext(f)[1] == '.domtblout':
                    file_path = os.path.join(i.output_dir, f)
                    new_path = os.path.join(self.align_result_path, os.path.basename(file_path))
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(file_path, new_path)
        self.cat_out.set_options({
            "hmmscan_out": self.align_result_path,
        })
        self.cat_out.on('end', self.run_anno)
        self.cat_out.run()

    def run_anno(self):
        self.anno.set_options({
            'hmmscan_result': os.path.join(self.cat_out.output_dir, 'hmmscan_result'),
        })
        self.anno.on('end', self.run_pfam_anno_stat)
        self.anno.run()

    def run_pfam_anno_stat(self):
        self.pfam_anno_stat_tool.set_options({
            "pfam_anno_table": self.anno.option('pfam_anno_result'),
            "reads_profile_table": self.option("reads_profile_table")
        })
        self.pfam_anno_stat_tool.on('end', self.set_output)
        self.pfam_anno_stat_tool.run()

    def set_output(self):
        self.option('pfam_result_dir', self.output_dir)
        anno_allfiles = os.listdir(self.anno.output_dir)
        stat_allfiles = os.listdir(self.pfam_anno_stat_tool.output_dir)
        out_oldfiles = [os.path.join(self.anno.output_dir, i) for i in anno_allfiles]
        stat_oldfiles = [os.path.join(self.pfam_anno_stat_tool.output_dir, i) for i in stat_allfiles]
        out_oldfiles.extend(stat_oldfiles)
        out_name = anno_allfiles + stat_allfiles
        output_newfiles = [os.path.join(self.output_dir, i) for i in out_name]
        for i in range(len(output_newfiles)):
            if os.path.exists(output_newfiles[i]):
                os.remove(output_newfiles[i])
            os.link(out_oldfiles[i], output_newfiles[i])
        self.end()

    def run(self):
        super(PfamAnnotationModule, self).run()
        self.run_split_fasta()

    def end(self):
        super(PfamAnnotationModule, self).end()


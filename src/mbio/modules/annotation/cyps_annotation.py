# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
#import unittest


class CypsAnnotationModule(Module):
    def __init__(self, work_id):
        super(CypsAnnotationModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},    #输入文件
            {"name": "lines", "type": "int", "default": 100000},    #将fasta序列文件进行拆分为多个文件
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},   #基因丰度表
            {"name": "evalue", "type":"float", "default": 1e-5}, #统计确定的evalue值
            {"name": "p450_result_dir", "type": "outfile", "format": "annotation.mg_anno_dir"}  #设置结果文件目录
        ]
        self.add_option(options)
        self.split_fasta = self.add_tool("sequence.split_fasta")
        self.cyps_align_tools = []
        self.cyps_anno_tool = self.add_tool("annotation.cyps_anno")
        self.cyps_anno_stat_tool = self.add_tool("annotation.cyps_anno_stat")
        self.align_result_path = ''

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="21202101")
        if not self.option("reads_profile_table").is_set:
            raise OptionError("必须设置参数reads_profile_table", code="21202102")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="21202103")
        return True

    def run_split_fasta(self):
        self.split_fasta.set_options({
            "fasta": self.option("query"),
            "lines": self.option("lines"),
        })
        self.split_fasta.on('end', self.cyps_align)
        self.split_fasta.run()

    def cyps_align(self):
        self.align_result_path = os.path.join(self.work_dir, "cyps_align")
        if os.path.exists(self.align_result_path):
            pass
        else:
            os.mkdir(self.align_result_path)
        for f in os.listdir(self.split_fasta.output_dir):
            file_path = os.path.join(self.split_fasta.output_dir, f)
            align_diamond = self.add_tool('align.meta_diamond')
            align_diamond.set_options({
                "query": file_path,
                "database": "p450",
                "evalue": self.option("evalue"),
                "sensitive": 0,
                "target_num": 1
            })
            self.cyps_align_tools.append(align_diamond)
        if len(self.cyps_align_tools) > 1:
            self.on_rely(self.cyps_align_tools, self.cyps_anno)
        else:
            self.cyps_align_tools[0].on('end', self.cyps_anno)
        for tool in self.cyps_align_tools:
            tool.run()

    def cyps_anno(self):
        xml_dir = self.align_result_path
        for i in self.cyps_align_tools:
            for f in os.listdir(i.output_dir):
                if os.path.splitext(f)[1] == '.xml_new':
                    file_path = os.path.join(i.output_dir, f)
                    new_path = os.path.join(self.align_result_path, os.path.basename(file_path))
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(file_path, new_path)
        self.cyps_anno_tool.set_options({
            "cyps_xml_dir": xml_dir
        })
        self.cyps_anno_tool.on('end', self.cyps_anno_stat)
        self.cyps_anno_tool.run()

    def cyps_anno_stat(self):
        self.cyps_anno_stat_tool.set_options({
            'cyps_anno_table': self.cyps_anno_tool.option('cyps_anno_result'),
            'reads_profile_table': self.option('reads_profile_table'),
        })
        self.cyps_anno_stat_tool.on('end', self.set_output)
        self.cyps_anno_stat_tool.run()

    def set_output(self):
        self.option("p450_result_dir", self.output_dir)
        anno_allfiles = os.listdir(self.cyps_anno_tool.output_dir)
        stat_allfiles = os.listdir(self.cyps_anno_stat_tool.output_dir)
        out_oldfiles = [os.path.join(self.cyps_anno_tool.output_dir, i) for i in anno_allfiles]
        stat_oldfiles = [os.path.join(self.cyps_anno_stat_tool.output_dir, i) for i in stat_allfiles]
        out_oldfiles.extend(stat_oldfiles)
        out_name = anno_allfiles + stat_allfiles
        output_newfiles = [os.path.join(self.output_dir, i ) for i in out_name]
        for i in range(len(output_newfiles)):
            if os.path.exists(output_newfiles[i]):
                os.remove(output_newfiles[i])
            os.link(out_oldfiles[i], output_newfiles[i])
        self.end()

    def run(self):
        super(CypsAnnotationModule, self).run()
        self.run_split_fasta()

    def end(self):
        super(CypsAnnotationModule, self).end()

"""
class TestFunction(unittest.TestCase):

    ##This is test for the tool. Just run this script to do test.

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Metagenomic" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "annotation.cyps_annotation",
            "instant": False,
            "options": dict(
                query="/mnt/ilustre/users/sanger-dev/sg-users/zhangqingchen/test/ArdbAnnotation/SplitFasta/output/fasta_2",
                reads_profile_table="/mnt/ilustre/users/sanger-dev/sg-users/zhangqingchen/test/pfam_test/tool/reads_number.xls",
                p450_result_dir="/mnt/ilustre/users/sanger-dev/sg-users/zhangqingchen/test/cyps_test/module/output",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
"""


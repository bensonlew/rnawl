# -*- coding: utf-8 -*-
# __author__ = 'caiping.shi'

import os
import unittest

from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class RepeatmaskerModule(Module):
    def __init__(self, work_id):
        super(RepeatmaskerModule, self).__init__(work_id)
        options = [
            {"name": "input_genome", "type": "infile", "format": "small_rna.fasta"},  # 参考基因组文件
            {"name": "gff", "type": "outfile", "format": "gene_structure.gff3"},  # 预测生成的GFF格式
            {"name": "repeat", "type": "outfile", "format": "small_rna.fasta"},  # 预测生成的重复fasta序列
            {"name": "repeat_result", "type": "string", "default": ""}
        ]
        self.add_option(options)
        self.splitfasta = self.add_tool("small_rna.srna.split_fasta_base")
        self.repeatmasker_merge = self.add_tool("small_rna.srna.repeatmasker_merge")
        self.repeatmasker_tools = []

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")

    def run_splitfasta(self):
        self.splitfasta.set_options({
            "fasta": self.option("input_genome"),
            "bases": 10000000,
        })
        self.splitfasta.on('end', self.run_repeatmasker)
        self.splitfasta.run()

    def run_repeatmasker(self):
        for f in os.listdir(self.splitfasta.output_dir):
            opts = {'input_genome': os.path.join(self.splitfasta.output_dir, f)}
            repeatmasker_tool = self.add_tool('small_rna.srna.repeatmasker')
            repeatmasker_tool.set_options(opts)
            self.repeatmasker_tools.append(repeatmasker_tool)
        if len(self.repeatmasker_tools) == 1:
            self.repeatmasker_tools[0].on("end", self.run_catrepeatmasker)
        else:
            self.on_rely(self.repeatmasker_tools, self.run_catrepeatmasker)
        for tool in self.repeatmasker_tools:
            tool.run()

    def run_catrepeatmasker(self):
        repeatmasker_merge = self.work_dir + "/repeatmasker_merge.out"
        i = 0
        with open(repeatmasker_merge, 'w') as w:
            for tool in self.repeatmasker_tools:
                i += 1
                for file in os.listdir(tool.work_dir):
                    if file.startswith("fasta") and file.endswith("out"):
                        with open(os.path.join(tool.work_dir, file), 'r') as r:
                            if i == 1:
                                lines = r.readlines()
                                if len(lines) > 3:
                                    for line in lines:
                                        w.write(line)
                                else:
                                    i = 1
                            else:
                                lines = r.readlines()
                                if len(lines) > 3:
                                    for line in lines[3:]:
                                        w.write(line)
        self.run_repeatmasker_merge()

    def run_repeatmasker_merge(self):
        self.repeatmasker_merge.set_options({
            "repeatmasker_merge": self.work_dir + "/repeatmasker_merge.out",
            "input_genome": self.option("input_genome"),
        })
        self.repeatmasker_merge.on('end', self.set_output)
        self.repeatmasker_merge.run()

    def set_output(self, event):
        if self.option("repeat_result") != None:
            for file in os.listdir(self.repeatmasker_merge.output_dir):
                file1 = os.path.join(self.output_dir, file)
                if os.path.exists(file1):
                    os.remove(file1)
                os.link(os.path.join(self.repeatmasker_merge.output_dir, file), file1)
                if file1.endswith(".gff"):
                    self.option('gff', file1)
                if file1.endswith(".fa"):
                    self.option('repeat', file1)
        else:
            for file in os.listdir(self.option("repeat_result")):
                file1 = os.path.join(self.option("repeat_result"), file)
                if os.path.exists(file1):
                    os.remove(file1)
                os.link(os.path.join(self.option("repeat_result"), file), file1)
                if file1.endswith(".gff"):
                    self.option('gff', file1)
                if file1.endswith(".fa"):
                    self.option('repeat', file1)
        self.end()

    def run(self):
        super(RepeatmaskerModule, self).run()
        if self.option("repeat_result") != None:
            self.run_splitfasta()
        else:
            self.set_output()

    def end(self):
        super(RepeatmaskerModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime

        data = {
            "id": "repeatmasker" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "module",
            "name": "small_rna.srna.repeatmasker",
            "instant": False,
            "options": dict(
                input_genome="/mnt/lustre/users/sanger/app/database/Genome_DB_finish/plants/Zea_mays"
                             "/ensemble_release_24/dna/Zea_mays.AGPv3.24.dna.toplevel.fa"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

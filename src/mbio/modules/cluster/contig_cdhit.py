# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify:2019.01.16

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class ContigCdhitModule(Module):
    def __init__(self, work_id):
        super(ContigCdhitModule, self).__init__(work_id)
        options = [
            {"name": "gene_tmp_fa", "type": "infile", "format": "sequence.fasta"},  # 样品单拼序列
            {"name": "number1", "type": "int", "default": 0},  # 样品单拼序列切分为几份，默认0表示按文件大小自动计算，指定某个整数时则按指定数量分割
            {"name": "uni_fasta", "type": "outfile", "format": "sequence.fasta"},  # 核酸序列
            {"name": "cdhit_identity", "type": "float", "default": 0.95},  # 给出cdhit的参数identity
            {"name": "cdhit_coverage", "type": "float", "default": 0.9},  # 给出cdhit的参数coverage
        ]
        self.add_option(options)
        self.cdhit = self.add_module("cluster.cdhit_contigs_sample")
        self.merge = self.add_tool("cluster.cdhit_merge_contig")


    def check_options(self):
        if not self.option("gene_tmp_fa").is_set:
            raise OptionError("必须提供样品单拼序列")
        if not 0.75 <= self.option("cdhit_identity") <= 1:
            raise OptionError("cdhit identity必须在0.75，1之间")
        if not 0 <= self.option("cdhit_coverage") <= 1:
            raise OptionError("cdhit coverage必须在0,1之间")
        if self.option("number1") < 0:
            raise OptionError("number1必须大于等于0")

    def cd_hit(self):
        self.cdhit.set_options({
            "gene_tmp_fa": self.option("gene_tmp_fa"),
            "number": self.number1,
            "out_dir": self.work_dir,
            "identity": self.option("cdhit_identity"),
            "coverage": self.option("cdhit_coverage"),
        })
        self.cdhit.on('end',self.merge_run)
        self.cdhit.run()

    def merge_run(self):
        opts = {
            "compare_dir": self.work_dir + '/gene.uniGeneset.fa.cd-hit-para-tmp',
            "clstr": 0,
            "num1":self.number1,
        }
        self.merge.set_options(opts)
        self.merge.on('end', self.set_output)
        self.merge.run()

    def set_output(self):
        if os.path.exists(self.output_dir + "/contigs.uniContigs.fa"):
            os.remove(self.output_dir + "/contigs.uniContigs.fa")
        os.link(self.merge.output_dir + "/contigs.uniContigs.fa",self.output_dir + "/contigs.uniContigs.fa")
        self.end()

    def run(self):
        super(ContigCdhitModule, self).run()
        if os.path.exists(self.work_dir + '/gene.uniGeneset.fa.cd-hit-para-tmp'):
            pass
        else:
            os.mkdir(self.work_dir + '/gene.uniGeneset.fa.cd-hit-para-tmp')
        if self.option("number1") == 0:
            self.number1 = os.path.getsize(self.option("gene_tmp_fa").prop['path']) / 300000000 + 1
        self.cd_hit()

    def end(self):
        super(ContigCdhitModule, self).end()
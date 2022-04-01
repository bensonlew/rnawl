#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == gaohao

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir


class MetagbinAssembleSeqModule(Module):
    """
    metabin的组装序列处理模块
    """
    def __init__(self, work_id):
        super(MetagbinAssembleSeqModule, self).__init__(work_id)
        options = [
            {"name": "assemble_dir", "type": "infile", "format": "metagbin.assemble_dir"}, ##组装文件夹
            {"name": "number1", "type": "int", "default": 0},  # 样品单拼序列切分为几份，默认0表示按文件大小自动计算，指定某个整数时则按指定数量分割
            {"name": "uni_fasta", "type": "outfile", "format": "sequence.fasta"},  # 核酸序列
            {"name": "cdhit_identity", "type": "float", "default": 0.95},  # 给出cdhit的参数identity
            {"name": "cdhit_coverage", "type": "float", "default": 0.9},  # 给出cdhit的参数coverage
        ]
        self.add_option(options)
        self.metagbin_seq = self.add_module('sequence.metagbin_seq')
        self.contig_cdhit = self.add_module("cluster.contig_cdhit")
        self.contig_split = self.add_tool("metagbin.split_large_contig")
        #self.mmseq = self.add_tool("metagbin.mmseq_seq")
        self.gz = self.add_tool("sequence.zip")

    def check_options(self):
        """
        检查参数
        """
        if not self.option("assemble_dir").is_set:
            raise OptionError("必须输入组装序列文件夹！")
        else:
            return True

    def run(self):
        super(MetagbinAssembleSeqModule, self).run()
        self.run_metagbin_seq()

    def run_metagbin_seq(self):
        self.metagbin_seq.set_options({
            "assemble_dir": self.option("assemble_dir"),
            })
        self.metagbin_seq.on('end', self.run_split)
        self.metagbin_seq.run()

    def run_split(self):
        self.contig_split.set_options({
          "sca_fa": self.metagbin_seq.option("out_fa"),
        })
        self.contig_split.on('end', self.run_contig_cdhit)
        self.contig_split.run()

    #def run_mmseqs(self):
    #    self.mmseq.set_options({
    #        "contig_fa": self.metagbin_seq.option("out_fa"),
    #    })
    #    self.mmseq.on("end",self.gzip_file)
    #    self.mmseq.run()

    def run_contig_cdhit(self):
        self.contig_cdhit.set_options({
            "gene_tmp_fa":self.contig_split.option('out'),
            "number1": self.option("number1"),
            "cdhit_identity":self.option("cdhit_identity"),
            "cdhit_coverage":self.option("cdhit_coverage"),
        })
        self.contig_cdhit.on('end', self.gzip_file)
        self.contig_cdhit.run()

    def run_cat(self):
       if not self.contig_split.option('out_large').is_set:
           if os.path.exists(self.work_dir + '/all.uniContigs.fa'):
               os.remove(self.work_dir + '/all.uniContigs.fa')
           os.link(self.contig_cdhit.output_dir + '/contigs.uniContigs.fa', self.work_dir + '/all.uniContigs.fa')
       else:
            os.system("cat %s %s >%s" %(self.contig_split.option('out_large').prop['path'],self.contig_cdhit.output_dir + '/contigs.uniContigs.fa', self.work_dir + '/all.uniContigs.fa'))

    def gzip_file(self):
        self.run_cat()
        self.gz.set_options({
            "file_path": self.work_dir + '/all.uniContigs.fa',
            "method": "gz",
            })
        self.gz.on('end', self.set_output)
        self.gz.run()

    def set_output(self):
        if os.path.exists(self.output_dir + '/all.uniContigs.fa'):
            os.remove(self.output_dir + '/all.uniContigs.fa')
        os.link(self.work_dir + '/all.uniContigs.fa' ,self.output_dir + '/all.uniContigs.fa')
        if os.path.exists(self.output_dir + '/all.uniContigs.fa.gz'):
            os.remove(self.output_dir + '/all.uniContigs.fa.gz')
        os.link(self.gz.output_dir + '/all.uniContigs.fa.gz' ,self.output_dir + '/all.uniContigs.fa.gz')
        self.end()

    def end(self):
        super(MetagbinAssembleSeqModule, self).end()
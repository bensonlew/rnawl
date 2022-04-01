#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == gaohao

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir


class AssembleSeqModule(Module):
    """
    解压改名字
    """
    def __init__(self, work_id):
        super(AssembleSeqModule, self).__init__(work_id)
        options = [
            {"name": "fasta", "type": "string"},
            {"name": "sample_name", "type": "string"},
            {"name": "out_fa", "type": "outfile", "format": "sequence.fasta"},
        ]
        self.add_option(options)
        self.gunzip_fasta = self.add_tool('sequence.fasta_ungz')
        self.rename = self.add_tool("metagbin.assemble_rename")

    def check_options(self):
        """
        检查参数
        """
        if not self.option("fasta"):
            raise OptionError("必须输入组装序列文件夹！")
        else:
            return True

    def run(self):
        super(AssembleSeqModule, self).run()
        self.run_ungiz()

    def run_ungiz(self):
        sample = self.option("sample_name")
        fasta = self.option("fasta")
        self.gunzip_fasta.set_options({
            "fasta": fasta,
            "sample_name": sample,
            })
        self.gunzip_fasta.on('end', self.run_rename)
        self.gunzip_fasta.run()

    def run_rename(self):
        sample = self.option("sample_name")
        self.rename.set_options({
            "assmle_fa": self.gunzip_fasta.option("out_fa"),
            "sample_name": sample,
        })
        self.rename.on('end', self.set_output)
        self.rename.run()

    def set_output(self):
        if os.path.exists(self.output_dir + '/' + self.option('sample_name') + '.scaf.fa'):
            os.remove(self.output_dir + '/' + self.option('sample_name') + '.scaf.fa')
        if os.path.exists(self.rename.output_dir + '/' + self.option('sample_name') + '.scaf.fa'):
            os.link(self.rename.output_dir + '/' + self.option('sample_name') + '.scaf.fa', self.output_dir + '/' + self.option('sample_name') +'.scaf.fa')
            self.option("out_fa",self.output_dir + '/' + self.option('sample_name') + '.scaf.fa')
        self.end()

    def end(self):
        super(AssembleSeqModule, self).end()
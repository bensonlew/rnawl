#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == gaohao

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir


class MetagbinSeqModule(Module):
    """
    对数据进行解压统计碱基质量
    """
    def __init__(self, work_id):
        super(MetagbinSeqModule, self).__init__(work_id)
        options = [
            {"name": "assemble_dir", "type": "infile", "format": "metagbin.assemble_dir"},
            {"name": "out_fa", "type": "outfile", "format": "sequence.fasta"},
        ]
        self.add_option(options)
        self.cat = self.add_tool("metagbin.cat_file")
        self.samplesi = ""
        self.tools = []

    def check_options(self):
        """
        检查参数
        """
        if not self.option("assemble_dir").is_set:
            raise OptionError("必须输入组装序列文件夹！")
        else:
            return True

    def run(self):
        super(MetagbinSeqModule, self).run()
        self.run_ungiz()

    def get_list(self):
        list_path = os.path.join(self.option("assemble_dir").prop['path'], "list.txt")
        if os.path.exists(list_path):
            self.logger.info(list_path)
        sample = {}
        with open(list_path, "rb") as l:
            for line in l.readlines()[1:]:
                line = line.strip().split()
                if len(line) == 2:
                    sample_path = os.path.join(self.option("assemble_dir").prop['path'], line[1])
                    sample[line[0]] = sample_path
                else:
                     raise OptionError('list.txt文件格式有误')
        return sample

    def run_ungiz(self):
        self.samples = self.get_list()
        samples = self.samples
        for sample in samples:
            assemble_seq = self.add_module('metagbin.assemble_seq')
            assemble_seq.set_options({
                "fasta": samples[sample],
                "sample_name": sample,
            })
            self.tools.append(assemble_seq)
        if len(self.tools) > 1:
            self.on_rely(self.tools, self.run_cat)
        else:
            self.tools[0].on('end', self.set_output)
        for tool in self.tools:
            tool.run()

    def run_cat(self):
        if not os.path.exists(self.cat.work_dir + '/seq_dir'):
            os.mkdir(self.cat.work_dir + '/seq_dir')
        for tool in self.tools:
            if os.path.exists(tool.option("out_fa").prop['path']):
                self.logger.info(tool.option("out_fa").prop['path'])
                file = os.path.basename(tool.option("out_fa").prop['path'])
                if os.path.exists(self.cat.work_dir + '/seq_dir/' + file):
                    os.remove(self.cat.work_dir + '/seq_dir/' + file)
                os.link(tool.option("out_fa").prop['path'], self.cat.work_dir + '/seq_dir/' + file)
        self.cat.set_options({
            "fa_dir":self.cat.work_dir + '/seq_dir',
        })
        self.cat.on("end",self.set_output)
        self.cat.run()

    def set_output(self):
        if len(self.tools) > 1:
            if os.path.exists(self.output_dir + '/all.scaf.fa'):
                os.remove(self.output_dir + '/all.scaf.fa')
            os.link(self.cat.output_dir + '/all.scaf.fa',self.output_dir + '/all.scaf.fa')
        else:
            if os.path.exists(self.output_dir + '/all.scaf.fa'):
                os.remove(self.output_dir + '/all.scaf.fa')
            os.link(self.tools[0].output_dir + '/' + self.samples.keys()[0] + '.scaf.fa' ,self.output_dir + '/all.scaf.fa')
        self.option("out_fa",self.output_dir + '/all.scaf.fa')
        self.end()

    def end(self):
        super(MetagbinSeqModule, self).end()
# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 2020619

import os
import re
import math
import time
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class FastqToFastaWorkflow(Workflow):
    """
    相关性小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FastqToFastaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fastq_input1", "type": "infile", "format": "sequence.fastq"},
            {"name": "fastq_input2", "type": "infile", "format": "sequence.fastq"},
            {"name": "reverse_complement", "type": "bool"},
            {"name": "trim_5", "type": "int", "default": 3},
            {"name": "trim_3", "type": "int", "default": 5},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.fastq_ungz = self.add_tool("tool_lab.fastq_ungz")
        self.fq2fa = self.add_tool("tool_lab.pair_fastq_to_fasta")
        self.sequence_treat = self.add_tool("tool_lab.sequence_treat")

    def check_options(self):
        if not self.option('fastq_input1'):
            raise OptionError("必须设置输入的fastq1文件")
        if not self.option('fastq_input2'):
            raise OptionError("必须设置输入的fastq2文件")

    def run_ungz(self):
        self.fastq_ungz.set_options({
            "fastq1": self.option("fastq_input1"),
            "fastq2": self.option("fastq_input2"),
        })
        self.fastq_ungz.on('end', self.set_output, 'fastq_ungz_dir')
        self.fastq_ungz.on('end', self.run_fq2fa)
        self.fastq_ungz.run()

    def run_fq2fa(self):
        self.logger.info("???????????????????????????????")
        self.logger.info(self.path_fast1)
        self.logger.info(self.path_fast2)
        self.fq2fa.set_options({
            # "fastq_input1": self.option("fastq_input1"),
            # "fastq_input2": self.option("fastq_input2"),
            "fastq_input1": self.path_fast1,
            "fastq_input2": self.path_fast2,
        })
        self.fq2fa.on('end', self.set_output, 'fq2fa_dir')
        self.logger.info((not self.option("reverse_complement")) & (self.option("trim_5") == 0) & (self.option("trim_3") == 0))
        self.logger.info((not self.option("reverse_complement")))
        self.logger.info(self.option("trim_5") == 0)
        self.logger.info(self.option("trim_3") == 0)
        if (not self.option("reverse_complement")) & (self.option("trim_5") == 0) & (self.option("trim_3") == 0):
            self.logger.info("不需要运行 seq_treat tool!")
            self.fq2fa.on('end', self.end)
        else:
            self.logger.info("需要运行 seq_treat tool!")
            self.fq2fa.on('end', self.run_seq_treat)
        self.fq2fa.run()

    def run_seq_treat(self):
        self.sequence_treat.set_options({
            "fasta": self.output_dir + "/fq2fa_dir/fasta",
            "reverse_complement": self.option("reverse_complement"),
            "trim_5": self.option("trim_5"),
            "trim_3": self.option("trim_3"),

        })
        self.sequence_treat.on('end', self.set_output, 'seq_treat_dir')
        self.sequence_treat.on('end', self.end)
        self.sequence_treat.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'fastq_ungz_dir':
            self.logger.info("------------ungz-------")
            self.linkdir(obj.output_dir, 'fastq_ungz_dir')
        if event['data'] == 'fq2fa_dir':
            self.logger.info("------------fq2fa-------")
            self.linkdir(obj.output_dir, 'fq2fa_dir')
        if event['data'] == 'seq_treat_dir':
            self.logger.info("------------seq_treat-------")
            self.linkdir(obj.output_dir, 'seq_treat_dir')

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    # def set_db(self):
    #     self.logger.info("保存结果到mongo")
    #     api_primer = self.api.api("tool_lab.primer_design")
    #     primer_result = os.path.join(self.output_dir, "primer_design/variation.result")
    #     api_primer.add_sg_primer_detail(self.option('main_id'), primer_result)
    #     self.end()

    def run(self):
        self.path_fast1 = self.option("fastq_input1").prop['path']
        self.path_fast2 = self.option("fastq_input2").prop['path']
        if re.search('\.gz$', self.option("fastq_input1").prop['path']):
            self.logger.info("需要解压缩")
            path1 = os.path.basename(self.path_fast1)
            path2 = os.path.basename(self.path_fast2)
            self.path_fast1 = self.work_dir + "/FastqUngz/unzip_dir/" + ".".join(path1.split('.')[0:-2]) + ".fq"
            self.path_fast2 = self.work_dir + "/FastqUngz/unzip_dir/" + ".".join(path2.split('.')[0:-2]) + ".fq"
            self.logger.info(self.path_fast1)
            self.logger.info(self.path_fast2)
            self.run_ungz()
        else:
            self.run_fq2fa()
        super(FastqToFastaWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(FastqToFastaWorkflow, self).end()

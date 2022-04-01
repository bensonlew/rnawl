# -*- coding: utf-8 -*-

# __author__ = 'linfang.jin'
# time: 2017/1/15 16:47

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.ref_rna.gene_structure.hdrs import HdrsFile
import subprocess
import re
import shutil
import json
import os
import glob
from mbio.files.sequence.fasta import FastaFile
from subprocess import Popen, PIPE

class AsprofileAsAgent(Agent):
    '''
    每个样本运行asprofile的extract-as,summarize_as,extract-as-fpkm三个程序
    '''

    def __init__(self, parent):
        super(AsprofileAsAgent, self).__init__(parent)
        options = [{"name": "sample_gtf", "type": "infile", "format": "ref_rna.assembly.gtf", "default": None},
                   {"name": "hdrs", "type": "infile", "format": "ref_rna.gene_structure.hdrs", "default": None},
                   {"name": "ref_gtf", "type": "infile", "format": "ref_rna.assembly.gtf", "default": None},
                   {"name": "sample_tmap", "type": "infile", "format": "ref_rna.gene_structure.tmap",
                    "default": None},
                   # {'name': "sample_as", "type": "outfile", "format": 'ref_rna.gene_structure.as'},
                   # {"name": "sample_as_nr", "type": "outfile", "format": "ref_rna.gene_structure.as"},
                   ]
        self.add_option(options)
        self.step.add_steps('asprofile_as')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def check_options(self):

        if self.option('sample_gtf') is None:
            raise OptionError("必须设置样本转录本文件")
        if self.option('hdrs') is None:
            raise OptionError("必须设置hdrs文件")
        if self.option('ref_gtf') is None:
            raise OptionError('必须设置参考基因组文件')
        if self.option('sample_tmap') is None:
            raise OptionError('必须设置样本转录本tmap文件')

    def step_start(self):
        self.step.asprofile_as.start()
        self.step.update()

    def step_end(self):
        self.step.asprofile_as.finish()
        self.step.update()

    def set_resource(self):
        '''
        所需资源
        :return:
        '''
        self._cpu = 10
        self._memory = '100G'

    def end(self):
        """
        agent结束后一些文件的操作

        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            [".*as.*", "gtf", "样本as信息文件"]
        ])

        super(AsprofileAsAgent, self).end()


class AsprofileAsTool(Tool):
    def __init__(self, config):
        super(AsprofileAsTool, self).__init__(config)
        self._version = "vb-1.0.4"
        self.script_path = self.config.SOFTWARE_DIR + "/bioinfo/rna/ASprofile.b-1.0.4/"
        self.python_path = 'program/Python/bin/python'
        self.perl_path = 'program/perl-5.24.0/bin/perl'

    def run_extract_as(self):
        """
        :return:
        """
        cmd = "{}extract-as  {}  {} -r {} {} ".format(self.script_path, self.option('sample_gtf').path,
                                                      self.option('hdrs').path, self.option('sample_tmap').path,
                                                      self.option('ref_gtf').path)
        self.logger.info('开始运行run_extract_as')
        self.logger.info(cmd)
        sample_name = os.path.basename(self.option('sample_gtf').path).split('_')[0]
        sample_as_path = os.path.join(self.output_dir, sample_name + ".as")
        try:
            pro = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            result = "".join(pro.stdout.readlines())
            with open(sample_as_path, "w") as w:
                w.write(result)
            pro.communicate()
            self.logger.info("run_extract_as运行完成")
            return True
        except subprocess.CalledProcessError:
            self.logger.info("run_extract_as运行失败")
            return False

    def run_sum_as(self):
        """
        :return:
        """
        sample_name = os.path.basename(self.option('sample_gtf').path).split('_')[0]
        prefix = os.path.join(self.output_dir, sample_name)
        sample_as_path = os.path.join(self.output_dir, sample_name + ".as")
        cmd = "{} {}summarize_as.pl  {}  {} -p {}".format(self.perl_path, self.script_path,
                                                          self.option('sample_gtf').path, sample_as_path,
                                                          prefix)
        self.logger.info('开始运行run_sum_as')
        command = self.add_command('run_sum_as', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_sum_as运行完成")
        else:
            self.set_error("run_sum_as运行出错!")


    def run(self):
        super(AsprofileAsTool, self).run()
        bo = self.run_extract_as()
        if bo:
            self.run_sum_as()
            self.end()
        else:
            self.end()

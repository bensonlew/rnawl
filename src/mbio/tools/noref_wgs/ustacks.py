# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last_modify: 20181220

import os
import re
import math
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class UstacksAgent(Agent):
    """
    fastq均一化
    """
    def __init__(self, parent=None):
        super(UstacksAgent, self).__init__(parent)
        options = [
            # {"name": "fastq", "type": "infile", "format": "sequence.fastq"},  # fastq文件
            {"name": "fastq", "type": "string"},
            {"name": "id", "type": "int"},  # a unique integer ID for this sample
            {"name": "max_distance", "type": "int", 'default': 4},  # Maximum distance allowed between stacks
            {"name": "min_depth", "type": "int", 'default': 2},  # Minimum depth of coverage
            {"name": "num_threads", "type": "int", 'default': 12}  # num_threads threads
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fastq"):
            raise OptionError("please input sample's fastq", code="35500105")
        if not self.option("id"):
            raise OptionError("please input a unique integer ID for this sample", code="35500106")

    def set_resource(self):
        self._cpu = self.option('num_threads') + 2
        # self._memory = "{}G".format(os.path.getsize(self.option("fastq").prop['path']) / 1024 / 1024 / 1024 + 5)
        self._memory = '50G'
        self.logger.info("设置的内存：{}".format(self._memory))

    def end(self):
        super(UstacksAgent, self).end()


class UstacksTool(Tool):
    def __init__(self, config):
        super(UstacksTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/bin")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64")
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.tag_stat = self.config.PACKAGE_DIR + "/noref_wgs/tag-stat.pl"
        self.ustacks = "/bioinfo/noRefWGS/stacks-2.2/bin/ustacks"

    def run_ustacks(self):
        """
        ustacks -f 1080.fastq.gz -o 03.ustacks/ -i 1 -M 6 -m 2 -p 8 --deleverage
        modified by hd 修改了参数的默认值-i 1 -M 4 -m 2 -p 12
        :return:
        """
        cmd = "{} -f {} -o {} -i {} -M {} -m {} -p {} --deleverage".format(self.ustacks,
                                                                           # self.option("fastq").prop['path'],
                                                                           self.option("fastq"),
                                                                           self.output_dir,
                                                                           self.option("id"),
                                                                           self.option("max_distance"),
                                                                           self.option("min_depth"),
                                                                           self.option('num_threads'))
        command = self.add_command("ustacks", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ustacks运行完成")
        else:
            self.set_error("ustacks运行失败", code="35500105")

    def run_tag_stat(self):
        """
        perl tag-stat.pl -i 03.ustacks/1080 -o 3.ustacks/1080.tags.stat -k 1080
         touch /03.ustacks/1080.check
        :return:
        """
        sample_name = os.path.basename(self.option("fastq")).split('.')[0]
        cmd = '{} {} -i {} -o {} -k {}'.format(self.perl_path, self.tag_stat,
                                               self.output_dir + '/{}'.format(sample_name),
                                               self.output_dir + "/{}.tags.stat".format(sample_name), sample_name)
        command = self.add_command("tag_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("tag_stat运行完成")
        else:
            self.set_error("tag_stat运行失败", code="35500106")
        os.system('touch {}/{}.check'.format(self.output_dir, sample_name))

    def run(self):
        super(UstacksTool, self).run()
        self.run_ustacks()
        self.run_tag_stat()
        self.end()

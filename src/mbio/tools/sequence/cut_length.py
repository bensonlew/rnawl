# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import subprocess


class CutLengthAgent(Agent):
    """
    cut_length:split_length.pl根据序列长度阈值对输入路径中的序列进行拆分，然后用combine_contig.pl将短序列合并为一个文件
    version 1.0
    author: guhaidong
    last_modify: 2017.09.13
    """

    def __init__(self, parent):
        super(CutLengthAgent, self).__init__(parent)
        options = [
            {"name": "contig", "type": "infile", "format": "sequence.fasta_dir"},  # 输入contig文件路径
            {"name": "cut_length", "type": "float", "default": "1000"},  # 拆分序列长度标准，默认1000
            {"name": "cut_contig", "type": "outfile", "format": "sequence.fasta_dir"},  # 输出contig文件路径
            {"name": "short_contig", "type": "outfile", "format": "sequence.fasta"},  # 输出供newbler拼接使用的contig文件
        ]
        self.add_option(options)
        self.step.add_steps('cut_length')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cut_length.start()
        self.step.update()

    def step_end(self):
        self.step.cut_length.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("contig").is_set:
            raise OptionError("请传入contig序列路径", code="34000401")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '2G'


class CutLengthTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(CutLengthTool, self).__init__(config)
        self.perl_path = '/miniconda2/bin/perl '
        # self.cut_length_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/split_length.pl '
        # self.combine_contig_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/combine_contig.pl '
        self.cut_length_path = self.config.PACKAGE_DIR + '/sequence/scripts/split_length.pl '
        self.combine_contig_path = self.config.PACKAGE_DIR + '/metagenomic/scripts/combine_contig.pl '
        self._version = 'v1.0'

    def cut_length(self):
        """
        合并序列
        :return:
        """
        files = os.listdir(self.option('contig').prop['path'])
        file_rout = self.option('contig').prop['path']
        n = 0
        fw1 = open(self.output_dir + '/more.list', 'w')
        fw2 = open(self.output_dir + '/less.list', 'w')
        for file in files:
            n += 1
            sample = file.split('.contig.fa')[0]
            fw1.write("{}/{}.more{}.fa\t{}\n".format(self.output_dir, sample, self.option('cut_length'), sample))
            fw2.write("{}/{}.less{}.fa\t{}\n".format(self.output_dir, sample, self.option('cut_length'), sample))
            cmd = self.perl_path + self.cut_length_path + " %s %s %s/%s" % (
                os.path.join(file_rout, file), self.option('cut_length'), self.output_dir, sample
            )
            self.logger.info('正在按%s长度拆分%s' % (self.option('cut_length'), file))
            command = self.add_command("cut_length_cmd_{}".format(n), cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("序列{}拆分成功".format(file))
            else:
                self.set_error("序列%s拆分失败！", variables=(file), code="34000401")
        fw1.close()
        fw2.close()

    def combine_contig(self):
        """
        将小于1000bp的reads合并在一起
        """
        cmd = self.perl_path + self.combine_contig_path + " -list %s -o %s -s Y" % (
            os.path.join(self.output_dir, "less.list"),
            os.path.join(self.output_dir, "newbler_input.fa"),
        )
        command = self.add_command("combine_contig", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("合并短fasta成功")
        else:
            self.set_error("合并短fasta失败", code="34000402")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        self.logger.info("设置结果目录")
        self.option('short_contig').set_path(os.path.join(self.output_dir, "newbler_input.fa"))
        self.option('cut_contig').set_path(self.output_dir)
        self.logger.info("设置结果目录成功")

    def run(self):
        super(CutLengthTool, self).run()
        self.cut_length()
        self.combine_contig()
        self.set_output()
        self.end()

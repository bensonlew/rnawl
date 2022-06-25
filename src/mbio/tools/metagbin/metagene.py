# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os
import re
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.ref_rna.trans_step import step_count

class MetageneAgent(Agent):
    """
    使用metagene软件进行基因预测
    version: Metagene
    author: wangzhaoyue
    last_modify: 2017.06.19
    """
    def __init__(self, parent):
        super(MetageneAgent, self).__init__(parent)
        options = [
            {"name": "cut_more_scaftig", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample_name", "type": "string"},  # 样本的名称
            {"name": "min_gene", "type": "string", "default": "100"},  # 输入最短基因长度，如100
            {"name": "fna", "type": "outfile", "format": "sequence.fasta"},  # 输出文件，样本的核酸序列
            {"name": "faa", "type": "outfile", "format": "sequence.fasta"},  # 输出文件，样本的氨基酸序列
            {"name": "gff", "type": "outfile", "format": "gene_structure.gff3,metagbin.file_table"},  # 输出文件，样本的gff
            {"name": "cut_more_fna", "type": "outfile", "format": "sequence.fasta"},  # 输出文件，样本去除最小值后的核酸序列
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('cut_more_scaftig'):
            raise OptionError('必须输入样本去掉小于最短contig长度的序列文件')
        if not self.option('sample_name'):
            raise OptionError('必须输入样本的名称')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        last modified: 20180320
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(MetageneAgent, self).end()

class MetageneTool(Tool):
    def __init__(self, config):
        super(MetageneTool, self).__init__(config)
        self._version = "metagene"
        # self.sh_path = 'bioinfo/metaGenomic/scripts/'
        # self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/gene_structure/scripts/'
        self.sh_path = self.config.PACKAGE_DIR + '/gene_structure/scripts/'
        self.metagene_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/metagene/'
        self.perl_path = '/miniconda2/bin/perl '
        self.metagene_seqs_path = self.config.PACKAGE_DIR + '/gene_structure/scripts/metagene_seqs.pl '
        self.cut_more_path = self.config.PACKAGE_DIR + '/sequence/scripts/cut_more.pl '
        self.emboss_path = "bioinfo/seq/EMBOSS-6.6.0/emboss/"
        self.LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/bioinfo/seq/EMBOSS-6.6.0/lib"
        self.set_environ(LD_LIBRARY_PATH=self.LD_LIBRARY_PATH)
        self.rerun_num = 0

    def run(self):
        """
        运行
        :return:
        """
        super(MetageneTool, self).run()
        self.run_metagene()
        self.end()

    def run_metagene(self):
        """
        metagene [input] -m > [output]
        :return:
        """
        cmd = self.sh_path + 'metagene.sh %s %s %s' % (self.metagene_path, self.option('cut_more_scaftig').prop['path'],
                                                       self.work_dir + '/' + self.option('sample_name') + '.metagene.csv')
        command = self.add_command("metagene", cmd, ignore_error=True, shell=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行metagene的cmd完成")
            self.run_metageneseqs()
        else:
            self.set_error("运行metagene的cmd出错!", code="32201102")

    def run_metageneseqs(self):
        """
        MetageneSeqs  -m [csv] -f [scaftig] -o [fna]
        :return:
        """
        cmd = self.perl_path + self.metagene_seqs_path + '-m %s -f %s -o %s' % (self.work_dir + '/' + self.option('sample_name') + '.metagene.csv', self.option('cut_more_scaftig').prop['path'], self.work_dir + '/' + self.option('sample_name') + '.metagene.fna')
        command = self.add_command("metageneseqs运行", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行MetageneSeqs的cmd完成")
            self.run_cut_more()
        else:
            self.set_error("运行MetageneSeqs的cmd运行出错!", code="32201103")

    def run_cut_more(self):
        """
        perl cut_more.pl [run_MetageneSeqs的输出文件] [最短contig长度] [输出文件的名称前缀]
        :return:
        """
        cmd = self.perl_path + self.cut_more_path + self.work_dir + '/' + self.option('sample_name') + '.metagene.fna' + ' ' + self.option('min_gene') + ' ' + self.option('sample_name') + '.metagene.fna'
        command = self.add_command("cut_more", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行cut_more完成")
            self.run_transeq()
        else:
            self.set_error("运行cut_more运行出错!", 32201104)

    def run_transeq(self):
        """
        transeq -sequence [fna.more] -table 11 -trim -outseq [faa]
        :return:
        """
        if os.path.exists(self.work_dir + '/' + self.option('sample_name') + '.metagene.faa'):
            os.remove(self.work_dir + '/' + self.option('sample_name') + '.metagene.faa')
        cmd = self.emboss_path + 'transeq -sequence %s -table 11 -trim -outseq %s' % (
            self.work_dir + '/' + self.option('sample_name') + '.metagene.fna.more' + self.option('min_gene'),
            self.work_dir + '/' + self.option('sample_name') + '.metagene.faa')
        command = self.add_command("transeq运行", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.set_output()
            self.logger.info("运行transeq的cmd完成")
        else:
            self.set_error("运行transeq的cmd运行出错!")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        outfasta1 = self.output_dir + '/' + self.option('sample_name') + '.metagene.fna'
        outfasta2 = self.output_dir + '/' + self.option('sample_name') + '.metagene.more' + self.option('min_gene') + '.fa'
        if os.path.exists(outfasta1):
            os.remove(outfasta1)
        if os.path.exists(outfasta2):
            os.remove(outfasta2)
        os.link(self.work_dir + '/' + self.option('sample_name') + '.metagene.fna.more' + self.option('min_gene'),
                outfasta2)
        os.link(self.work_dir + '/' + self.option('sample_name') + '.metagene.faa',
                self.output_dir + '/' + self.option('sample_name') + '.metagene.faa')
        self.option('cut_more_fna').set_path(outfasta2)
        self.option('faa').set_path(self.output_dir + '/' + self.option('sample_name') + '.metagene.faa')
        self.option("gff", self.work_dir + '/' + self.option('sample_name') + '.metagene.csv')
        self.logger.info("设置Metagene分析结果目录成功")
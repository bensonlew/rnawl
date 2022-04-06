# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import os
import re
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.ref_rna.trans_step import step_count

class MetageneStatAgent(Agent):
    """
    统计基因预测结果
    version: 1
    author: guhaidong
    last_modify: 2017.08.30
    """
    def __init__(self, parent):
        super(MetageneStatAgent, self).__init__(parent)
        options = [
            {"name": "contig_dir", "type": "infile", "format": "sequence.fasta_dir"},
            # 输入文件，预测后的序列路径
            {"name": "sample_stat", "type": "outfile", "format": "sequence.profile_table"},  # 输出文件，对各基因预测结果进行统计
            {"name": "fasta", "type": "outfile", "format": "sequence.fasta"},  # 输出核酸文件，输出序列用于构建非冗余基因集
            {"name": "faa", "type": "outfile", "format": "sequence.fasta"},  # 输出蛋白文件，输出序列用于构建非冗余基因集
            {"name": "fasta_sample", "type": "outfile", "format": "sequence.fasta"},  # 输出核酸文件，输出序列用于构建非冗余基因集
            {"name": "fasta_mix", "type": "outfile", "format": "sequence.fasta"},  # 输出核酸文件，输出序列用于构建非冗余基因集
            {"name": "faa_sample", "type": "outfile", "format": "sequence.fasta"},  # 输出单拼样本蛋白文件，输出序列用于构建非冗余基因集
            {"name": "faa_mix", "type": "outfile", "format": "sequence.fasta"},  # 输出混拼样本蛋白文件，输出序列用于构建非冗余基因集
        ]
        self.add_option(options)
        self.step.add_steps("MetageneStat")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.MetageneStat.start()
        self.step.update()

    def stepfinish(self):
        self.step.MetageneStat.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('contig_dir'):
            raise OptionError('必须输入基因预测结果文件路径', code="32201201")
        if not os.path.exists(self.option('contig_dir').prop['path']):
            raise OptionError('基因预测结果文件夹不存在', code="32201202")
        if not os.listdir(self.option('contig_dir').prop['path']):
            raise OptionError('基因预测结果文件夹为空: %s', variables=(self.option('contig_dir').prop['path']), code="32201203")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 5
        #self._memory = "50G"
        self._memory = str((50 + self._rerun_time * 20)) + 'G' # by xieshichang 20200424

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(MetageneStatAgent, self).end()


class MetageneStatTool(Tool):
    def __init__(self, config):
        super(MetageneStatTool, self).__init__(config)
        self._version = "1"
        self.python_path = '/miniconda2/bin/python '
        # self.gene_stat_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/gene_stat.py'
        self.gene_stat_path = self.config.PACKAGE_DIR + '/gene_structure/gene_stat.py'
        self.set_environ(LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + '/bioinfo/seq/EMBOSS-6.6.0/lib')
        self.emboss_path = "/bioinfo/seq/EMBOSS-6.6.0/emboss/"

    def run(self):
        """
        运行
        :return:
        """
        super(MetageneStatTool, self).run()
        self.run_metagenestat()
        self.set_output()
        self.end()

    def run_metagenestat(self):
        """
        gene_stat -gene_dir  gene.directory -output_stat  stat_file -output_fa  fasta_file
        :return:
        """
        cmd = self.python_path + ' %s -gene_dir %s -output_stat %s -output_fa %s' % (self.gene_stat_path,
         self.option('contig_dir').prop['path'], self.work_dir + '/sample.metagene.stat', self.work_dir + '/Total.metagene.fa' )
        command = self.add_command("metagenestat", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行metagenestat的cmd完成")
            self.run_transeq()
        else:
            self.set_error("运行metagenestat的cmd运行出错!", code="32201201")

    def run_transeq(self):
        """
        transeq -sequence [fna.more] -table 11 -trim -outseq [faa]
        :return:
        """
        if os.path.exists(self.work_dir + '/Total.metagene.fa'):
            cmd = self.emboss_path + 'transeq -sequence %s -table 11 -trim -outseq %s' % (
                self.work_dir + '/Total.metagene.fa', self.work_dir + '/Total.metagene.faa')
            command = self.add_command("transeq运行", cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行transeq的cmd完成")
            elif command.return_code in [-7]:
                self.logger.info("return code: %s" % command.return_code)
                self.add_state('memory_limit', 'memory is low!')
            else:
                self.set_error("运行transeq的cmd运行出错!")
        if os.path.exists(self.work_dir + '/Total.metagene.fa_sample'):
            cmd1 = self.emboss_path + 'transeq -sequence %s -table 11 -trim -outseq %s' % (
                self.work_dir + '/Total.metagene.fa_sample', self.work_dir + '/Total.metagene_sample.faa')
            command1 = self.add_command("sample_transeq运行", cmd1)
            command1.run()
            self.wait(command1)
            if command1.return_code == 0:
                self.logger.info("运行transeq的cmd1完成")
            else:
                self.set_error("运行transeq的cmd1运行出错!")
        if os.path.exists(self.work_dir + '/Total.metagene.fa_mix'):
            cmd2 = self.emboss_path + 'transeq -sequence %s -table 11 -trim -outseq %s' % (
                self.work_dir + '/Total.metagene.fa_mix', self.work_dir + '/Total.metagene_mix.faa')
            command2 = self.add_command("mix_transeq运行", cmd2)
            command2.run()
            self.wait(command2)
            if command2.return_code == 0:
                self.logger.info("运行transeq的cmd2完成")
            else:
                self.set_error("运行transeq的cmd2运行出错!")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir + '/sample.metagene.stat'):
            os.remove(self.output_dir + '/sample.metagene.stat')
        if os.path.exists(self.output_dir + '/Total.metagene.fa'):
            os.remove(self.output_dir + '/Total.metagene.fa')
        if os.path.exists(self.output_dir + '/Total.metagene.faa'):
            os.remove(self.output_dir + '/Total.metagene.faa')
        if os.path.exists(self.work_dir + '/Total.metagene.faa'):
            os.link(self.work_dir + '/Total.metagene.faa', self.output_dir + '/Total.metagene.faa')
            self.option('faa').set_path(self.output_dir + '/Total.metagene.faa')
        if os.path.exists(self.work_dir + '/sample.metagene.stat'):
            os.link(self.work_dir + '/sample.metagene.stat', self.output_dir + '/sample.metagene.stat')
            self.option('sample_stat').set_path(self.output_dir + '/sample.metagene.stat')
        if os.path.exists(self.work_dir + '/Total.metagene.fa'):
            os.link(self.work_dir + '/Total.metagene.fa', self.output_dir + '/Total.metagene.fa')
            self.option('fasta').set_path(self.output_dir + '/Total.metagene.fa')
        if os.path.exists(self.work_dir + '/Total.metagene_sample.faa'):
            self.option('faa_sample').set_path(self.work_dir + '/Total.metagene_sample.faa')
            self.option('fasta_sample').set_path(self.work_dir + '/Total.metagene.fa_sample')
        if os.path.exists(self.work_dir + '/Total.metagene_mix.faa'):
            self.option('faa_mix').set_path(self.work_dir + '/Total.metagene_mix.faa')
            self.option('fasta_mix').set_path(self.work_dir + '/Total.metagene.fa_mix')
        self.logger.info("设置Metagene分析结果目录成功")

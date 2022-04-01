# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.files.sequence.fasta import FastaFile
import os
import subprocess
from biocluster.core.exceptions import OptionError
import shutil


class CdhitMergeContigAgent(Agent):
    """
    version v1.0
    author: gaohao
    last modified:2019.01.16
    """

    def __init__(self, parent):
        super(CdhitMergeContigAgent, self).__init__(parent)
        options = [
            {"name": "compare_dir", "type": "infile", "format": "sequence.cdhit_cluster_dir"},  # 输入cd-hit比对后的文件夹
            {"name": "fa", "type": "outfile", "format": "sequence.fasta"},  # 非冗余组装序列
            {"name": "clstr", "type": "int", "default": 1},  # 是否生成cluster文件，0不生成，1生成。
            {"name": "num1", "type": "int", "default": 0},  # 单拼去冗余生成的.o文件个数。
        ]
        self.add_option(options)
        self.step.add_steps('cdhitmerge')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 50  # 每次重运行增加内存50G

    def step_start(self):
        self.step.cdhitmerge.start()
        self.step.update()

    def step_end(self):
        self.step.cdhitmerge.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("compare_dir").is_set:
            raise OptionError("必须设置参数compare_dir", code="31600301")
        if not self.option("clstr") in [0, 1]:
            raise OptionError("clstr参数必须为0或1", code="31600302")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = str(len(os.listdir(self.option("compare_dir").prop['path'])) / 6 + 1) + 'G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(CdhitMergeContigAgent, self).end()


class CdhitMergeContigTool(Tool):
    def __init__(self, config):
        super(CdhitMergeContigTool, self).__init__(config)
        self._version = '1.0'
        self.perl_path = '/program/perl-5.24.0/bin/perl '
        self.merge_path = self.config.SOFTWARE_DIR + '/bioinfo/uniGene/cd-hit-v4.6.1-2012-08-27/clstr_merge.pl'
        self.sort_clstr_path = self.config.SOFTWARE_DIR + '/bioinfo/uniGene/scripts/sort_clstr.pl'

    def run(self):
        super(CdhitMergeContigTool, self).run()
        self.fasta_merge()
        self.end()

    def fasta_merge(self):
        if os.path.isfile(self.work_dir + '/gene.uniGeneset.bak.clstr'):
            os.remove(self.work_dir + '/gene.uniGeneset.bak.clstr')
        else:
            pass
        if 'clstr' not in self.get_option_object().keys() or self.option("clstr") == 1:
            # 为保证工作流不出错而增加的判断 ^^^^^^^^^^^^
            for i in range(0, self.option("num1")):
                cmd2 = self.config.SOFTWARE_DIR + self.perl_path + self.merge_path + ' '
                clstr = self.option("compare_dir").prop['path'] + '/gene.geneset.tmp.fa.div-' + str(i) + '-/o.clstr '
                if i < self.option("num1") - 1:
                    for j in range(i + 1, self.option("num1")):
                        clstr += self.option("compare_dir").prop['path'] + '/gene.geneset.tmp.fa.div-' + str(
                            j) + '-/vs.' + str(i) + '.clstr '
                cmd2 += clstr + '>> ' + self.work_dir + '/gene.uniGeneset.bak.clstr'
                try:
                    subprocess.check_output(cmd2, shell=True)
                    self.logger.info("clstr" + str(i) + "succeed")
                except subprocess.CalledProcessError:
                    self.set_error("clstr %s failed", variables=(i), code="31600301")
        cmd1 = 'cat '
        if self.option("num1") > 0:
            for i in range(0, self.option("num1")):
                cmd1 += self.option("compare_dir").prop['path'] + '/gene.geneset.tmp.fa.div-' + str(i) + '-/o '

        cmd1 += ' > ' + self.work_dir + '/contigs.uniContigs.fa'
        self.logger.info(cmd1)
        try:
            subprocess.check_output(cmd1, shell=True)
            self.successful('fa')
        except subprocess.CalledProcessError:
            self.set_error("fasta failed", code="31600302")
            raise Exception("fasta failed")
        if 'clstr' not in self.get_option_object().keys() or self.option("clstr") == 1:
            # 为保证现有工作流不出错而加的判断
            cmd4 = '%s %s %s %s' % (self.perl_path, self.sort_clstr_path, self.work_dir + '/gene.uniGeneset.bak.clstr',
                                    self.work_dir + '/gene.uniGeneset.clstr')
            command4 = self.add_command('cmd_4', cmd4)
            command4.run()
            self.wait(command4)
            if command4.return_code == 0:
                self.logger.info("clstr succeed")
                os.remove(self.work_dir + '/gene.uniGeneset.bak.clstr')
            else:
                self.set_error("clstr failed", code="31600304")

    def successful(self, type):
        self.logger.info(type + " succeed")
        name = "contigs.uniContigs." + type
        filename = os.path.join(self.work_dir, name)
        linkfile = os.path.join(self.output_dir, name)
        if os.path.exists(linkfile):
            os.remove(linkfile)
        os.link(filename, linkfile)
        if type in ["fa"]:
            self.option(type, linkfile)
        else:
            pass


# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20190106

import os, sys
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class MetagbinSpadesAgent(Agent):
    """
    宏基因组binning用spades软件组装
    """
    def __init__(self, parent):
        super(MetagbinSpadesAgent, self).__init__(parent)
        options = [
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},  # 输入文件,cd-hit去冗余后的read1序列
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},  # 输入文件,cd-hit去冗余后的read2序列
            {"name": "fastqs", "type": "infile", "format": "sequence.fastq"},  # 输入文件,cd-hit去冗余后的reads序列
            {"name": "kmer", "type": "int", "default": "21,23,25,27,29,31,33,35,37,39,41"},  # k_mer值，例"39"
            {"name": "sample_name", "type": "string"},  #基因组名称
            {"name": "scaffold", "type": "outfile", "format": "sequence.fasta"},  # 输出文件
        ]
        self.add_option(options)
        self.step.add_steps("SPAdes")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self._memory_increase_step = 20

    def stepstart(self):
        self.step.SPAdes.start()
        self.step.update()

    def stepfinish(self):
        self.step.SPAdes.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fastq1'):
            raise OptionError('必须输入组装配置fastq1文件')
        if not self.option('fastq2'):
            raise OptionError('必须输入组装配置fastq2文件')
        if not self.option('sample_name'):
            raise OptionError('必须输入基因组名称')
        if not self.option('kmer'):
            raise OptionError('必须输入kmer值')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 8
        self._memory = "50G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(MetagbinSpadesAgent, self).end()

class MetagbinSpadesTool(Tool):
    def __init__(self, config):
        super(MetagbinSpadesTool, self).__init__(config)
        self.sample_name = self.option('sample_name')
        self.kmer = self.option('kmer')
        self.spades_path = '/bioinfo/metaGenomic/SPAdes-3.11.0-Linux/bin/spades.py'
        self.fastq1_path = self.option('fastq1').prop['path']
        self.fastq2_path = self.option('fastq2').prop['path']

    def run(self):
        """
        运行
        :return:
        """
        super(MetagbinSpadesTool, self).run()
        have_result = self.run_spdes()  # 如果已有拼接结果，则stat为1，如果没有，则stat为0
        self.set_output(have_result)
        self.end()


    def run_spdes(self):
        """
        用软件spades组装序列
        :return:
        """
        self.logger.info(self.spades_path)
        self.logger.info("正在进行SPAdes组装")
        if os.path.exists(self.work_dir + "/tmp_assemble"):
            shutil.rmtree(self.work_dir + "/tmp_assemble")
        os.mkdir(self.work_dir + "/tmp_assemble")
        out_path = (self.work_dir + "/tmp_assemble")
        if not self.option('fastqs').is_set:
            cmd = '{} --meta --only-assembler -1 {} -2 {} -t 8 -k {} -o {}'.format(self.spades_path, self.fastq1_path, self.fastq2_path, self.kmer, out_path)
            self.logger.info(cmd)
            if os.path.exists(self.work_dir + "/tmp_assemble" + "/" + "scaffolds.fasta"):
                self.logger.info("%s.%s.kmer.scafSeq已存在，跳过拼接" % (self.sample_name, self.kmer))
                result_stat = 1
            else:
                command = self.add_command("spades", cmd)
                command.run()
                self.wait(command)
                if command.return_code == 0:
                    self.logger.info("运行cmd完成")
                else:
                    self.set_error("运行cmd运行出错!")
                result_stat = 0
        else:
            self.fastqs_path = self.option('fastqs').prop['path']
            cmd1 = '{} --meta --only-assembler -1 {} -2 {} -s {} -t 8 -k {} -o {}'.format(self.spades_path, self.fastq1_path, self.fastq2_path, self.fastqs_path, self.kmer, out_path)
            self.logger.info(cmd1)
            if os.path.exists(self.work_dir + "/tmp_assemble" + "/" + "scaffolds.fasta"):
                self.logger.info("%s.%s.kmer.scafSeq已存在，跳过拼接" % (self.sample_name, self.kmer))
                result_stat = 1
            else:
                command = self.add_command("spades", cmd1)
                command.run()
                self.wait(command)
                if command.return_code == 0:
                    self.logger.info("运行cmd完成")
                else:
                    self.set_error("运行cmd运行出错!")
                result_stat = 0
        return result_stat

    def set_output(self, status):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if status == 0:
            self.logger.info("现在复制结果")
            if os.path.exists(self.output_dir + "/" + self.sample_name + "_"  + 'scaffolds.fasta'):
                os.remove(self.output_dir + "/" + self.sample_name +  'scaffolds.fasta')
            os.link(self.work_dir + "/" + "tmp_assemble" + "/" + "scaffolds.fasta",self.output_dir + "/"  + self.sample_name + "_"  + 'scaffolds.fasta')
        self.option('scaffold').set_path(self.output_dir + "/" + self.sample_name + "_"  + 'scaffolds.fasta')
        self.logger.info("设置SPAdes组装结果目录成功")
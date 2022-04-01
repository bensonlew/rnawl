# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 20119.05.09

import os
import re, shutil
from Bio import SeqIO
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir

class AddSampleAgent(Agent):
    """
    细菌比较基因组 泛基因组分析对上传文件或者基因预测的结果文件夹进行添加样本名称
    """
    def __init__(self, parent):
        super(AddSampleAgent, self).__init__(parent)
        options = [
            {"name": "is_merge", "type": "int", "default": 0},  # 是否需要合并
            {"name": "fasta_dir", "type": "infile", "format": "bac_comp_genome.input_dir"},  # 输入文件夹fna文件和faa文件
            {"name": "out", "type": "outfile", "format": "sequence.fasta"},#输出PAGP所需要格式数据
            {"name": "out_dir", "type": "outfile", "format": "bac_comp_genome.input_dir"},#输出PAGP所需要格式数据
        ]
        self.add_option(options)

    def check_options(self):
        # if not self.option("is_merge"):
        #     raise OptionError("必须设置参数is_merge!")
        if not self.option("fasta_dir").is_set:
            raise OptionError("必须设置参数fasta_dir文件夹!")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(AddSampleAgent, self).end()

class AddSampleTool(Tool):
    def __init__(self, config):
        super(AddSampleTool, self).__init__(config)
        self.merge = self.option("is_merge")
        self.fasta = self.option("fasta_dir").prop['path']
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'
        self.out = self.output_dir

    def change_name(self):
        """
        对文件的序列改名，加上样本名称
        :return:
        """
        self.logger.info("开始对序列进行改名")
        all_files = os.listdir(self.fasta)
        for file in all_files:
            file_path = os.path.join(self.fasta, file)
            if file.endswith('_CDS.faa'):
                file_name = file.split('_CDS.faa')[0]
                outfile_path = os.path.join(self.out, file_name + '.faa')
                with open(outfile_path, 'w') as outf:
                    for seq_record in SeqIO.parse(file_path, 'fasta'):
                        seq_id = seq_record.id
                        seq = seq_record.seq
                        outf.write('>{}|{}\n{}\n'.format(file_name,seq_id, seq))
            elif file.endswith('_CDS.fna'):
                file_name = file.split('_CDS.fna')[0]
                outfile_path = os.path.join(self.out, file_name + '.fna')
                with open(outfile_path, 'w') as outm:
                    for seq_record in SeqIO.parse(file_path, 'fasta'):
                        seq_id = seq_record.id
                        seq = seq_record.seq
                        outm.write('>{}|{}\n{}\n'.format(file_name,seq_id, seq))
        self.logger.info("序列改名完成")

    def run_merge(self):
        """
        开始运行整理
        :return:
        """
        self.logger.info("开始对样本进行合并")
        file_list = os.listdir(self.out)
        cmd = self.sh_path + 'cat_seq.sh'
        for file in file_list:
            file_path = self.out + '/' + file
            if file.endswith('.faa'):
                cmd += ' ' + file_path
        cmd += ' ' + self.work_dir + '/all_fasta.faa'
        self.logger.info('运行cat_seq，将table进行合并')
        self.logger.info(cmd)
        command = self.add_command("merge", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("合并文件完成啦")
        else:
            self.set_error("合并文件失败！")

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        self.logger.info("开始设置tool结果文件目录")
        if os.path.exists(self.output_dir):
            self.option("out_dir", self.output_dir)
        if os.path.exists(self.work_dir + '/all_fasta.faa'):
            self.option("out", self.work_dir + '/all_fasta.faa')


    def run(self):
        """
        运行
        :return:
        """
        self.logger.info("开始运行")
        super(AddSampleTool, self).run()
        self.change_name()
        if self.merge == 1:
            self.run_merge()
        self.set_output()
        self.end()

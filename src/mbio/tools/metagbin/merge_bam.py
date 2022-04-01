# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'@20190114
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from Bio import SeqIO


class MergeBamAgent(Agent):
    """
    宏基因组binning将sam文件转为fasta格式进行去冗余，将比对结果进行合并，并转为fastq和fasta文件
    """
    def __init__(self, parent):
        super(MergeBamAgent, self).__init__(parent)
        options = [
            {"name": "sort_file", "type": "infile", "format": "align.bwa.bam_dir"}, #输入sort后的文件夹
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # 输入文件,参考基因组bin的fasta序列
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},    #用于兼容s端序列
            #{"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 用于兼容s端序列
            {"name": "analysis", "type": "string", "default": "metagbin"},
            #{"name": "fasta1", "type": "outfile", "format": "sequence.fasta"}, #输出文件read1的fasta
            #{"name": "fasta2", "type": "outfile", "format": "sequence.fasta"}, #输出文件read2的fasta
            #{"name": "fasta1", "type": "outfile", "format": "sequence.fasta"}, #输出文件read1的fasta
            #{"name": "fasta2", "type": "outfile", "format": "sequence.fasta"}, #输出文件read2的fasta
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sort_file").is_set:
            raise OptionError("必须添加sort_file文件！")

    def set_resource(self):
        sort_bam = self.option('sort_file').prop['path']
        all_files = os.listdir(sort_bam)
        if self.option("analysis") in ['metagbin']:
            self._cpu = 8
            self._memory = '50G'
        else:
            self._cpu = 8
            num = int(len(all_files)/3)
            self._memory = str(num) + "G"

    def end(self):
        super(MergeBamAgent, self).end()

class MergeBamTool(Tool):
    def __init__(self, config):
        super(MergeBamTool, self).__init__(config)
        self.samtools = "bioinfo/align/samtools-1.4/bin/samtools"
        #self.depth = self.work_dir + "/" + 'cov.depth.txt'

    def run_merge(self):
        """
        将sort结果进行merge
        :return:
        """
        if os.path.exists(self.work_dir + '/merge.bam'):
            os.remove(self.work_dir + '/merge.bam')
        self.logger.info(self.option('sort_file'))
        self.logger.info(self.option('sort_file').prop.keys())
        sort_bam = self.option('sort_file').prop['path']
        all_files = os.listdir(sort_bam)
        file_list = []
        for file in all_files:
            if os.path.splitext(file)[1] == '.bam':
                file_path = os.path.join(sort_bam, file)
                file_list.append(file_path)
                self.logger.info('%s'%(file_list))
        merge_file = " ".join(str(i) for i in file_list)
        self.logger.info('开始对sort文件进行merge')
        self.merge_bam = self.work_dir + '/merge.bam'
        self.logger.info(self.merge_bam)
        cmd_bam = "{} merge --threads 8 {} {}".format(self.samtools, self.merge_bam, merge_file)
        self.logger.info(cmd_bam)
        to_bam = 'to_bam'
        command1 = self.add_command(to_bam, cmd_bam).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("%s运行完成"%to_bam)
        elif command1.return_code == 1:
            self.add_state('merge.bam文件已经存在', 'memory is low!')
        else:
            self.set_error("%s运行失败", code="")

    def run_sort_bam(self):
        """
        将bam文件进行排序
        :return:
        """
        self.merge_bam = self.work_dir + '/merge.bam'
        cmd_sort = "{} sort {} -o {} --threads 8 ".format(self.samtools,self.merge_bam,self.work_dir + '/merge.sort.bam')
        self.logger.info(cmd_sort)
        command = self.add_command('run_sort_bam', cmd_sort).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_sort_bam运行完成")
        else:
            self.set_error("run_sort_bam运行失败")

    def run_index_bam(self):
        """
        将bam文件构建索引
        :return:
        """
        cmd_index = "{} index {}".format(self.samtools,self.work_dir + '/merge.sort.bam')
        self.logger.info(cmd_index)
        command = self.add_command('run_index_bam', cmd_index).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_index_bam运行完成")
        else:
            self.set_error("run_index_bam运行失败")

    def run_convert_fasta(self):
        """
        将bam文件转为fasta格式
        :return:
        """
        merge_bam = self.work_dir + '/merge.sort.bam'
        mk = 2
        for file in os.listdir(self.option("fastq_dir").prop['path']):
            if re.search(r's.fq', file):
                mk =1
                break
        if mk ==1:
            cmd_fasta = "{} fasta -1 {} -2 {} -s {} {}".format(self.samtools, self.output_dir + "/read.1.fa", self.output_dir + "/read.2.fa", self.output_dir + "/read.s.fa", merge_bam)
        else:
            cmd_fasta = "{} fasta -1 {} -2 {} {}".format(self.samtools, self.output_dir + "/read.1.fa", self.output_dir + "/read.2.fa", merge_bam)
        self.logger.info(cmd_fasta)
        to_fasta = 'to_fasta'
        command = self.add_command(to_fasta, cmd_fasta).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("%s运行完成"%to_fasta)
        else:
            self.set_error("%s运行失败")

    def run_convert_fastq(self):
        """
        将bam文件转为fastq格式
        :return:
        """
        merge_bam = self.work_dir + '/merge.sort.bam'
        mk = 2
        for file in os.listdir(self.option("fastq_dir").prop['path']):
            if re.search(r's.fq', file):
                mk =1
                break
        if mk ==1:
            cmd_fastq = "{} fastq -1 {} -2 {} -s {} {}".format(self.samtools, self.output_dir + "/read.1.fastq", self.output_dir + "/read.2.fastq", self.output_dir + "/read.s.fastq", merge_bam)
        else:
            cmd_fastq = "{} fastq -1 {} -2 {} {}".format(self.samtools, self.output_dir + "/read.1.fastq", self.output_dir + "/read.2.fastq", merge_bam)
        self.logger.info(cmd_fastq)
        to_fastq = 'to_fastq'
        command = self.add_command(to_fastq, cmd_fastq).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("%s运行完成"%to_fastq)
        else:
            self.set_error("运行失败")

    def run_depth(self):
        """
        计算每个位点的测序深度
        :return:
        """
        self.logger.info('开始计算覆盖度')
        self.tmp_sh = self.work_dir + '/tmp.sh'
        self.depth = self.work_dir + '/depth.txt'
        if len(os.listdir(self.option('sort_file').prop['path'])) >= 2:
            merge_bam = self.work_dir + '/merge.sort.bam'
        else:
            merge_bam = self.work_dir + '/1_sort.bam'
        self.tmp_sh = "{} {} {} {}".format(self.config.PACKAGE_DIR + '/bacgenome/samtools_depth.sh', self.config.SOFTWARE_DIR ,merge_bam, self.depth)
        cmd1 = '/program/sh {}'.format(self.tmp_sh)
        self.logger.info('1111111111')
        command1 = self.add_command('depth', cmd1).run()
        self.logger.info('bbbbbbbbb')
        self.wait(command1)
        self.logger.info('ccccccccc')
        if command1.return_code == 0:
            self.logger.info("%s运行成功"%cmd1)
        else:
            self.set_error("%s运行失败", variables=(cmd1))

    def cal_coverage(self):
        """
        计算整个基因组的测序深度
        :return:
        """
        self.depth = self.work_dir + '/depth.txt'
        total_base = 0
        value_r = 0
        coverage_path = self.output_dir + '/Bin_coverage.xls'
        with open(self.depth, 'r') as f, open(coverage_path, 'wb') as w:
            w.write("Coverage\n")
            for line in f:
                line = line.strip().split('\t')
                coverage = int(line[2])
                total_base = total_base + coverage
            for seq_record in SeqIO.parse(self.option('ref_fa').prop['path'], 'fasta'):
                length = int(len(seq_record.seq))
                value_r += length
            final_coverage = int(total_base/value_r)
            w.write('{}\n'.format(final_coverage))

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        if self.option("analysis") in ['metagbin']:
            if os.path.exists(self.output_dir + "/" + "read.1.fastq"):
                self.logger.info('fastq结果文件生成成功')
            else:
                self.set_error('fastq结果文件未生成')
        else:
            for i in ['merge.bam']:
                if os.path.exists(self.output_dir + '/' + i):
                    os.remove(self.output_dir + '/' + i)
                os.link(self.work_dir + '/' + i, self.output_dir + '/' + i)

    def run(self):
        super(MergeBamTool, self).run()
        if self.option("analysis") in ['metagbin']:
            self.run_merge()
            self.run_sort_bam()
            self.run_index_bam()
            self.run_convert_fasta()
            self.run_convert_fastq()
            self.run_depth()
            self.cal_coverage()
            self.set_output()
            self.end()
        else:
            self.run_merge()
            self.set_output()
            self.end()




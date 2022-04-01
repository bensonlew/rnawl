# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from Bio import SeqIO
import re
import os


class ChooseFastqAgent(Agent):
    """
    宏基因组用于从去冗余的fa文件中挑选id，再从fastq文件中挑选出序列，最后统计出来coverage
    """
    def __init__(self, parent):
        super(ChooseFastqAgent, self).__init__(parent)
        options = [
            {"name": "fasta1", "type": "infile", "format": "sequence.fasta"},#去冗余后的fasta1
            {"name": "fasta2", "type": "infile", "format": "sequence.fasta"},#去冗余后的fasta2
            {"name": "fastas", "type": "infile", "format": "sequence.fasta"},#去冗余后的fastas
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},#提取比对上的fastq1
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},#提取比对上的fastq2
            {"name": "fastqs", "type": "infile", "format": "sequence.fastq"},#提取比对上的fastqs
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # 输入文件,参考基因组bin的fasta序列
            {"name": "coverage", "type": "outfile", "format": "sequence.profile_table"},#计算的coverage结果
        ]
        self.add_option(options)
        self.step.add_steps('choose_fastq')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 50  # 每次重运行增加内存50G by hao.gao @ 2019.10.31

    def step_start(self):
        self.step.choose_fastq.start()
        self.step.update()

    def step_end(self):
        self.step.choose_fastq.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查函数
        :return:
        """
        if not self.option("fasta1").is_set:
            raise OptionError("必须设置参数fasta1")
        if not self.option("fasta2").is_set:
            raise OptionError("必须设置参数fasta2")
        if not self.option("fastq1").is_set:
            raise OptionError("必须设置参数fastq1")
        if not self.option("fastq2").is_set:
            raise OptionError("必须设置参数fastq2")

    def set_resource(self):
        """
        设置所需资源
        :return:
        """
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(ChooseFastqAgent, self).end()


class ChooseFastqTool(Tool):
    def __init__(self,config):
        super(ChooseFastqTool, self).__init__(config)
        #self._version = '1.0'
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        self.script = self.config.PACKAGE_DIR + '/metagbin/'

    def extract_fasta_name(self):
        """
        从fasta序列中提取出index，

        将相同部分的提取出来去fastq中抽取read1和read2，
        剩下不相同的部分合并到singleton，用于组装
        :return:
        """
        self.logger.info('正在提取序列id')
        read_1_list = self.work_dir + '/read.1.list'
        read_2_list = self.work_dir + '/read.2.list'
        read_s_list = self.work_dir + '/read.s.list'
        with open(read_1_list, 'w+') as w:
            for seq_record in SeqIO.parse(self.option('fasta1').prop['path'], 'fasta'):
                id = seq_record.id
                w.write('{}\n'.format(id))
        with open(read_2_list, 'w+') as w1:
            for seq_record in SeqIO.parse(self.option('fasta2').prop['path'], 'fasta'):
                id = seq_record.id
                w1.write('{}\n'.format(id))
        if self.option('fastas').is_set:
            with open(read_s_list, 'w+') as w2:
                for seq_record in SeqIO.parse(self.option('fastas').prop['path'], 'fasta'):
                    id = seq_record.id
                    w2.write('{}\n'.format(id))
        else:
            os.system("touch %s" %(read_s_list))
        self.logger.info('正在生成list的id')

    def compare_fasta_name(self):
        """
        比较read1和read2的id
        :return:
        """
        self.logger.info('正在比较序列的id')
        fasta1 = self.work_dir + '/read.1.list'
        fasta2 = self.work_dir + '/read.2.list'
        result_1_list = self.work_dir + '/result.1.list'
        result_2_list = self.work_dir + '/result.2.list'
        result_same_list = self.work_dir + '/result.same.list'
        cmd1 = '{} {}compare_reads_name.pl {} {} {} {}'.format(self.perl_path, self.script, fasta1, fasta2, result_same_list, result_2_list)
        to_list2 = 'to_list2'
        self.logger.info(cmd1)
        command1 = self.add_command(to_list2, cmd1).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("%s运行完成"%to_list2)
        else:
            self.set_error("%s运行失败")
        cmd2 = '{} {}compare_reads_name.pl {} {} {} {}'.format(self.perl_path, self.script, fasta2,fasta1, result_same_list, result_1_list)
        to_list1 = 'to_list1'
        self.logger.info(cmd2)
        command2 = self.add_command(to_list1, cmd2).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("%s运行完成"%to_list1)
        else:
            self.set_error("%s运行失败")
        self.logger.info('正在生成最后的list文件')

    def choose_seq(self):
        """
        根据list文件从去冗余的fastq文件中抽取序列,
        并将single1和single2合并到singleton中
        :return:
        """
        result_same_list = self.work_dir + '/result.same.list'
        result_1_list = self.work_dir + '/result.1.list'
        result_2_list = self.work_dir + '/result.2.list'
        fastq1 = self.option('fastq1').prop['path']
        fastq2 = self.option('fastq2').prop['path']
        cmd1 = '{} {}choose_fastq.pl {} {} {}'.format(self.perl_path, self.script, fastq1, result_same_list, self.output_dir + '/last.read.1.fastq')
        to_fastq1 = 'to_fastq1'
        self.logger.info(cmd1)
        command1 = self.add_command(to_fastq1, cmd1, ignore_error=True).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("%s运行完成"%to_fastq1)
        elif command1.return_code in [1, -9, -7, 250,137,245, 255, 247, 249, -1]:
            self.logger.info("return code: %s" % command1.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by gaohao @20191227
        else:
            self.set_error("%s运行失败")

        cmd2 = '{} {}choose_fastq.pl {} {} {}'.format(self.perl_path, self.script, fastq2, result_same_list, self.output_dir + '/last.read.2.fastq')
        to_fastq2 = 'to_fastq2'
        self.logger.info(cmd2)
        command2 = self.add_command(to_fastq2, cmd2, ignore_error=True).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("%s运行完成"%to_fastq2)
        elif command2.return_code in [1, -9, -7, 250,137,245, 255, 247, 249, -1]:
            self.logger.info("return code: %s" % command2.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by gaohao @20191227
        else:
            self.set_error("%s运行失败")

        if self.option('fastas').is_set:
            result_s_list = self.work_dir + '/read.s.list'
            fastqs = self.option('fastqs').prop['path']
            cmds = '{} {}choose_fastq.pl {} {} {}'.format(self.perl_path, self.script, fastqs, result_s_list, self.output_dir + '/last.read.s.fastq')
            to_fastqs = 'to_fastqs'
            self.logger.info(cmds)
            commands = self.add_command(to_fastqs, cmds, ignore_error=True).run()
            self.wait(commands)
            if commands.return_code == 0:
                self.logger.info("%s运行完成"%to_fastqs)
            elif commands.return_code in [1, -9, -7, 250,137,245, 255, 247, 249, -1]:
                self.logger.info("return code: %s" % commands.return_code)
                self.add_state('memory_limit', 'memory is low!')   # add memory limit error by gaohao @20191227
            else:
                self.set_error("%s运行失败")
        else:
            os.system('touch %s'%(self.output_dir + '/last.read.s.fastq'))

        cmd3 = '{} {}choose_fastq.pl {} {} {}'.format(self.perl_path, self.script, fastq1, result_1_list, self.output_dir + '/left.read.1.fastq')
        to_fastq3 = 'to_fastq3'
        self.logger.info(cmd3)
        command3 = self.add_command(to_fastq3, cmd3, ignore_error=True).run()
        self.wait(command3)
        if command3.return_code == 0:
            self.logger.info("%s运行完成"%to_fastq3)
        elif command3.return_code in [1, -9, -7, 250,137,245, 255, 247, 249, -1]:
            self.logger.info("return code: %s" % command3.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by gaohao @20191227
        else:
            self.set_error("%s运行失败")

        cmd4 = '{} {}choose_fastq.pl {} {} {}'.format(self.perl_path, self.script, fastq1, result_2_list, self.output_dir + '/left.read.2.fastq')
        to_fastq4 = 'to_fastq4'
        self.logger.info(cmd4)
        command4 = self.add_command(to_fastq4, cmd4, ignore_error=True).run()
        self.wait(command4)
        if command4.return_code == 0:
            self.logger.info("%s运行完成"%to_fastq4)
        elif command4.return_code in [1, -9, -7, 250,137,245, 255, 247, 249, -1]:
            self.logger.info("return code: %s" % command4.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by gaohao @20191227
        else:
            self.set_error("%s运行失败")
        cmd5 = 'cat %s >> %s'%(self.output_dir + '/left.read.1.fastq', self.output_dir + '/last.read.s.fastq')
        os.system(cmd5)
        self.logger.info('将未匹配上的read1合并到reads上完成,命令：%s' % cmd5)
        cmd6 = 'cat %s >> %s'%(self.output_dir + '/left.read.2.fastq', self.output_dir + '/last.read.s.fastq')
        os.system(cmd6)
        self.logger.info('将未匹配上的read2合并到reads上完成，命令：%s' % cmd6)

    def cal_coverage(self): #没用
        """
        统计 reads的测序深度
        :return:
        """
        value_1 = 0
        value_2 = 0
        value_s = 0
        value_r = 0
        coverage_path = self.output_dir + '/coverage.txt'
        with open(coverage_path, 'wb') as w:
            #w.write("#Coverage\n")
            for seq_record in SeqIO.parse(self.output_dir + '/last.read.1.fastq', 'fastq'):
                length = int(len(seq_record.seq))
                value_1 += length
            self.logger.info("value_1: %s" %value_1)
            for seq_record in SeqIO.parse(self.output_dir + '/last.read.2.fastq', 'fastq'):
                length = int(len(seq_record.seq))
                value_2 += length
            self.logger.info("value_2: %s" %value_2)
            for seq_record in SeqIO.parse(self.output_dir + '/last.read.s.fastq', 'fastq'):
                length = int(len(seq_record.seq))
                value_s += length
            self.logger.info("value_s: %s" %value_s)
            for seq_record in SeqIO.parse(self.option('ref_fa').prop['path'], 'fasta'):
                length = int(len(seq_record.seq))
                value_r += length
            self.logger.info("value_r: %s" %value_r)
            coverage = float((value_1 + value_2 + value_s)/value_r)
            self.logger.info("coverage: %s" %coverage)
            w.write('#Coverage\n{}\n'.format(coverage))

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        self.logger.info("正在生成结果目录")
        if os.path.exists(self.output_dir + '/last.read.1.fastq'):
            pass
        self.logger.info("生成结果目录成功")

    def run(self):
        super(ChooseFastqTool, self).run()
        self.extract_fasta_name()
        self.compare_fasta_name()
        self.choose_seq()
        #self.cal_coverage()
        self.set_output()
        self.end()
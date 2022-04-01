# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import shutil
"""
Trimmomatic-0.35主要针对Illumina的fastq文件进行trim和去接头，认为接头只在3'-端出现。
由于Trimmomatic-0.35在获取adapter文件路径时，各种参数使用“:”隔开，故软件在windows上使用时避免使用“C:”等路径。
各个步骤之间存在先后顺序remove_adapter,run_head_cut,run_start_reserve,run_start_quality,run_end_quality,
run_slidingwindow,run_len_fliter,run_quality_fliter
"""


class FastqTrimAgent(Agent):
    """
    Trimmomatic
    version 0.35
    author shenghe
    last_modified:2016.1.7
    """
    def __init__(self, parent):
        super(FastqTrimAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq"},
            {"name": "fastq_re", "type": "infile", "format": "sequence.fastq"},
            {"name": "remove_adapter", "type": "bool", "default": False},
            {"name": "mode", "type": "string", "default": "simple"},  # palindrome
            {"name": "adapter_mode", "type": "string", "default": "TruSeq3-PE"},  # TruSeq3-PE,TruSeq2-PE,
            # simple_custom,palindrome_custom
            {"name": "pre_adaptor", "type": "string", "default": "TACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
            {"name": "pre_adaptor_re", "type": "string", "default": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"},
            {"name": "phred", "type": "string", "default": "Detect"},  # phred33, phred64
            {"name": "seedmismatch", "type": "int", "default": 1},
            {"name": "palindrome_matchscore", "type": "int", "default": 30},
            {"name": "simple_matchscore", "type": "int", "default": 10},
            {"name": "normal_adaptor", "type": "string", "default":
             "TACACTCTTTCCCTACACGACGCTCTTCCGATCT,GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"},
            # {"name": "normalmismatch", "type": "int", "default": -1},
            {"name": "run_head_cut", "type": "bool", "default": False},
            {"name": "head_cut", "type": "int", "default": 0},
            {"name": "run_start_reserve", "type": "bool", "default": False},
            {"name": "start_reserve", "type": "int", "default": 999},
            {"name": "run_start_quality", "type": "bool", "default": False},
            {"name": "start_quality", "type": "int", "default": 0},
            {"name": "run_end_quality", "type": "bool", "default": False},
            {"name": "end_quality", "type": "int", "default": 0},
            {"name": "run_slidingwindow", "type": "bool", "default": False},
            {"name": "window", "type": "int", "default": 3},
            {"name": "window_quality", "type": "int", "default": 0},
            {"name": "run_len_fliter", "type": "bool", "default": False},
            {"name": "len_fliter", "type": "int", "default": 0},
            {"name": "run_quality_fliter", "type": "bool", "default": False},
            {"name": "quality_fliter", "type": "int", "default": 0},
            {"name": "outfastq1", "type": "outfile", "format": "sequence.fastq"},
            {"name": "outfastq2", "type": "outfile", "format": "sequence.fastq"}
        ]
        self.add_option(options)
        self.step.add_steps('fastq_trim')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.fastq_trim.start()
        self.step.update()

    def step_end(self):
        self.step.fastq_trim.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if self.option('remove_adapter'):
            pass
        elif self.option('run_head_cut'):
            pass
        elif self.option('run_start_reserve'):
            pass
        elif self.option('run_start_quality'):
            pass
        elif self.option('run_end_quality'):
            pass
        elif self.option('run_slidingwindow'):
            pass
        elif self.option('run_len_fliter'):
            pass
        elif self.option('run_quality_fliter'):
            pass
        elif self.option('run_len_fliter'):
            pass
        else:
            raise OptionError('没有选择任何操作')
        if not self.option('fastq').is_set:
            raise OptionError('没有提供fastq文件')
        if self.option('remove_adapter'):
            if self.option('simple_matchscore') < 0 or self.option('simple_matchscore') > 30:
                raise OptionError('接头匹配分值设定不在正常范围内：%s' % self.option('simple_matchscore'))
            if self.option('mode') == 'palindrome':
                if self.option('seedmismatch') < 0:
                    raise OptionError('错误的种子错配容忍数:%s' % self.option('seedmismatch'))
                if 300 > self.option('palindrome_matchscore') > 0:
                    pass
                else:
                    raise OptionError('正反序列匹配的分值设定不在正常范围内：%s' % self.option('palindrome_matchscore'))
                if self.option('phred') not in ['phred33', 'phred64', 'Detect']:
                    raise OptionError('错误的phred设置：%s,如果不知道phred，可以使用Detect，由程序自己判断' % self.option('phred'))
                if not self.option('fastq_re').is_set:
                    raise OptionError('在模式palindrome下，必须提供前一个文件的反向序列文件')
                else:
                    if self.option('adapter_mode') not in ['simple_custom', 'palindrome_custom',
                                                           'TruSeq3-PE', 'TruSeq2-PE']:
                        raise OptionError('错误的adapter_mode模式选择：%s' % self.option('adapter_mode'))
                    if self.option('adapter_mode') == 'palindrome_custom':
                        if not self.option('pre_adaptor') or not self.option('pre_adaptor_re'):
                            raise OptionError('在palindrome模式下自定义adapter时，必须提供双向测序的两个prefix序列')
                        else:
                            if self.check_seq(self.option('pre_adaptor')):
                                pass
                            elif self.check_seq(self.option('pre_adaptor_re')):
                                pass
                            else:
                                raise OptionError('提供的adapters不是完全由ATGC组成:%s/%s' %
                                                  (self.option('pre_adaptor'), self.option('pre_adaptor_re')))
            else:
                if self.option('mode') != 'simple':
                    raise OptionError('去接头的模式选择错误：%s' % self.option('mode'))
        if self.option('run_start_quality'):
            if 64 > self.option('start_quality') > 0:
                pass
            else:
                raise OptionError('错误的序列起始质量控制设置:%s' % self.option('start_quality'))
        if self.option('run_end_quality'):
            if 64 > self.option('end_quality') > 0:
                pass
            else:
                raise OptionError('错误的序列结尾质量控制设置:%s' % self.option('end_quality'))
        if self.option('run_head_cut'):
            if self.option('head_cut') < 0:
                raise OptionError('序列从头去除数设置错误:%s' % self.option('head_cut'))
        if self.option('run_slidingwindow'):
            if self.option('window') < 1:
                raise OptionError('错误的质量检测窗口大小设置：%s' % self.option('window'))
            if 64 > self.option('window_quality') > 0:
                pass
            else:
                raise OptionError('错误的质量检测窗口质量设置：%s' % self.option('window_quality'))
        if self.option('run_quality_fliter'):
            if self.option('quality_fliter') < 0:
                raise OptionError('错误的序列过滤质量设定：%s' % self.option('quality_fliter'))
        if self.option('run_len_fliter'):
            if self.option('len_fliter') < 0:
                raise OptionError('错误的序列长度过滤设定：%s' % self.option('len_fliter'))

    def check_seq(self, seq):
        basecount = seq.count('A') + seq.count('T') + seq.count('G') + seq.count('C')
        if basecount == len(seq):
            return True
        else:
            return False

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 6
        self._memory = '5G'


class FastqTrimTool(Tool):

    def __init__(self, config):
        super(FastqTrimTool, self).__init__(config)
        self._version = '0.35'
        self.trimmomatic_path = ('program/sun_jdk1.8.0/bin/java -jar ' + self.config.SOFTWARE_DIR +
                                 '/bioinfo/seq/trimmomatic-0.36/trimmomatic-0.36.jar')
        # self.cutadapt_path = self.config.SOFTWARE_DIR + 'Python/bin/cutadapt'
        self.adapter = ''

    def run(self):
        """
        运行
        """
        super(FastqTrimTool, self).run()
        self.run_trimmomatic()

    def run_trimmomatic(self):
        """
        运行trimmomatic
        """
        cmd = ''
        if self.option('fastq_re').is_set:
            self.logger.info(self.work_dir)
            self.logger.info('forward1.fastq ')
            sequence1 = self.work_dir + '/forward1.fastq '
            sequence2 = self.work_dir + '/forward2.fastq '
            sequence3 = self.work_dir + '/reverse1.fastq '
            sequence4 = self.work_dir + '/reverse2.fastq '
            if self.option('remove_adapter'):
                if self.option('phred') != 'Detect':
                    cmd = self.trimmomatic_path + ' PE ' + '-' + self.option('phred') + ' '
                else:
                    cmd = self.trimmomatic_path + ' PE '
                adapter = self.get_adapeter()
                self.logger.info(self.option('seedmismatch'))
                self.logger.info(sequence1)
                self.logger.info(adapter)
                cmd = cmd + (self.option('fastq').prop['path'] + ' ' + self.option('fastq_re').prop['path'] +
                             ' ' + sequence1 + sequence2 + sequence3 + sequence4 + 'ILLUMINACLIP:' + adapter +
                             ':' + str(self.option('seedmismatch')) + ':' +
                             str(self.option('palindrome_matchscore')) + ':' + str(self.option('simple_matchscore')) +
                             ':1:true ')
            else:
                if self.option('phred') != 'Detect':
                    cmd = self.trimmomatic_path + ' PE ' + '-' + self.option('phred') + ' '
                else:
                    cmd = self.trimmomatic_path + ' PE '
                cmd = cmd + (self.option('fastq').prop['path'] + ' ' + self.option('fastq_re').prop['path'] +
                             ' ' + sequence1 + sequence2 + sequence3 + sequence4)
        else:
            sequence1 = self.work_dir + '/forward1.fastq '
            if self.option('phred') != 'Detect':
                cmd = cmd + self.trimmomatic_path + ' SE ' + '-' + self.option('phred') + ''
            else:
                cmd = self.trimmomatic_path + ' SE '
            if self.option('remove_adapter'):
                adapter = self.get_adapeter()
                cmd = cmd + (self.option('fastq').prop['path'] + ' ' + sequence1 + 'ILLUMINACLIP:' +
                             adapter + ':' + str(self.option('seedmismatch')) +
                             ':30:' + str(self.option('simple_matchscore')))
            else:
                cmd = cmd + self.option('fastq').prop['path'] + ' ' + sequence1
        if self.option('run_head_cut'):
            cmd = cmd + ' HEADCROP:' + str(self.option('head_cut'))
        if self.option('run_start_reserve'):
            cmd = cmd + ' CROP:' + str(self.option('start_reserve'))
        if self.option('run_start_quality'):
            cmd = cmd + ' LEADING:' + str(self.option('start_quality'))
        if self.option('run_end_quality'):
            cmd = cmd + ' TRAILING:' + str(self.option('end_quality'))
        if self.option('run_slidingwindow'):
            cmd = cmd + ' SLIDINGWINDOW:' + str(self.option('window')) + ':' + str(self.option('window_quality'))
        if self.option('run_len_fliter'):
            cmd = cmd + ' MINLEN:' + str(self.option('len_fliter'))
        if self.option('run_quality_fliter'):
            cmd = cmd + ' AVGQUAL:' + str(self.option('quality_fliter'))
        self.logger.info('运行Trimmomatic-0.36程序')
        self.logger.info(cmd)
        trimmomatic_command = self.add_command('trimmomatic', cmd)
        trimmomatic_command.run()
        self.wait()
        if trimmomatic_command.return_code == 0:
            self.logger.info('Trimmomatic-0.36正常结束')
            self.set_output()
            self.end()
        else:
            self.logger.info('程序返回值不为0：%s' % trimmomatic_command.return_code)
            self.set_error('返回值不为0：%s' % trimmomatic_command.return_code)

    def get_adapeter(self):
        """
        获取或者生成正确的adapter文件
        """
        self.logger.info(self.option('adapter_mode'))
        if self.option('adapter_mode') in ['simple_custom', 'palindrome_custom']:
            self.make_adapterfile()
            self.adapter = self.work_dir + '/adaptertemp.fasta'
        elif self.option('adapter_mode') == 'TruSeq3-PE':
            self.adapter = self.config.SOFTWARE_DIR + '/trim/Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa'
        elif self.option('adapter_mode') == 'TruSeq2-PE':
            self.adapter = self.config.SOFTWARE_DIR + '/trim/Trimmomatic-0.35/adapters/TruSeq2-PE.fa'
        else:
            self.set_error('adapter_mode 错误，无法选定正确的adapter文件！')
        return self.adapter

    def make_adapterfile(self):
        """
        生成一个adapter文件，依据使用的mode和提供的序列
        """
        tempadapter = open(self.work_dir + '/adaptertemp.fasta', 'w')
        if self.option('mode') == 'palindrome':
            templine = []
            templine.append('>PrefixPE/1\n')
            templine.append(self.option('pre_adaptor') + '\n')
            templine.append('>PrefixPE/2\n')
            templine.append(self.option('pre_adaptor_re') + '\n')
            templine.append('>PE1\n')
            templine.append(self.option('pre_adaptor') + '\n')
            templine.append('>PE1_rc\n')
            templine.append(self.comple_seq(self.option('pre_adaptor'), True) + '\n')
            templine.append('>PE2\n')
            templine.append(self.option('pre_adaptor_re') + '\n')
            templine.append('>PE2_rc\n')
            templine.append(self.comple_seq(self.option('pre_adaptor_re'), True) + '\n')
            tempadapter.writelines(templine)
        else:
            templine = []
            i = 1
            for one_adapter in self.option('normal_adaptor').split(','):
                templine.append('>PE' + str(i) + '\n')
                i += 1
                templine.append(one_adapter + '\n')
            tempadapter.writelines(templine)
        tempadapter.close()

    def comple_seq(self, sequence, reverse=False):
        """
        互补或者反向互补一条序列
        """
        basedict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        newseq = ''
        for i in sequence:
            newseq += basedict[i]
        if reverse:
            newseq = newseq[::-1]  # 把序列倒过来
        return newseq

    def set_output(self):
        """
        设置输出文件夹
        """
        self.clean_dir(self.output_dir)
        if self.option('fastq_re').is_set:
            os.link(self.work_dir + '/forward1.fastq', self.output_dir + '/forward1.fastq')
            os.link(self.work_dir + '/reverse1.fastq', self.output_dir + '/reverse1.fastq')
            self.option('outfastq1', self.output_dir + '/forward1.fastq')
            self.option('outfastq2', self.output_dir + '/reverse1.fastq')
        else:
            os.link(self.work_dir + '/forward1.fastq', self.output_dir + '/forward1.fastq')
            self.option('outfastq1', self.output_dir + '/forward1.fastq')

    def clean_dir(self, onedir):
        """
        清空一个文件夹
        """
        shutil.rmtree(onedir, ignore_errors=True)
        if not os.path.exists(self.output_dir):
            os.mkdir(onedir)

# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import shutil
import threading
from mbio.packages.sequence.standard_qualname import standard_qualname
from mbio.packages.sequence.fastq_to_fasta import convertfastq



class FastaTrimAgent(Agent):
    """
    cutadapt
    version 1.9.1
    author shenghe
    last_modified:2016.1.12
    """
    def __init__(self, parent):
        super(FastaTrimAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "qual", "type": "infile", "format": "sequence.qual"},
            {"name": "phred", "type": "string", "default": "phred33"},  # phred64
            {"name": "remove_adapter", "type": "bool", "default": False},
            {"name": "mismatch_rate", "type": "float", "default": 0.1},
            {"name": "border_minmatch", "type": "int", "default": 5},  # 边界最低匹配数
            {"name": "indel", "type": "bool", "default": True},
            {"name": "mode", "type": "string", "default": "5\'-end"},  # 3\'-end
            {"name": "adapter", "type": "string", "default":
             "TACACTCTTTCCCTACACGACGCTCTTCCGATCT,GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"},
            {"name": "run_head_cut", "type": "bool", "default": False},
            {"name": "head_cut", "type": "int", "default": 0},
            {"name": "run_end_cut", "type": "bool", "default": False},
            {"name": "end_cut", "type": "int", "default": 0},
            {"name": "run_start_quality", "type": "bool", "default": False},
            {"name": "start_quality", "type": "int", "default": 0},
            {"name": "run_end_N", "type": "bool", "default": False},
            {"name": "run_end_quality", "type": "bool", "default": False},
            {"name": "end_quality", "type": "int", "default": 0},
            {"name": "run_minlen_fliter", "type": "bool", "default": False},
            {"name": "minlen_fliter", "type": "int", "default": 0},
            {"name": "run_maxlen_fliter", "type": "bool", "default": False},
            {"name": "maxlen_fliter", "type": "int", "default": 999},
            {"name": "run_count_N", "type": "bool", "default": False},
            {"name": "max_N", "type": "float", "default": 10},
            {"name": "fastq_return", "type": "bool", "default": False},
            {"name": "outfastq", "type": "outfile", "format": "sequence.fastq"},
            {"name": "outfasta", "type": "outfile", "format": "sequence.fasta"},
            {"name": "outqual", "type": "outfile", "format": "sequence.qual"}
        ]
        self.add_option(options)
        self.step.add_steps('fasta_trim')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.fasta_trim.start()
        self.step.update()

    def step_end(self):
        self.step.fasta_trim.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if self.option('remove_adapter'):
            pass
        elif self.option('run_head_cut'):
            pass
        elif self.option('run_end_cut'):
            pass
        elif self.option('run_start_quality'):
            pass
        elif self.option('run_end_quality'):
            pass
        elif self.option('run_end_N'):
            pass
        elif self.option('run_minlen_fliter'):
            pass
        elif self.option('run_maxlen_fliter'):
            pass
        elif self.option('run_count_N'):
            pass
        else:
            raise OptionError('没有选择任何操作')
        if not self.option('fasta').is_set:
            raise OptionError('没有提供fasta文件')
        else:
            if os.path.splitext(self.option('fasta').prop['path'])[1] not in ['.fa', '.fasta',
                                                                              '.fna', '.csfasta', '.csfa']:
                raise OptionError('提供的fasta文件必须以fasta相关扩展名:.fa, .fasta, .fna, .csfasta, .csfa')
        if self.option('run_start_quality') or self.option('run_end_quality'):
            if self.option('qual').is_set:
                if os.path.splitext(self.option('qual').prop['path'])[1] != '.qual':
                    raise OptionError('提供的fasta文件的qual质量文件必须以.qual为扩展名:%s' %
                                      self.option('qual').prop['path'])
            else:
                raise OptionError('没有提供质量文件，无法进行质量过滤')
            if self.option('start_quality') < 0:
                raise OptionError('不在范围内的5\'端质量过滤设置：%s' % self.option('start_quality'))
            if self.option('end_quality') < 0:
                raise OptionError('不在范围内的3\'端质量过滤设置：%s' % self.option('end_quality'))
            if self.option('phred') not in ['phred33', 'phred64']:
                raise OptionError('phred设置错误：%s' % self.option('phred'))
        if self.option('remove_adapter'):
            if 0 <= self.option('mismatch_rate') < 1:
                pass
            else:
                raise OptionError('错配比率应在[0-1)的范围内：%s' % self.option('mismatch_rate'))
            if self.option('border_minmatch') < 3:
                raise OptionError('边界最小的碱基匹配数目过小：%s' % self.option('border_minmatch'))
            if self.option('mode') not in ['5\'-end', '3\'-end']:
                raise OptionError('去接头的模式必须为5端或者3端：%s' % self.option('mode'))
            adapter = self.option('adapter').split(',')
            adapter_just = [self.check_seq(i) for i in adapter]
            if False in adapter_just:
                raise OptionError('存在错误的接头序列：%s' % self.option('adapter'))
        if self.option('run_head_cut'):
            if self.option('head_cut') < 0:
                raise OptionError('序列从头去除数设置错误:%s' % self.option('head_cut'))
        if self.option('run_end_cut'):
            if self.option('end_cut') < 0:
                raise OptionError('错误的序列结尾去除设定：%s' % self.option('end_cut'))
        if self.option('run_start_quality'):
            if self.option('start_quality') < 0:
                raise OptionError('错误的序列起始质量过滤设定：%s' % self.option('start_quality'))
        if self.option('run_end_quality'):
            if self.option('end_quality') < 0:
                raise OptionError('错误的序列结尾质量过滤设定：%s' % self.option('end_quality'))
        if self.option('run_minlen_fliter'):
            if self.option('minlen_fliter') < 0:
                raise OptionError('错误的序列最小长度过滤设定：%s' % self.option('minlen_fliter'))
        if self.option('run_maxlen_fliter'):
            if self.option('maxlen_fliter') < 0:
                raise OptionError('错误的序列最大长度过滤设定：%s' % self.option('maxlen_fliter'))
        if self.option('run_count_N'):
            if self.option('max_N') < 0:
                raise OptionError('错误的最大N碱基设定：%s' % self.option('max_N'))

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
        self._memory = ''


class FastaTrimTool(Tool):

    def __init__(self, config):
        super(FastaTrimTool, self).__init__(config)
        self._version = '1.9.1'
        self.cutadapt_path = '/miniconda2/bin/python ' + self.config.SOFTWARE_DIR + '/miniconda2/bin/cutadapt '
        self.adapter = ''

    def run(self):
        """
        运行
        """
        super(FastaTrimTool, self).run()
        self.run_cutadapt()

    class result_thread(threading.Thread):
        def __init__(self, func, *argu, **kwargu):
            super(FastaTrimTool.result_thread, self).__init__()
            self.func = func
            self.argu = argu
            self.kwargu = kwargu
            self.result = None

        def run(self):
            self.result = self.func(*self.argu, **self.kwargu)

    def run_cutadapt(self):
        """
        运行cutadapt
        """
        cmd = self.cutadapt_path
        if self.option('qual').is_set:
            if self.option('remove_adapter'):
                if self.option('mode') == '5\'-end':
                    for i in self.option('adapter').split(','):
                        cmd = cmd + ' -g ' + i
                else:
                    for i in self.option('adapter').split(','):
                        cmd = cmd + ' -a ' + i
                cmd = cmd + ' -e ' + str(self.option('mismatch_rate'))
                cmd = cmd + ' -O ' + str(self.option('border_minmatch'))
                if not self.option('indel'):
                    cmd = cmd + ' --no-indels '
            if self.option('phred') == 'phred64':
                cmd = cmd + ' --quality-base 64 '
            if self.option('run_head_cut'):
                cmd = cmd + ' -u ' + str(self.option('head_cut'))
            if self.option('run_end_cut'):
                cmd = cmd + ' -u ' + '-' + str(self.option('head_cut'))
            if self.option('run_start_quality'):
                if self.option('run_end_quality'):
                    cmd = cmd + ' -q ' + str(self.option('start_quality')) + ',' + str(self.option('end_quality'))
                else:
                    cmd = cmd + ' -q ' + str(self.option('start_quality')) + ',0'
            else:
                if self.option('run_end_quality'):
                    cmd = cmd + ' -q ' + str(self.option('end_quality'))
            if self.option('run_end_N'):
                cmd = cmd + ' --trim-n '
            if self.option('run_minlen_fliter'):
                cmd = cmd + ' -m ' + str(self.option('minlen_fliter'))
            if self.option('run_maxlen_fliter'):
                cmd = cmd + ' -M ' + str(self.option('maxlen_fliter'))
            if self.option('run_count_N'):
                cmd = cmd + ' --max-n ' + str(self.option('max_N'))
            cmd = cmd + ' -o ' + self.work_dir + '/trimed.fastq'
            self.logger.info('核对并转换质量文件序列名称')
            qual_thread = FastaTrimTool.result_thread(standard_qualname, self.option('fasta').prop['path'],
                                                      self.option('qual').prop['path'], self.work_dir + '/new.qual')
            qual_thread.setDaemon(True)
            qual_thread.start()
            qual_thread.join()
            return_code = qual_thread.result
            if return_code[0] == 0:
                pass
            else:
                self.set_error('核对并转换质量文件序列名称出错，错误为：%s' % return_code[1])
            cmd = cmd + ' ' + self.option('fasta').prop['path'] + ' ' + self.work_dir + '/new.qual'
        else:
            if self.option('remove_adapter'):
                if self.option('mode') == '5\'-end':
                    for i in self.option('adapter').split(','):
                        cmd = cmd + ' -g ' + i
                else:
                    for i in self.option('adapter').split(','):
                        cmd = cmd + ' -a ' + i
                cmd = cmd + ' -e ' + str(self.option('mismatch_rate'))
                cmd = cmd + ' -O ' + str(self.option('border_minmatch'))
                if not self.option('indel'):
                    cmd = cmd + ' --no-indels '
            if self.option('run_head_cut'):
                cmd = cmd + ' -u ' + str(self.option('head_cut'))
            if self.option('run_end_cut'):
                cmd = cmd + ' -u ' + '-' + str(self.option('head_cut'))
            if self.option('run_end_N'):
                cmd = cmd + ' --trim-n '
            if self.option('run_minlen_fliter'):
                cmd = cmd + ' -m ' + str(self.option('minlen_fliter'))
            if self.option('run_maxlen_fliter'):
                cmd = cmd + ' -M ' + str(self.option('maxlen_fliter'))
            if self.option('run_count_N'):
                cmd = cmd + ' --max-n ' + str(self.option('max_N'))
            cmd = cmd + ' -o ' + self.work_dir + '/trimed.fasta'
            cmd = cmd + ' ' + self.option('fasta').prop['path']
        self.logger.info('运行cutadapt程序')
        self.logger.info(cmd)
        cutadapt_command = self.add_command('cutadapt', cmd)
        cutadapt_command.run()
        self.wait()
        if cutadapt_command.return_code == 0:
            self.logger.info('cutadapt_command正常结束')
            self.set_output()
            self.end()
        else:
            self.logger.info('程序返回值不为0：%s' % cutadapt_command.return_code)
            self.set_error('返回值不为0：%s' % cutadapt_command.return_code)

    def set_output(self):
        """
        设置输出文件夹
        """
        self.clean_dir(self.output_dir)
        if self.option('qual').is_set:
            if self.option('fastq_return'):
                os.link(self.work_dir + '/trimed.fastq', self.output_dir + '/trimed.fastq')
                self.option('outfastq', self.output_dir + '/trimed.fastq')
            else:
                convert_thread = FastaTrimTool.result_thread(convertfastq, self.work_dir + '/trimed.fastq',
                                                             self.output_dir + '/trimed.fasta',
                                                             self.output_dir + '/trimed.qual')
                convert_thread.setDaemon(True)
                convert_thread.start()
                convert_thread.join()
                return_code = convert_thread.result
                if not return_code:
                    self.logger.info('fastq转换fasta和qual正常')
                else:
                    self.set_error('fastq转换fasta和qual文件出错，错误为：%s' % return_code[1])
                self.option('outfasta', self.output_dir + '/trimed.fasta')
                self.option('outqual', self.output_dir + '/trimed.qual')
        else:
            os.link(self.work_dir + '/trimed.fasta', self.output_dir + '/trimed.fasta')
            self.option('outfasta', self.output_dir + '/trimed.fasta')

    def clean_dir(self, onedir):
        """
        清空一个文件夹
        """
        shutil.rmtree(onedir, ignore_errors=True)
        if not os.path.exists(self.output_dir):
            os.mkdir(onedir)

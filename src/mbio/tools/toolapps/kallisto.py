# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from Bio import SeqIO
import numpy as np
import shutil
import os


class KallistoAgent(Agent):
    '''
    宏转录组 kallisto计算基因表达量
    '''
    def __init__(self, parent):
        super(KallistoAgent, self).__init__(parent)
        options = [
            {'name': 'fq_type', 'type': 'string'}, # fastq类型 PE SE
            {'name': 'ref_fa', 'type': 'infile', 'format': 'sequence.fasta'}, ##参考组装拼接序列
            {'name': 'threads', 'type': 'int', 'default': 10}, ## cpu线程数
            {'name': 'fq_l', 'type': 'infile', 'format': 'sequence.fastq'}, ## 左端序列
            {'name': 'fq_r', 'type': 'infile', 'format': 'sequence.fastq'}, ## 右端序列
            {'name': 'fq_s', 'type': 'infile', 'format': 'sequence.fastq'}, ## 单端序列
            # {'name': 'strand_specific', 'type': 'bool', 'default': False},
            # {'name': 'stranded', 'type': 'string', 'default': None}, # rf (workflow forward) fr (workflow reverse)
            {'name': 'sample', 'type': 'string'}, ## 样本名称
            {'name': 'abundance_tsv', 'type': 'outfile', 'format': 'sequence.profile_table'}, #输出结果文件目录
        ]
        self.add_option(options)
        self.step.add_steps('kallisto')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.kallisto.start()
        self.step.update()

    def stepfinish(self):
        self.step.kallisto.finish()
        self.step.update()

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option('fq_type'):
            raise OptionError('必须设置测序类型：PE OR SE', code="34401101")
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError('测序类型不在所给范围内', code="34401102")
        if (self.option("fq_type") == "PE") and (not self.option("fq_r").is_set) and (not self.option("fq_l").is_set):
            raise OptionError("PE测序时需设置左端序列和右端序列输入文件", code="34401103")
        if (self.option("fq_type") == "SE") and (not self.option("fq_s").is_set):
            raise OptionError("SE测序时需设置序列输入文件", code="34401104")
        if not self.option("ref_fa").is_set:
            raise OptionError("需设置ref_fa输入文件", code="34401105")
        if (not self.option('sample')):
            raise OptionError('必须设置样本名称', code="34401106")
        return True

    def set_resource(self):
        """
        设置资源
        :return:
        """
        self._cpu = self.option('threads')
        if self.option('fq_type') == 'PE':
            size = os.path.getsize(self.option('fq_l').path) + os.path.getsize(self.option('fq_r').path)
            self._memory = '{}G'.format(int(float(size) / 1024 ** 3 + 20))
        elif self.option('fq_type') == 'SE':
            size = os.path.getsize(self.option('fq_s').path)
            self._memory = '{}G'.format(int(float(size) / 1024 ** 3 + 20))


    def end(self):
        """
        end 结束
        :return:
        """
        super(KallistoAgent, self).end()

class KallistoTool(Tool):
    """
    Tool kallisto计算基因表达量
    """
    def __init__(self, config):
        super(KallistoTool, self).__init__(config)
        self.kallisto = 'bioinfo/rna/kallisto_linux-v0.43.1/kallisto'
        self.index = os.path.join(self.work_dir, 'index')
        self.output = os.path.join(self.work_dir, self.option('sample'))
        self.abundance_tsv = os.path.join(self.output, 'abundance.tsv')

    def run(self):
        """
        运行
        :return:
        """
        self.logger.info("开始运行kallisto软件的tool！")
        super(KallistoTool, self).run()
        self.run_kallisto_index()
        if self.option('fq_type') == 'PE':
            self.run_kallisto_quant_pe()
        elif self.option('fq_type') == 'SE':
            self.run_kallisto_quant_se()
        self.set_output()
        self.end()

    def run_kallisto_index(self):
        """
        建索引
        :return:
        """
        self.logger.info("开始运行kallisto软件建立索引!")
        cmd = '{} index'.format(self.kallisto)
        cmd += ' -i {}'.format(self.index)
        cmd += ' {}'.format(self.option('ref_fa').path)
        cmd_name = 'run_kallisto_index'
        self.logger.info(cmd)
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('kallisto软件建立索引成功')
        elif command.return_code is None:
            self.logger.warn('kallisto软件建立索引失败, 再次运行')
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('kallisto软件建立索引成功')
            else:
                self.set_error('%skallisto软件建立索引失败', variables=(cmd_name), code="34401101")
        else:
            self.set_error('kallisto软件建立索引失败，返回码为：%s', variables=(command.return_code), code="34401102")

    def run_kallisto_quant_pe(self):
        """
        运行kallisto软件计算PE数据
        :return:
        """
        self.logger.info("开始运行kallisto软件计算表达量!")
        cmd = '{} quant'.format(self.kallisto)
        cmd += ' -i {}'.format(self.index)
        cmd += ' -o {}'.format(self.output)
        cmd += ' -t {}'.format(self.option('threads'))
        cmd += ' {} {}'.format(self.option('fq_l').prop['path'], self.option('fq_r').prop['path'])
        cmd_name = 'run_kallisto_quant_pe'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('kallisto软件计算成功')
        elif command.return_code is None:
            self.logger.warn('kallisto软件计算失败 {}, 再次运行')
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('kallisto软件计算成功')
            else:
                self.set_error('kallisto软件运行失败'.format(cmd_name), code="34401103")
        else:
            self.set_error('kallisto软件计算失败，返回码为：%s', variables=(command.return_code), code="34401104")

    def run_kallisto_quant_se(self):
        """
        运行kallisto软件计算SE数据
        注意 -l 和-s 参数 在计算single时是必须的
        :return:
        """
        self.logger.info("开始运行kallisto软件计算表达量!")
        kallisto_data = self.prepare_kallisto()
        cmd = '{} quant --single'.format(self.kallisto)
        cmd += ' -i {}'.format(self.index)
        cmd += ' -o {}'.format(self.output)
        cmd += ' -l {}'.format(kallisto_data['l'])
        cmd += ' -s {}'.format(kallisto_data['s'])
        cmd += ' -t {}'.format(self.option('threads'))
        cmd += ' {}'.format(self.option('fq_s').prop['path'])
        cmd_name = 'run_kallisto_quant_se'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('kallisto软件计算成功')
        elif command.return_code is None:
            self.logger.warn('kallisto软件计算失败 {}, 再次运行')
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('kallisto软件计算成功')
            else:
                self.set_error('kallisto软件运行失败'.format(cmd_name), code="34401105")
        else:
            self.set_error('kallisto软件计算失败，返回码为：%s', variables=(command.return_code), code="34401106")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info('开始设置结果文件目录')
        source = self.abundance_tsv
        link_name = os.path.join(self.output_dir, self.option("sample") + '.abundance.xls')
        if os.path.exists(link_name):
            os.remove(link_name)
        if os.path.exists(source):
            os.link(source, link_name)
            self.option('abundance_tsv').set_path(link_name)
        self.logger.info('设置结果文件目录完成')

    def prepare_kallisto(self):
        lengths = [len(seq_record) for seq_record in SeqIO.parse(self.option("fq_s").prop['path'], 'fastq')]
        kallisto = {'l': np.mean(lengths), 's': np.std(lengths)}
        return kallisto


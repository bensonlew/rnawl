# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import os


class SalmonAgent(Agent):
    '''
    宏转录组 salmon软件计算表达量
    '''
    def __init__(self, parent):
        super(SalmonAgent, self).__init__(parent)
        options = [
            {'name': 'fq_type', 'type': 'string', 'default': None},  ## fastq文件类型PE SE
            {'name': 'ref_fa', 'type': 'infile', 'format': 'sequence.fasta'}, ##参考组装拼接序列
            {'name': 'lib_type', 'type': 'string', 'default': 'IU'},##
            {'name': 'fq_l', 'type': 'infile', 'format': 'sequence.fastq'},##左端序列
            {'name': 'fq_r', 'type': 'infile', 'format': 'sequence.fastq'},##右端序列
            {'name': 'fq_s', 'type': 'infile', 'format': 'sequence.fastq'},##单端序列
            {'name': 'threads', 'type': 'int', 'default': 10},##cpu线程数
            {'name': 'sample', 'type': 'string'},##样品名称
        ]
        self.add_option(options)
        self.step.add_steps('salmon')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.salmon.start()
        self.step.update()

    def stepfinish(self):
        self.step.salmon.finish()
        self.step.update()

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option('fq_type'):
            raise OptionError('必须设置测序类型：PE OR SE', code="34400501")
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError('测序类型不在所给范围内', code="34400502")
        if (self.option("fq_type") == "PE") and (not self.option("fq_r").is_set) and (not self.option("fq_l").is_set):
            raise OptionError("PE测序时需设置左端序列和右端序列输入文件", code="34400503")
        if (self.option("fq_type") == "SE") and (not self.option("fq_s").is_set):
            raise OptionError("SE测序时需设置序列输入文件", code="34400504")
        if not self.option("ref_fa").is_set:
            raise OptionError("需设置ref_fa输入文件", code="34400505")
        if (not self.option('sample')):
            raise OptionError('必须设置样本名称', code="34400506")
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
        运行结束
        :return:
        """
        super(SalmonAgent, self).end()

class SalmonTool(Tool):
    """
    Tool salmon软件计算表达量
    """
    def __init__(self, config):
        super(SalmonTool, self).__init__(config)
        self.salmon = 'bioinfo/rna/Salmon-0.8.2_linux_x86_64/bin/salmon'
        self.index = os.path.join(self.work_dir, 'index')
        self.output = os.path.join(self.work_dir, '{}'.format(self.option('sample')))
        self.quant_sf = os.path.join(self.output, 'quant.sf')
        self.quant_genes_sf = os.path.join(self.output, 'quant.genes.sf')

    def run(self):
        """
        运行
        :return:
        """
        self.logger.info("开始运行tool！")
        super(SalmonTool, self).run()
        self.run_salmon_index()
        if self.option('fq_type') == 'PE':
            self.run_salmon_quant_pe()
        elif self.option('fq_type') == 'SE':
            self.run_salmon_quant_se()
        self.set_output()
        self.end()

    def run_salmon_index(self):
        """
        建立索引
        :return:
        """
        self.logger.info("开始建立索引!")
        cmd = '{} index'.format(self.salmon)
        cmd += ' -t {}'.format(self.option('ref_fa').path)
        cmd += ' -i {}'.format(self.index)
        cmd += ' -p {}'.format(self.option('threads'))
        cmd += ' --type quasi'
        cmd_name = 'run_salmon_index'
        self.logger.info(cmd)
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('salmon软件建立索引成功')
        elif command.return_code is None:
            self.logger.warn('salmon软件建立索引失败, 再次运行')
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('salmon软件建立索引成功')
            else:
                self.set_error('%ssalmon软件建立索引失败', variables=(cmd_name), code="34400501")
        else:
            self.set_error('salmon软件建立索引失败，返回码为：%s', variables=(command.return_code), code="34400502")

    def run_salmon_quant_pe(self):
        """
        运行salmon软件计算PE数据
        :return:
        """
        self.logger.info("开始用salmon软件计算表达量!")
        cmd = '{} quant --meta '.format(self.salmon)
        cmd += ' -i {}'.format(self.index)
        cmd += ' -l {}'.format(self.option('lib_type'))
        cmd += ' -1 {}'.format(self.option('fq_l').prop['path'])
        cmd += ' -2 {}'.format(self.option('fq_r').prop['path'])
        cmd += ' -o {}'.format(self.output)
        cmd += ' -p {}'.format(self.option('threads'))
        cmd_name = 'run_salmon_quant_pe'
        self.logger.info(cmd)
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('salmon软件计算成功')
        elif command.return_code is None:
            self.logger.warn('salmon软件计算失败, 再次运行')
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('salmon软件计算成功')
            else:
                self.set_error('%ssalmon软件计算失败', variables=(cmd_name), code="34400503")
        else:
            self.set_error('salmon软件计算失败，返回码为：%s', variables=(command.return_code), code="34400504")

    def run_salmon_quant_se(self):
        """
        运行salmon软件计算SE数据
        :return:
        """
        self.logger.info("开始用salmon软件计算表达量!")
        cmd = '{} quant --gcBias'.format(self.salmon)
        cmd += ' -i {}'.format(self.index)
        cmd += ' -l {}'.format(self.option('lib_type'))
        cmd += ' -r {}'.format(self.option('unmated_reads').path)
        cmd += ' -o {}'.format(self.output)
        cmd += ' -p {}'.format(self.option('threads'))
        cmd += ' -g {}'.format(self.option('t2g').path)
        cmd_name = 'run_salmon_quant_se'
        self.logger.info(cmd)
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('salmon软件计算成功')
        elif command.return_code is None:
            self.logger.warn('salmon软件计算失败, 再次运行')
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('salmon软件计算成功')
            else:
                self.set_error('%ssalmon软件计算失败', variables=(cmd_name), code="34400505")
        else:
            self.set_error('salmon软件计算失败，返回码为：%s', variables=(command.return_code), code="34400506")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info('开始设置结果文件目录')
        source = self.quant_sf
        link_name = os.path.join(self.output_dir, self.option("sample") + '.quant.xls')
        if os.path.exists(link_name):
            os.remove(link_name)
        if os.path.exists(source):
            os.link(source, link_name)
        source_genes = self.quant_genes_sf
        link_name_genes = os.path.join(self.output_dir, self.option("sample") + '.quant.genes.xls')
        if os.path.exists(link_name_genes):
            os.remove(link_name_genes)
        if os.path.exists(source_genes):
            os.link(source_genes, link_name_genes)
        self.logger.info("设置结果文件目录完成")

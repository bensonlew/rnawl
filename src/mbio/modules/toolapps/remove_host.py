# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
## @20200226

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.toolapps.common import link_dir


class RemoveHostModule(Module):
    """
    小工具：宏转录组bwa去宿主
    """
    def __init__(self, work_id):
        super(RemoveHostModule, self).__init__(work_id)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq序列文件夹
            {'name': 'ref_database', 'type': 'string', 'default': ''},  # 宿主参考序列库中对应的物种名，eg：
            {'name': 'second_ref', 'type': 'string', 'default': ''}, #参考数据库二级下拉框数据
            {'name': 'ref_undefined', "type": 'infile', 'format': 'sequence.fasta_dir'}, # 自定义上传的fq序列文件夹
            {'name': 'ref_undefined_name', 'type': 'string', 'default': 'undefined'},  # 自定义参考宿主名称，适应页面参数
        ]
        self.add_option(options)
        self.step.add_steps('rm_host')
        self.sequence = self.add_module('metagenome.reads_unzip')
        self.rm_host = self.add_module('meta.qc.bwa_remove_host')

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option('in_fastq').is_set:
            raise OptionError('请输入原始序列文件夹或者质控优质序列文件夹！', code="24400401")
        if self.option('ref_database') == '' and not self.option('ref_undefined').is_set:
            raise OptionError('需输入参考数据库或参考序列', code="24400402")
        if self.option('ref_database') not in ['', 'Custom'] and self.option('ref_undefined').is_set:
            raise OptionError('去宿主不可同时提供参考数据库及参考序列', code="24400403")

    def run_sequence(self):
        """
        对数据进行解压，上传序列为压缩格式
        :return:
        """
        opts = {
            'fastq_dir': self.option('in_fastq'),
        }
        self.sequence.set_options(opts)
        self.sequence.run()

    def run_rm_host(self):
        """
        去宿主
        :return:
        """
        ref_database = self.option('ref_database') + "," + self.option('second_ref')
        opts = {
            'fq_type': 'PE',
            'ref_database': ref_database,
            'ref_undefined': self.option('ref_undefined'),
        }
        if self.option('ref_database') == 'Custom':
            opts['ref_database'] = ""
        if self.option("in_fastq").is_set:
            opts['fastq_dir'] = self.sequence.output_dir + '/data'
        self.rm_host.set_options(opts)
        self.rm_host.run()

    def run(self):
        """
        运行
        :return:
        """
        super(RemoveHostModule, self).run()
        self.sequence.on("end", self.run_rm_host)
        self.rm_host.on("end", self.set_output)
        self.run_sequence()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("开始将结果文件复制到结果文件下")
        remove_host = os.path.join(self.output_dir, "romove_host")
        if os.path.exists(remove_host):
            shutil.rmtree(remove_host)
        else:
            os.mkdir(remove_host)
        link_dir(self.rm_host.output_dir, remove_host)
        self.logger.info("链接结果文件完成")
        self.end()

    def end(self):
        """
        运行结束
        :return:
        """
        result_dir = self.add_upload_dir(os.path.join(self.output_dir, 'romove_host'))
        self.logger.info("开始上传结果文件")
        result_dir.add_relpath_rules([
            [".", "", "bwa去宿主结果文件目录"],
        ])
        super(RemoveHostModule, self).end()
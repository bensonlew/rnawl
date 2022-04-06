# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class CircosPreAgent(Agent):
    """

    author: zouxuan
    last_modify: 20180402
    """

    def __init__(self, parent):
        super(CircosPreAgent, self).__init__(parent)
        options = [
            {"name": "assemble", "type": "infile", "format": "sequence.fasta"},  # 序列文件
            {"name": "location", "type": "string", "default": "Scaffold"},  # 位置信息
            {"name": "trna", "type": "infile", "format": "gene_structure.gff3"},  # trna.gff文件
            {"name": "rrna", "type": "infile", "format": "gene_structure.gff3"},  # rrna.gff文件
            {"name": "gene", "type": "infile", "format": "gene_structure.gff3"},  # gene.gff文件
            {"name": "anno_cog", "type": "infile", "format": "sequence.profile_table"}  # cog注释结果文件
        ]
        self.add_option(options)

    def check_options(self):
        # if not self.option("sequence").is_set:
        #     raise OptionError("必须设置输入序列文件")
        # if not self.option("table").is_set:
        #     raise OptionError("必须设置待修改的二维表格")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '3G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(CircosPreAgent, self).end()


class CircosPreTool(Tool):
    def __init__(self, config):
        super(CircosPreTool, self).__init__(config)
        self._version = "1.0"
        self.gc_script_path = self.config.PACKAGE_DIR + '/bacgenome/gc_count_and_skew.pl'
        self.ncrna_script_path = self.config.PACKAGE_DIR + '/bacgenome/ncRNA_legend.pl'
        self.cog_script_path = self.config.PACKAGE_DIR + '/bacgenome/gene_cog_type.py'
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        self.python_path = '/miniconda2/bin/python'

    def run(self):
        """
        运行
        :return:
        """
        super(CircosPreTool, self).run()
        self.option('assemble').merge_all_seq(self.work_dir + '/' + self.option('location'))
        self.get_karyotype()
        self.gc()
        self.cog()
        if not self.option('location').startswith('p'):
            self.ncrna()
        self.end()

    def get_karyotype(self):
        self.option('assemble').get_info()
        base = self.option('assemble').prop['bases']
        if self.option('location'):
            self.location = self.option('location')
        else:
            self.location = 'Scaffold'
        with open(self.work_dir + '/karyotype.txt', 'w') as wf, open(self.work_dir + '/temp.txt', 'w') as wt:
            wf.write("chr - " + self.location + " 1 1 " + base + " deepskyblue")
            wt.write(self.location + ' 1 ' + base)
        if os.path.exists(self.output_dir + '/karyotype.txt'):
            os.remove(self.output_dir + '/karyotype.txt')
        os.link(self.work_dir + '/karyotype.txt', self.output_dir + '/karyotype.txt')
        if os.path.exists(self.output_dir + '/temp.txt'):
            os.remove(self.output_dir + '/temp.txt')
        os.link(self.work_dir + '/temp.txt', self.output_dir + '/temp.txt')

    def gc(self):
        if self.option('location').startswith('p'):
            cmd = '{} {} {} 100 50'.format(self.perl_path, self.gc_script_path,
                                             self.work_dir + '/' + self.option('location') + '.fasta')
        else:
            cmd = '{} {} {} 1000 500'.format(self.perl_path, self.gc_script_path,
                                             self.work_dir + '/' + self.option('location') + '.fasta')
        command = self.add_command('gc_file', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gc文件生成成功")
        else:
            self.set_error("gc文件生成失败", code="31400901")
            self.set_error("gc文件生成失败", code="31400902")
        if os.path.exists(self.output_dir + '/negative_gc_count.txt'):
            os.remove(self.output_dir + '/negative_gc_count.txt')
        os.link(self.work_dir + '/negative_gc_count.txt', self.output_dir + '/negative_gc_count.txt')
        if os.path.exists(self.output_dir + '/negative_gc_skew.txt'):
            os.remove(self.output_dir + '/negative_gc_skew.txt')
        os.link(self.work_dir + '/negative_gc_skew.txt', self.output_dir + '/negative_gc_skew.txt')
        if os.path.exists(self.output_dir + '/positive_gc_count.txt'):
            os.remove(self.output_dir + '/positive_gc_count.txt')
        os.link(self.work_dir + '/positive_gc_count.txt', self.output_dir + '/positive_gc_count.txt')
        if os.path.exists(self.output_dir + '/positive_gc_skew.txt'):
            os.remove(self.output_dir + '/positive_gc_skew.txt')
        os.link(self.work_dir + '/positive_gc_skew.txt', self.output_dir + '/positive_gc_skew.txt')

    def ncrna(self):
        cmd = '{} {} -t {} -r {} -l {}'.format(self.perl_path, self.ncrna_script_path, self.option('trna').prop['path'],
                                               self.option('rrna').prop['path'], self.option('location'))
        command = self.add_command('ncrna', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("ncrna文件生成成功")
        else:
            self.set_error("ncrna文件生成失败", code="31400903")
            self.set_error("ncrna文件生成失败", code="31400904")
        if os.path.exists(self.output_dir + '/ncRNA.txt'):
            os.remove(self.output_dir + '/ncRNA.txt')
        os.link(self.work_dir + '/ncRNA.txt', self.output_dir + '/ncRNA.txt')

    def cog(self):
        cmd = '{} {} -g {} -c {} -l {}'.format(self.python_path, self.cog_script_path, self.option('gene').prop['path'],
                                               self.option('anno_cog').prop['path'], self.option('location'))
        command = self.add_command('cog', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cog文件生成成功")
        else:
            self.set_error("cog文件生成失败", code="31400905")
            self.set_error("cog文件生成失败", code="31400906")
        if os.path.exists(self.output_dir + '/sense_strand_cog.txt'):
            os.remove(self.output_dir + '/sense_strand_cog.txt')
        os.link(self.work_dir + '/sense_strand_cog.txt', self.output_dir + '/sense_strand_cog.txt')
        if os.path.exists(self.output_dir + '/antisense_strand_cog.txt'):
            os.remove(self.output_dir + '/antisense_strand_cog.txt')
        os.link(self.work_dir + '/antisense_strand_cog.txt', self.output_dir + '/antisense_strand_cog.txt')

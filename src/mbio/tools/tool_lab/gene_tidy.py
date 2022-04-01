# -*- coding: utf-8 -*-

from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
from biocluster.core.exceptions import OptionError
import subprocess


class GeneTidyAgent(Agent):
    """
    progigal小工具基因预测
    """

    def __init__(self, parent):
        super(GeneTidyAgent, self).__init__(parent)
        options = [
            {"name": "prodigal", "type": "infile", "format": "gene_structure.gff3"},
            {"name": "genome", "type": "infile", "format":"sequence.fasta"},
            {'name': 'sample_name', "type": "string","default":"out"},  # 样本名
            {"name": "orf_prefix", "type" : "string", "default": "ORF"},
            {"name": "trans_code", "type" : "string","default":"11"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("prodigal").is_set:
            raise OptionError("必须有一个软件的预测结果！")
        if not self.option('genome').is_set:
            raise OptionError("必须输入基因组的fna文件")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '2G'

    def end(self):
        super(GeneTidyAgent, self).end()

class GeneTidyTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(GeneTidyTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.script_tidy = self.config.PACKAGE_DIR + '/bacgenome/ '
        self.python_path = '/program/Python/bin/python'
        self.transeq = "/bioinfo/seq/EMBOSS-6.6.0/emboss/transeq"
        self.sample_name = self.option("sample_name")
        self.trans_table = self.config.SOFTWARE_DIR + '/database/bacgenome/trans_table/code.dic'
        self.trans_py = self.config.PACKAGE_DIR + '/bacgenome/trans_table.py'
        self.LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/bioinfo/seq/EMBOSS-6.6.0/lib"
        self.set_environ(LD_LIBRARY_PATH=self.LD_LIBRARY_PATH)

        self.prodigal_gff =  self.option("prodigal").prop['path']
        self.genemark_gff = ''
        self.glimmer_gff = ''
        self.trna_gff = ''
        self.rrna_gff = ''
        self.genome_fasta = self.option('genome').prop['path']

    def run_tidy(self):
        cmd = '{} {}/predict_result_remove_overlap.py'.format(self.python_path,self.config.PACKAGE_DIR +'/tool_lab')
        if self.prodigal_gff:
            cmd += ' -p {} '.format(self.prodigal_gff)
        if self.genemark_gff:
            cmd += ' -m {} '.format(self.genemark_gff)
        if self.glimmer_gff:
            cmd += ' -g {} '.format(self.glimmer_gff)
        if self.trna_gff:
            cmd += ' -t {} '.format(self.trna_gff)
        if self.rrna_gff:
            cmd += ' -r {} '.format(self.rrna_gff)

        command = self.add_command("run_tidy", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_tidy运行完成")
        else:
            self.set_error("run_tidy运行出错!")


    def run_fna2ffn(self):
        predict_file = self.work_dir + "/pos.predict"
        cmd = '{} {}glimmer2ffn.pl {} {} {} 0 {}'.format(self.perl_path, self.config.PACKAGE_DIR + "/tool_lab/", self.genome_fasta,
                                                         predict_file, self.option("orf_prefix"), self.sample_name)
        command = self.add_command("run_fna2ffn", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取fnn运行完成")
            if os.path.getsize(self.work_dir + "/" + self.sample_name + ".fnn"):
                self.run_transeq()
            else:
                self.end()
        else:
            self.set_error("提取fnn运行出错!")

    def run_transeq(self):
        nul_seq = self.work_dir + "/" + self.sample_name + ".fnn"
        port_seq = self.work_dir + "/" + self.sample_name + ".faa"
        cmd = '{} -trim -table 11 -sequence {} -outseq {}'.format(self.transeq, nul_seq, port_seq)
        command = self.add_command("transeq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("翻译蛋白序列运行完成")
            self.run_get_information()
        else:
            self.set_error("翻译蛋白序列运行出错!")

    def run_get_information(self):
        nul_seq = self.work_dir + "/" + self.sample_name + ".fnn"
        port_seq = self.work_dir + "/" + self.sample_name + ".faa"
        cmd = '{} {}trim_gene_info.pl {} {} {} {} {}'.format(self.perl_path, self.perl_script, self.genome_fasta,
                                                             nul_seq, port_seq, self.option("orf_prefix"),
                                                             self.output_dir)
        command = self.add_command("get_information", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("统计预测结果运行完成")
        else:
            self.set_error("统计预测结果运行出错!")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        self.logger.info("设置结果目录")
        gff = '{}.predict.gff'.format(self.sample_name)
        if os.path.exists(self.output_dir + '/' + gff):
            os.remove(self.output_dir + '/' + gff)
        os.link(self.work_dir + '/' + gff, self.output_dir + '/' + gff)

        self.logger.info("设置结果目录成功")

    def run(self):
        super(GeneTidyTool, self).run()
        self.run_tidy()
        self.run_fna2ffn()
        self.set_output()
        self.end()

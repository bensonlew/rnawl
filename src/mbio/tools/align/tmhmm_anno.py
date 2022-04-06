# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import subprocess


class TmhmmAnnoAgent(Agent):
    """
    跨膜蛋白注释结果 v1.0
    author: zouxuan
    last_modify: 20180203
    """

    def __init__(self, parent):
        super(TmhmmAnnoAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample", "type": "string"} , # 样品名称
            {"name": "analysis", "type": "string", "default":""} #complete or  uncomplete
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置输入序列文件", code="31102401")
        if not self.option("sample"):
            raise OptionError("必须设置样品名称", code="31102402")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '3G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(TmhmmAnnoAgent, self).end()


class TmhmmAnnoTool(Tool):
    def __init__(self, config):
        super(TmhmmAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.tmhmm_path = self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/tmhmm-2.0c/bin/tmhmm'
        self.tmhmm_sh = 'bioinfo//Genomic/Sofware/tmhmm-2.0c/bin/tmhmm.sh'
        self.perl_path = 'program/perl-5.24.0/bin/perl '
        self.script_path = self.config.PACKAGE_DIR + '/annotation/signalp_tmhmmtxt2xls.pl'
        self.python_script_path = self.config.PACKAGE_DIR + '/bacgenome/add_gene_info.py '
        self.python_path = '/miniconda2/bin/python'

    def run(self):
        """
        运行
        :return:
        """
        super(TmhmmAnnoTool, self).run()
        self.tmhmm_predict()
        self.end()

    def tmhmm_predict(self):
        linkfile = self.work_dir + '/' + os.path.basename(self.option('query').prop['path'])
        if os.path.exists(linkfile):
            os.remove(linkfile)
        os.link(self.option('query').prop['path'], linkfile)
        cmd = '{} {} {} {}'.format(self.tmhmm_sh, self.tmhmm_path, linkfile, self.work_dir + '/tmhmm.txt')
        command = self.add_command('predict', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("tmhmm预测成功")
        else:
            self.set_error("tmhmm预测失败", code="31102401")
            self.set_error("tmhmm预测失败", code="31102402")
        cmd1 = '{} {} -faa {} -{} {} -outfile {}'.format(self.perl_path, self.script_path,
                                                         self.option('query').prop['path'],
                                                         'tmhmm', self.work_dir + '/tmhmm.txt',
                                                         self.work_dir + '/' + 'tmhmm_stat.xls')
        command1 = self.add_command('stat', cmd1).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("tmhmm注释结果统计成功")
        else:
            self.set_error("tmhmm注释结果统计失败", code="31102403")
            self.set_error("tmhmm注释结果统计失败", code="31102404")
        cmd2 = "awk -F \"\t\" \'{if((NR >13 && $7==\"YES\")||($1==\"Gene ID\")){print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$8}}\' " \
               + self.work_dir + "/tmhmm_stat.xls > " + self.output_dir + "/" + "gene_tmhmm_anno.xls"
        try:
            subprocess.check_output(cmd2, shell=True)
            self.logger.info("数据整理成功")
        except subprocess.CalledProcessError:
            self.set_error("数据整理失败", code="31102405")

        if self.option("analysis") in ["complete"]:
            split = True
        else:
            split = False
        anno_file = self.output_dir + "/" + "gene_tmhmm_anno.xls"
        lines_counts = 0
        if os.path.exists(anno_file):
            with open(anno_file, "r") as f:
                lines_counts = len(f.readlines())
        if lines_counts >= 2:
            cmd3 = '{} {} -i {} -f {} -s {} -d {} '.format(self.python_path, self.python_script_path,
                                                           self.output_dir + "/" + "gene_tmhmm_anno.xls",
                                                           self.option('query').prop['path'],
                                                           self.option('sample'), split)
            command3 = self.add_command('get_file', cmd3).run()
            self.wait(command3)
            if command3.return_code == 0:
                self.logger.info("文件生成成功")
            else:
                self.set_error("文件生成失败", code="31102406")
                self.set_error("文件生成失败", code="31102407")
        self.rename("gene_tmhmm_anno.xls", self.option('sample') + "_tmhmm_anno.xls")

    def rename(self, old, new):
        if os.path.exists(self.output_dir + '/' + old):
            os.rename(self.output_dir + '/' + old, self.output_dir + '/' + new)

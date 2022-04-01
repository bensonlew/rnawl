# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os
import pandas as pd


class ContributeAgent(Agent):
    """
    宏基因组贡献度分析
    author: shaohua.yuan
    last_modify: haidong.gu
    """

    def __init__(self, parent):
        super(ContributeAgent, self).__init__(parent)
        options = [
            {"name": "taxon_file", "type": "infile", "format": "sequence.profile_table"},
            {"name": "function_file", "type": "infile", 'format': "sequence.profile_table"},
            {"name": "gene_profile", "type": "infile", 'format': "sequence.profile_table"},
            {"name": "tax_level", "type": "string"},
            {"name": "fun_level", "type": "string"},
            {"name": "top_tax", "type": "int", "default": 10},
            {"name": "top_fun", "type": "int", "default": 10},
            {"name": "tax_fun_abu", "type": "outfile", 'format': "sequence.profile_table"},
            {"name": "fun_tax_abu", "type": "outfile", 'format': "sequence.profile_table"},
        ]
        self.add_option(options)
        self.result_name = ''

    def check_options(self):
        if not self.option("taxon_file").is_set:
            raise OptionError("必须设置输入物种注释文件", code="32701001")
        if not self.option("function_file").is_set:
            raise OptionError("必须设置输入功能注释文件", code="32701002")
        if not self.option("gene_profile").is_set:
            raise OptionError("必须设置输入基因丰度文件", code="32701003")
        return True

    def set_resource(self):
        self._cpu = 2
        # memory = 5 + 10 * self._rerun_time  # 每次重运行加10G  by guhaidong @ 20180416
        self._memory = '5G'  # 改回 by guhaidong @ 20180427
        # self._memory = "%sG" % memory

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(ContributeAgent, self).end()


class ContributeTool(Tool):
    def __init__(self, config):
        super(ContributeTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        self.script = self.config.PACKAGE_DIR + '/metagenomic/scripts/tax_fun.pl'

    def run(self):
        """
        运行
        :return:
        """
        super(ContributeTool, self).run()
        self.run_tax_fun()
        self.set_output()
        self.end()

    def run_tax_fun(self):
        tax_f = self.option("taxon_file").prop["path"]
        fun_f = self.option("function_file").prop["path"]
        gene_profile = self.option("gene_profile").prop["path"]
        cmd = '{} {} -tax {} -fun {} -tl {} -fl {} -p {} -Tt {} -Tf {} -o {}'.format(self.perl_path, self.script, tax_f,
                                                                                    fun_f, self.option("tax_level"),
                                                                                    self.option("fun_level"),gene_profile,
                                                                                    self.option("top_tax"),
                                                                                    self.option("top_fun"),
                                                                                    self.output_dir)
        command = self.add_command("contribute", cmd, ignore_error=True).run()
        self.logger.info(cmd)
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("contribute select succeed")
        elif command.return_code == -9:  # 检测超出内存，重新运行并自动添加内存 by guhaidong @ 20180416
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("contribute select failed", code="32701001")
            raise Exception("contribute select failed")


    def set_output(self):
        tax_fun_file = os.path.join(self.output_dir,"Taxon_function_abundance.xls")
        fun_tax_file = os.path.join(self.output_dir,"Function_taxon_abundance.xls")
        table1 = pd.read_table(tax_fun_file, sep='\t', header=0)
        if len(table1) > 0:
            self.option("tax_fun_abu", tax_fun_file)
        else:
            self.set_error("筛选的TOP物种与功能间无共同基因！", code="32701002")
        table2 = pd.read_table(fun_tax_file, sep='\t', header=0)
        if len(table1) > 0:
            self.option("tax_fun_abu", fun_tax_file)
        else:
            self.set_error("筛选的TOP物种与功能间无共同基因！", code="32701003")
        self.logger.info("contribute output right")



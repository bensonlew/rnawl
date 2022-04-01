#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re


class RenameGbkAgent(Agent):
    """
    用于修改gbk文件生成
    version 1.0
    author: gaohao
    last_modify: 2018.04.27
    """

    def __init__(self, parent):
        super(RenameGbkAgent, self).__init__(parent)
        options = [
            {"name": "gbk", "type": "infile", "format": "gene_structure.gbk"},  #
            {"name": "des", "type": "string"},  #
            {"name": "source", "type": "string"},  #
            {"name": "organism", "type": "string"},  #
            {"name": "author", "type": "string"},  #
            {"name": "title", "type": "string"},  #
            {"name": "journal", "type": "string"},  #
            {"name": "sample_name", "type": "string"},
            {"name": "seq_type", "type": "string"},
        ]
        self.add_option(options)

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('gbk').is_set:
            raise OptionError("请设置基因组原始gbk文件！", code="31402801")
        if not self.option('des'):
            raise OptionError("请设置gbk文件的des信息！", code="31402802")
        if not self.option('source'):
            raise OptionError("请设置gbk文件的source信息！", code="31402803")
        if not self.option('organism'):
            raise OptionError("请设置gbk文件的organism信息！", code="31402804")
        if not self.option('author'):
            raise OptionError("请设置gbk文件的author信息！", code="31402805")
        if not self.option('title'):
            raise OptionError("请设置gbk文件的title信息！", code="31402806")
        if not self.option('journal'):
            raise OptionError("请设置gbk文件的journal信息！", code="31402807")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '2G'

    def end(self):
        super(RenameGbkAgent, self).end()


class RenameGbkTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(RenameGbkTool, self).__init__(config)
        self.gbk = self.option('gbk').prop['path']
        self.des =self.option('des')
        self.source = self.option('source')
        self.organism = self.option('organism')
        self.author = self.option('author')
        self.title = self.option('title')
        self.journal = self.option('journal')
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        if not self.option('seq_type'):
            self.out = self.work_dir + '/' + self.option('sample_name') + '.gbk'
        else:
            self.out = self.work_dir + '/' + self.option('sample_name') + '_' + self.option('seq_type') + '.gbk'




    def run_gbk(self):
        cmd = "{} {}modify_gbk.pl -i {} -des '{}' -source '{}' -organism '{}' -author '{}' -title '{}' -journal '{}' -o {}".format(self.perl_path, self.perl_script,self.gbk ,self.des,self.source,self.organism ,self.author ,self.title,self.journal,self.out)
        self.logger.info(cmd)
        self.logger.info("开始运行run_gbk")
        command = self.add_command("run_gbk", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_gbk完成")
        else:
            self.set_error("运行run_gbk运行出错!", code="31402801")



    def set_output(self):
        self.logger.info('set output')
        if not self.option('seq_type'):
            if os.path.exists(self.output_dir + '/' + self.option('sample_name') + '.gbk'):
                os.remove(self.output_dir + '/' + self.option('sample_name') + '.gbk')
            os.link(self.out,self.output_dir + '/' + self.option('sample_name') + '.gbk')
        else:
            self.out = self.work_dir + '/' + self.option('sample_name') + '_' + self.option('seq_type') + '.gbk'
            if os.path.exists(self.output_dir + '/' + self.option('sample_name') + '_' + self.option('seq_type') + '.gbk'):
                os.remove(self.output_dir + '/' + self.option('sample_name') + '_' + self.option('seq_type') + '.gbk')
            os.link(self.out,self.output_dir + '/' + self.option('sample_name') + '_' + self.option('seq_type') + '.gbk')


    def run(self):
        """
        运行
        """
        super(RenameGbkTool, self).run()
        self.run_gbk()
        self.set_output()
        self.end()



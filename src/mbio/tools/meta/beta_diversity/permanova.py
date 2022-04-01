# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError
import re  # add re by guhaidong 20171025
from mbio.packages.taxon.mask_taxon import mask_taxon  # add by guhaidong 20171025


class PermanovaAgent(Agent):
    """
    脚本permanova.pl
    version v1.0
    author: zouxuan
    last_modified:2017.10.19
    """

    def __init__(self, parent):
        super(PermanovaAgent, self).__init__(parent)
        options = [
            {"name": "abu_table", "type": "infile",
             "format": "meta.otu.otu_table, meta.otu.tax_summary_dir, toolapps.table"},
            # modify by zhouxuan 20170623 小工具的模块是指定toolapps.table这个文件类型的
            # {"name": "level", "type": "string", "default": "otu"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "permutations", "type": "int", "default": 999},
            {"name": "dis_method", "type": "string", "default": "bray"},
            {"name": "binary", "type": "string", "default": "false"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('abu_table').is_set:
            raise OptionError('必须提供丰度表', code="32703101")
        if not self.option('envtable').is_set:
            raise OptionError('必须提供环境因子表', code="32703102")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
        ])
        # print self.get_upload_files()
        super(PermanovaAgent, self).end()


class PermanovaTool(Tool):
    def __init__(self, config):
        super(PermanovaTool, self).__init__(config)
        self._version = '1.0'
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.script1_path = self.config.PACKAGE_DIR + "/beta_diversity/permanova.pl"
        self.R_path = 'program/R-3.3.1/bin/Rscript'

    def run_permanova(self):
        otu_table = self.option("abu_table").prop['path']
        self.name_to_name = mask_taxon(otu_table, self.work_dir + "/tmp_mask_otu.xls")  # add by guhaidong 20171025
        otu_table = self.work_dir + '/tmp_mask_otu.xls'  # add by guhaidong 20171025
        cmd = self.perl_path
        cmd += ' %s -i %s -env %s -o %s  -pe %s -d_method %s -binary %s' % (
            self.script1_path, otu_table, self.option("envtable").prop['path'],
            self.work_dir + '/permanova.xls', self.option("permutations"), self.option("dis_method"),
            self.option("binary"))
        self.logger.info(cmd)
        command = self.add_command('cmd', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("permanova的r文件生成成功")
        else:
            self.set_error("permanova的r文件生成失败", code="32703101")
            raise Exception("permanova的r文件生成失败")
        cmd_1 = '%s %s' % (self.R_path, self.work_dir + "/cmd.r")
        self.logger.info(cmd_1)
        command1 = self.add_command('cmd_1', cmd_1)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("R程序计算permanova成功")
        else:
            self.set_error("R程序计算permanova失败", code="32703102")
            raise Exception("R程序计算permanova失败")

    def dashrepl(self, matchobj):
        """
        add func by guhaidong 20171025
        """
        return self.name_to_name[matchobj.groups()[0]]

    def add_taxon(self, old_result, taxon_result):
        """
        add func by guhaidong 20171025
        description: 将旧注释的名称，根据词典替换成新注释名称
        """
        with open(old_result, "r") as f, open(taxon_result, "w") as w:
            # w.write(old_result)
            for i in f.readlines():
                # line = i.strip()
                new_line = re.sub(r"(name\d+)", self.dashrepl, i)
                w.write(new_line)

    def run(self):
        super(PermanovaTool, self).run()
        self.run_permanova()
        self.set_output()

    def set_output(self):
        linksites = os.path.join(self.output_dir, 'permanova.xls')
        if os.path.exists(linksites):
            os.remove(linksites)
        os.link(self.work_dir + '/permanova.xls', linksites)
        self.end()

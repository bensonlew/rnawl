# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError
import shutil


class UnigeneProfileAgent(Agent):
    """
    version v1.0
    author: zouxuan
    last modified:2017.9.11
    """

    def __init__(self, parent):
        super(UnigeneProfileAgent, self).__init__(parent)
        options = [
            {"name": "map_dir", "type": "string", "default": ""},  # map结果
            {"name": "fafile", "type": "infile", "format": "sequence.fasta"},  # 非冗余基因集fasta文件
            {"name": "insertsize", "type": "infile", "format": "sample.insertsize_table"},  # 插入片段文件
            {"name": "rpkm_abundance", "type": "outfile", "format": "sequence.profile_table"},  # RPKM丰度
            {"name": "reads_abundance", "type": "outfile", "format": "sequence.profile_table"},  # reads丰度
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fafile").is_set:
            raise OptionError("必须提供非冗余基因集")
        if not self.option("map_dir"):
            raise OptionError("必须提供map结果")
        if not self.option("insertsize").is_set:
            raise OptionError("必须插入片段文件")

    def set_resource(self):
        self._cpu = 1
        if os.path.getsize(self.option("fafile").prop['path']) / 100000000 < 6:
            self._memory = '10G'
        else:
            self._memory = str(os.path.getsize(self.option("fafile").prop['path']) / 100000000 + 35) + 'G'
            if os.path.getsize(self.option("fafile").prop['path']) / 100000000 + 35 > 125:
                self._memory = '125G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["gene_profile.reads_number.xls", "xls", "reads丰度表"],
            ["gene_profile.reads_percent.xls", "xls", "reads相对丰度表"],
            ["gene_profile.base_percent.xls", "xls", "reads数除以基因长度丰度表"],
            ["gene_profile.base_number.xls", "xls", "reads数除以基因长度相对丰度表"],
            ["gene_profile.RPKM.xls", "xls", "RPKM丰度表"],
            ["gene_profile.RPKM_percent.xls", "xls", "RPKM相对丰度表"],
        ])
        super(UnigeneProfileAgent, self).end()


class UnigeneProfileTool(Tool):
    def __init__(self, config):
        super(UnigeneProfileTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.script1_path = self.config.PACKAGE_DIR + "/statistical/prepare_profile.pl"
        self.script2_path = self.config.PACKAGE_DIR + "/statistical/gene_profile.pl"

    def mapinfo(self):
        cmd1 = '%s %s %s %s %s %s' % (self.perl_path, self.script1_path, self.option("map_dir"),
                                      self.option("fafile").prop['path'], self.option("insertsize").prop['path'],
                                      self.work_dir)
        self.logger.info(cmd1)
        self.logger.info("生成soap_info")
        command1 = self.add_command("saop_info", cmd1)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("soap_info succeed")
        else:
            self.set_error("soap_info failed")
            raise Exception("soap_info failed")

    def get_profile(self):
        if not os.path.exists(self.work_dir + '/gene_profile'):
            os.mkdir(self.work_dir + '/gene_profile')
        cmd2 = '%s %s -i1 %s -i2 %s -o %s' % (
            self.perl_path, self.script2_path, self.option("fafile").prop['path'], self.work_dir + '/soap.info',
            self.work_dir + '/gene_profile')
        self.logger.info("生成profile文件")
        command2 = self.add_command("profile", cmd2)
        command2.run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("profile succeed")
        else:
            self.set_error("profile failed")
            raise Exception("profile failed")

    def set_output(self):
        self.linkdir(self.work_dir + '/gene_profile', "")
        self.option('rpkm_abundance', os.path.join(self.output_dir, "RPKM.xls"))
        self.option('reads_abundance', os.path.join(self.output_dir, "reads_number.xls"))
        self.end()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        """
        # allfiles = os.listdir(dirpath)
        allfiles = ["gene.uniGeneset.fa.length.txt", "reads_number.xls", "reads_number_relative.xls",
                    "reads_length_ratio.xls", "reads_length_ratio_relative.xls", "RPKM.xls", "TPM.xls",
                    "top100_reads_number_relative.xls", "top100_reads_number.xls"]
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def run(self):
        super(UnigeneProfileTool, self).run()
        self.mapinfo()
        self.get_profile()
        self.set_output()

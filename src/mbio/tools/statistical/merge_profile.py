# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError
import shutil


class MergeProfileAgent(Agent):
    """
    version v1.0
    author: zouxuan
    last modified:2018.0228
    """

    def __init__(self, parent):
        super(MergeProfileAgent, self).__init__(parent)
        options = [
            {"name": "profile_dir", "type": "string", "default": ""},  # 单样品丰度结果文件夹
            {"name": "fafile", "type": "infile", "format": "sequence.fasta"},  # 非冗余基因集fasta文件
            {"name": "samples", "type": "string"},  # 样品名称以逗号分隔
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("profile_dir"):
            raise OptionError("必须提供单样品丰度结果文件夹", code="34100701")
        if not self.option("fafile").is_set:
            raise OptionError("必须提供非冗余基因集", code="34100702")
        if not self.option("samples"):
            raise OptionError("必须提供样品名", code="34100703")

    def set_resource(self):
        self._cpu = 1
        if os.path.getsize(self.option("fafile").prop['path']) / 100000000 < 6:
            self._memory = '10G'
        else:
            self._memory = str(os.path.getsize(self.option("fafile").prop['path']) / 100000000 + 10) + 'G'
            if os.path.getsize(self.option("fafile").prop['path']) / 100000000 + 10 > 125:
                self._memory = '125G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(MergeProfileAgent, self).end()


class MergeProfileTool(Tool):
    def __init__(self, config):
        super(MergeProfileTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.script1_path = self.config.PACKAGE_DIR + "/statistical/merge_profile.pl"

    def merge(self):
        cmd1 = '%s %s -i %s -s %s -f %s -o %s' % (self.perl_path, self.script1_path, self.option("profile_dir"),
                                      self.option("samples"),self.option("fafile").prop['path'],
                                      self.work_dir+'/gene_profile')
        self.logger.info(cmd1)
        self.logger.info("开始merge丰度文件")
        command1 = self.add_command("merge", cmd1)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("merge profile succeed")
        else:
            self.set_error("merge profile failed", code="34100701")
            self.set_error("merge profile failed", code="34100702")

    def set_output(self):
        self.linkdir(self.work_dir + '/gene_profile', "")
        self.end()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        """
        # allfiles = os.listdir(dirpath)
        allfiles = ["gene.uniGeneset.fa.length.txt", "reads_number.xls", "reads_number_relative.xls",
                    "reads_length_ratio.xls", "reads_length_ratio_relative.xls", "RPKM.xls", "TPM.xls",
                    "top100_reads_number_relative.xls", "top100_reads_number.xls","PPM.xls"]
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            if os.path.exists(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])

    def run(self):
        super(MergeProfileTool, self).run()
        self.merge()
        self.set_output()

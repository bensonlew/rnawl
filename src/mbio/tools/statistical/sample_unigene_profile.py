# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError
import shutil


class SampleUnigeneProfileAgent(Agent):
    """
    version v1.0
    author: zouxuan
    last modified:20180604
    """

    def __init__(self, parent):
        super(SampleUnigeneProfileAgent, self).__init__(parent)
        options = [
            {"name": "map_result", "type": "string", "default": ""},  # map结果
            {"name": "fafile", "type": "infile", "format": "sequence.fasta"},  # 非冗余基因集fasta文件
            {"name": "insertsize", "type": "int"},  # 插入片段长度
            {"name": "sample", "type": "string"},  # 样品名称
            {"name": "out_dir", "type": "string"}  # 存放文件夹
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fafile").is_set:
            raise OptionError("必须提供非冗余基因集", code="34101101")
        if not self.option("map_result"):
            raise OptionError("必须提供map结果", code="34101102")
        if not self.option("insertsize"):
            raise OptionError("必须提供插入片段长度", code="34101103")
        if not self.option("sample"):
            raise OptionError("必须提供样品名称", code="34101104")

    def set_resource(self):
        self._cpu = 1
        if os.path.getsize(self.option("fafile").prop['path']) / 100000000 < 6:
            self._memory = '20G'  # 10G内存改为15G by guhaidong @ 20180413
        else:
            self._memory = str(os.path.getsize(self.option("fafile").prop['path']) / 100000000 + 20) + 'G'  # 基础内存改为20G guhaidong @ 20180413
            if os.path.getsize(self.option("fafile").prop['path']) / 100000000 + 20 > 125:
                self._memory = '125G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(SampleUnigeneProfileAgent, self).end()


class SampleUnigeneProfileTool(Tool):
    def __init__(self, config):
        super(SampleUnigeneProfileTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.script2_path = self.config.PACKAGE_DIR + "/statistical/sample_gene_profile.pl"

    def get_profile(self):
        if not os.path.exists(self.work_dir + '/gene_profile'):
            os.mkdir(self.work_dir + '/gene_profile')
        cmd2 = '%s %s -i1 %s -i2 %s -i3 %s -n %s -o %s' % (
            self.perl_path, self.script2_path, self.option("fafile").prop['path'], self.option("insertsize"),self.option("map_result"),self.option("sample"),
            self.work_dir + '/gene_profile')
        cmd2 += " -ppm T"  # add by shaohua.yuan 20181107, add ppm method
        self.logger.info("生成profile文件")
        command2 = self.add_command("profile", cmd2, ignore_error=True)  # add ignore_err by guhaidong @ 20180614
        command2.run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("profile succeed")
        elif command2.return_code == -9:
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180604
        else:
            self.set_error("profile failed", code="34101101")
            self.set_error("profile failed", code="34101102")

    def set_output(self):
        self.linkdir(self.work_dir + '/gene_profile', "")
        self.end()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.option("out_dir"), dirname)
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
        super(SampleUnigeneProfileTool, self).run()
        self.get_profile()
        self.set_output()

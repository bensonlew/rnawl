#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re,shutil
from mbio.packages.metagbin.common_function import link_dir



class IslandDimobAgent(Agent):
    """
    用于扫描图的gbk文件生成
    version 1.0
    author: gaohao
    last_modify: 2020.07.02
    """

    def __init__(self, parent):
        super(IslandDimobAgent, self).__init__(parent)
        options = [
            {"name": "dir", "type": "string"},  #输入文件的dir
            {"name": "sample_name", "type": "string"}, #样品名称
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self.list =[]


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('dir'):
            raise OptionError("请设置基因组文件夹不存在")
        else:
            for file in os.listdir(self.option('dir')):
                if not file.endswith((".ptt", '.faa', '.ffn')):
                    raise OptionError("请提供正确的文件,{}不是.ptt，.faa，.ffn结尾的文件".format(file))

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(IslandDimobAgent, self).end()


class IslandDimobTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(IslandDimobTool, self).__init__(config)
        self.dir =self.option('dir')
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.dimob = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/islandpath_dimob/Dimob.pl"

    def run_dimob(self):
        if os.path.exists(self.work_dir + "/temp"):
            shutil.rmtree(self.work_dir + "/temp")
        link_dir(self.dir, self.work_dir + "/temp")
        cmd = "{} {} {} {} {} ".format(self.perl_path, self.dimob, self.work_dir + "/temp", self.option("sample_name"), self.work_dir + "/" + self.option("sample_name") + ".island_dimob.xls")
        self.logger.info(cmd)
        command = self.add_command("run_dimob", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            if os.path.getsize(self.work_dir + "/" + self.option("sample_name") + ".island_dimob.xls") >0:
                self.get_location(self.work_dir + "/" + self.option("sample_name") + ".island_dimob.xls", self.output_dir + "/" + self.option("sample_name") + ".island_dimob.xls", self.option("sample_name"))
            self.logger.info("运行run_dimob完成")
        else:
            self.set_error("运行run_dimob运行出错!")

    def set_output(self):
        if os.path.exists(self.output_dir + "/" + self.option("sample_name") + ".island_dimob.xls"):
            self.option('out').set_path(self.output_dir + "/" + self.option("sample_name") + ".island_dimob.xls")

    def run(self):
        """
        运行
        """
        super(IslandDimobTool, self).run()
        self.run_dimob()
        self.set_output()
        self.end()

    def get_location(self, file, out, location):
        with open(file,"r") as f,open(out,"w") as g:
            lines = f.readlines()
            for line in lines:
                lin = line.strip().split("\t")
                g.write("{}\t{}\t{}\t{}\n".format(lin[0],location,lin[1], lin[2]))
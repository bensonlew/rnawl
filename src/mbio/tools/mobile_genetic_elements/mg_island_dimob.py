#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re,shutil
from mbio.packages.metagbin.common_function import link_dir



class MgIslandDimobAgent(Agent):
    """
    用于扫描图的gbk文件生成
    version 1.0
    author: gaohao
    last_modify: 2020.07.02
    """

    def __init__(self, parent):
        super(MgIslandDimobAgent, self).__init__(parent)
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
            raise OptionError("请设置文件夹不存在")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(MgIslandDimobAgent, self).end()


class MgIslandDimobTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(MgIslandDimobTool, self).__init__(config)
        self.dir =self.option('dir')
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.dimob = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/islandpath_dimob/Dimob.pl"

    def run_dimob(self):
        if os.path.exists(self.work_dir + "/" + self.option("sample_name")):
            shutil.rmtree(self.work_dir + "/" + self.option("sample_name"))
        os.mkdir(self.work_dir + "/" + self.option("sample_name"))
        n = 1
        for i in os.listdir(self.dir):
            if os.path.exists(self.work_dir + "/temp"+str(n)):
                shutil.rmtree(self.work_dir + "/temp"+str(n))
            link_dir(self.dir+"/" + i, self.work_dir + "/temp"+str(n))
            if os.path.getsize(self.work_dir + "/temp"+str(n)+"/"+i+".faa") > 0:
                cmd = "{} {} {} {} {} ".format(self.perl_path, self.dimob, self.work_dir + "/temp" + str(n),
                                               i, self.work_dir + "/" + self.option(
                        "sample_name") + "/" + "island_dimob." + str(n) + ".xls")
                self.logger.info(cmd)
                command = self.add_command("run_dimob" + str(n), cmd)
                command.run()
                self.wait(command)
                if command.return_code == 0:
                    self.logger.info("运行run_dimob完成")
                else:
                    self.set_error("运行run_dimob运行出错!")
                n += 1


    def diomb_stat(self):
        with open(self.output_dir + "/" + self.option("sample_name") + ".island_dimob.xls", "w") as g:
            n=1
            for i in os.listdir(self.work_dir + "/" + self.option("sample_name")):
                location = i.split(".")[0]
                if os.path.getsize(self.work_dir + "/" + self.option("sample_name") + "/" + i) > 0:
                    with open(self.work_dir + "/" + self.option("sample_name") + "/" + i, "r") as f:
                        lines = f.readlines()
                        for line in lines:
                            lin = line.strip().split("\t")
                            g.write("{}\t{}\t{}\t{}\n".format("GI"+str(n), location, lin[1], lin[2]))
                            n += 1

    def set_output(self):
        if os.path.exists(self.output_dir + "/" + self.option("sample_name") + ".island_dimob.xls"):
            self.option('out').set_path(self.output_dir + "/" + self.option("sample_name") + ".island_dimob.xls")

    def run(self):
        """
        运行
        """
        super(MgIslandDimobTool, self).run()
        self.run_dimob()
        self.diomb_stat()
        self.set_output()
        self.end()

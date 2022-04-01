#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
import shutil


class MgIslandIslanderAgent(Agent):
    """
    使用islander预测基因组岛
    version 1.0
    author: gaohao
    last_modify: 2020.07.02
    """

    def __init__(self, parent):
        super(MgIslandIslanderAgent, self).__init__(parent)
        options = [
            {"name": "fa_dir", "type": "infile", "format": "sequence.fasta_dir"},  #
            {"name": "out", "type": "outfile", "format": "bacgenome.island"},
        ]
        self.add_option(options)
        self.list =[]


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('fa_dir').is_set:
            raise OptionError("请设置fa_dir目录！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(MgIslandIslanderAgent, self).end()


class MgIslandIslanderTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(MgIslandIslanderTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl5path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/tRNAscan-SE-2.0/lib/tRNAscan-SE:" + self.config.SOFTWARE_DIR + "/program/perl-5.24.0/lib"
        self.set_environ(PERL5LIB=self.perl5path)
        self.islander = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/Islander_software/Islander.pl"


    def run_islander(self):
        if os.path.exists(self.work_dir + "/result"):
            shutil.rmtree(self.work_dir + "/result")
        os.mkdir(self.work_dir + "/result")
        n = 1
        for i in os.listdir(self.option("fa_dir").prop['path']):
            sample = os.path.basename(self.option("fa_dir").prop['path'] + "/" + i)
            if os.path.exists(self.work_dir + "/" + sample):
                os.remove(self.work_dir + "/" + sample)
            os.link(self.option("fa_dir").prop['path'] + "/" + i, self.work_dir + "/" + sample)
            fa = sample.split(".fna")[0]
            cmd = "{} {} {} --verbose --translate --trna --annotate --reisland --table 11 --nocheck".format(
                self.perl_path, self.islander, fa)
            self.logger.info(cmd)
            path = 'run_islander'+str(n)
            self.logger.info("开始运行%s " % path)
            command = self.add_command(path, cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                if os.path.exists(self.work_dir + '/all_final_arrayed_islands.gff'):
                    os.rename(self.work_dir + '/all_final_arrayed_islands.gff', self.work_dir + "/result/" + fa + ".all_final_arrayed_islands.gff")
                self.logger.info("运行%s完成" % path)
            else:
                self.set_error("运行%s运行出错!" % path)
            n += 1

    def get_stat(self):
        with open(self.output_dir + "/islander.xls", "w") as g:
            for module in os.listdir(self.work_dir + "/result"):
                    if self.get_num(self.work_dir + "/result/" + module) >= 2:
                        with open(self.work_dir + "/result/" + module, "r") as f:
                            lines = f.readlines()
                            for line in lines[1:]:
                                lin = line.strip().split("\t")
                                g.write("{}\t{}\t{}\t{}\n".format(lin[0], lin[1], lin[3], lin[4]))
    def set_output(self):
        if os.path.exists(self.output_dir + "/islander.xls") and os.path.getsize(self.output_dir + "/islander.xls") > 0:
            self.option('out').set_path(self.output_dir + "/islander.xls")

    def run(self):
        """
        运行
        """
        super(MgIslandIslanderTool, self).run()
        self.run_islander()
        self.get_stat()
        self.set_output()
        self.end()

    def get_num(self, file):
        with open(file, "r") as f:
            lines =f.readlines()
            num = len(lines)
        return num
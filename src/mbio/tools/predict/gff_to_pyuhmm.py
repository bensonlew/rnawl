# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'
# version 1.0
# last_modify: 2018.06.04

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.ref_rna.trans_step import step_count  # 对fasta序列做长度分布图输入表


class GffToPyuhmmAgent(Agent):
    """
    生成pyu.hmm文件
    """

    def __init__(self, parent):
        super(GffToPyuhmmAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"} ,# 组装拼接好的scaffold文件
            {"name": "maker_gff", "type": "string", "default": ""},
            {"name": "cegma_gff","type":"string","default":""},
            {"name": "type","type":"string","default":"cegma"}  #value : cegma or maker
        ]

        self.add_option(options)

    def check_options(self):
        if self.option("type") == "cegma":
            if not self.option("cegma_gff"):
                raise OptionError("必须设置参数cegma_gff", code="33300401")
            if not self.option("fasta").is_set:
                raise OptionError("必须设置参数fasta", code="33300402")
        else:
            if not self.option("maker_gff"):
                raise OptionError("必须设置参数maker_gff", code="33300403")

    def set_resource(self):
        self._cpu = 1
        self._memory = '1G'

    def end(self):
        super(GffToPyuhmmAgent, self).end()


class GffToPyuhmmTool(Tool):
    def __init__(self, config):
        super(GffToPyuhmmTool, self).__init__(config)
        self.maker2zff = "bioinfo/Genomic/Sofware/maker/bin/maker2zff"
        self.cegma2zff = "bioinfo/Genomic/Sofware/maker/bin/cegma2zff"
        self.fathom = "bioinfo/Genomic/Sofware/SNAP-master/SNAP-master/fathom"
        self.forge = "bioinfo/Genomic/Sofware/SNAP-master/SNAP-master/forge"
        self.hmm_assemble = "/bioinfo/Genomic/Sofware/SNAP-master/SNAP-master/hmm-assembler.pl"
        self.pyu_hmm = ""


    def run_cmd1(self):
        if self.option("type") == 'cegma':
            self.genome_fasta = self.option("fasta").prop['path']
            cmd1 = "{} {} {}".format(self.cegma2zff,self.option("cegma_gff"),self.genome_fasta)
        else:
            perl_lib="{0}/program/perl-5.24.0/lib/site_perl/5.24.0:{0}/bioinfo/Genomic/Sofware/maker-2.31.9/maker/lib:{0}/bioinfo/Genomic/Sofware/maker-2.31.9/maker/perl/lib:{0}/program/perl-5.24.0/lib/site_perl/5.24.0/x86_64-linux-thread-multi".format(self.config.SOFTWARE_DIR)
            self.set_environ(PERL5LIB=perl_lib)
            cmd1 = "{} -n -l 90 {}".format(self.maker2zff,self.option("maker_gff"))
        self.logger.info("开始运行{}".format(cmd1))
        command = self.add_command("run_cmd1", cmd1).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd1))
        else:
            self.set_error("%s运行出错!", variables=(cmd1), code="33300401")

    def run_cmd2(self):
        cmd2 = "{} -categorize 1000 genome.ann genome.dna".format(self.fathom)
        self.logger.info("开始运行{}".format(cmd2))
        command = self.add_command("run_cmd2", cmd2).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd2))
        else:
            self.set_error("%s运行出错!", variables=(cmd2), code="33300402")

    def run_cmd3(self):
        cmd3 = "{} -export 1000 -plus uni.ann uni.dna".format(self.fathom)
        self.logger.info("开始运行{}".format(cmd3))
        command = self.add_command("run_cmd3", cmd3).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd3))
        else:
            self.set_error("%s运行出错!", variables=(cmd3), code="33300403")


    def run_cmd4(self):
        cmd4 = "{}  export.ann export.dna".format(self.forge)
        self.logger.info("开始运行{}".format(cmd4))
        command = self.add_command("run_cmd4", cmd4).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd4))
        else:
            self.set_error("%s运行出错!", variables=(cmd4), code="33300404")


    def run_cmd5(self):
        cmd5 = self.config.SOFTWARE_DIR + "{} pyu . > pyu.hmm".format(self.hmm_assemble)
        self.logger.info("开始运行{}".format(cmd5))
        desc = os.system(cmd5)
        if desc == 0:
            self.logger.info("{}运行完成".format(cmd5))
            self.pyu_hmm = os.path.join(self.work_dir,"pyu.hmm")
            if os.path.exists(self.pyu_hmm):
                self.logger.info("生成{}".format(self.pyu_hmm))
            else:
                self.logger.info("没有生成{}".format(self.pyu_hmm))
        else :
            self.logger.info("{}运行失败".format(cmd5))

    def set_output(self):
        new_path = os.path.join(self.output_dir, "pyu.hmm")
        if os.path.exists(new_path):
            os.remove(new_path)
        os.link(self.pyu_hmm,new_path)



    def run(self):
        super(GffToPyuhmmTool, self).run()
        self.run_cmd1()
        self.run_cmd2()
        self.run_cmd3()
        self.run_cmd4()
        self.run_cmd5()
        self.set_output()
        self.end()

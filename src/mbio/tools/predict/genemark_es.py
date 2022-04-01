# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'
# version 1.0
# last_modify: 2018.06.04

import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import shutil


class GenemarkEsAgent(Agent):
    """
    Genemark-ES 进行基因预测或训练
    """

    def __init__(self, parent):
        super(GenemarkEsAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"}, # 组装拼接好的scaffold文件
            {"name": "genemark_options", "type": "string","default":"--ES "}
        ]

        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome", code="33300301")

    def set_resource(self):
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(GenemarkEsAgent, self).end()


class GenemarkEsTool(Tool):
    def __init__(self, config):
        super(GenemarkEsTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.sample_name = (self.genome_fasta).split('/')[-1].rsplit('.')[0]
        self.perl_path = "/program/perl-5.24.0/bin/perl"

        self.genemark_path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/gm_et_linux_64/gmes_petap/gmes_petap.pl"
        self.genemark_out_gtf = ""
        self.genemark_out_hmm = ""


    def run_genemark(self):
        tmp_list = ["data","info","run"]
        mark_out = [os.path.join(self.work_dir, i) for i in tmp_list]
        for rm in mark_out:
            if os.path.exists(rm):
                shutil.rmtree(rm)
        os.system("cp ~/.gm_key ./")
        path_str="{}/program/perl-5.24.0/bin".format(self.config.SOFTWARE_DIR)
        self.set_environ(PATH=path_str)

        perl_lib="{0}/program/perl-5.24.0/lib/site_perl/5.24.0:{0}/program/perl-5.24.0/lib:$PERL5LIB".format(self.config.SOFTWARE_DIR)
        self.set_environ(PERL5LIB=perl_lib)


        cmd = '{} {} {} --sequence {}'.format(self.perl_path, self.genemark_path,self.option("genemark_options"), self.genome_fasta)
        command = self.add_command("run_genemark", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("genemark运行完成")
            self.genemark_out_gtf = os.path.join(self.work_dir,"genemark.gtf")
            self.genemark_out_hmm = os.path.join(self.work_dir,"output/gmhmm.mod")  #genemark 软件 运行会生成一个output文件夹，所以gmhmm.mod 无需在set_output中再做处理

        else:
            self.set_error("genemark运行出错!", code="33300301")

    def set_output(self):
        gtf_new_path = os.path.join(self.output_dir,'genemark.gtf')
        if os.path.exists(gtf_new_path):
            os.remove(gtf_new_path)
        os.link(self.genemark_out_gtf,gtf_new_path)



    def run(self):
        super(GenemarkEsTool, self).run()
        self.run_genemark()
        self.set_output()
        self.end()

# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'
# version 1.0
# last_modify: 2018.06.05

import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

import shutil

class AugustusAgent(Agent):
    """
    Augustus 进行基因预测
    """

    def __init__(self, parent):
        super(AugustusAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"}, # 组装拼接好的scaffold文件
            {"name": "species", "type": "string", "default": ""},  #augustus软件运行参数
            {"name": "trainingset", "type": "string", "default": ""},  #augustus软件运行参数
            {"name": "hints","type": "string", "default": ""}, #augustus软件运行参数
            {"name": "estali","type": "string", "default": ""}, #augustus软件运行参数
            {"name": "cdna","type": "string", "default": ""},  #augustus软件运行参数
            {"name": "pasa","type": "string", "default": "0"}, #augustus软件运行参数 ,设置有这参数时值为非0
            {"name": "pasapolyAhints","type": "string", "default": "0"}, #augustus软件运行参数 ，设置有这参数时值为非0
            {"name": "augustus_other_options", "type": "string", "default": " -v --useexisting"}   #augustus软件运行参数
        ]

        self.add_option(options)
        self.essential_options = ["genome","species"]


    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数%s", variables=(self.essential_options[0]), code="33300101")
        if not self.option("species"):
            raise OptionError("必须设置参数%s", variables=(self.essential_options[1]), code="33300102")


    def set_resource(self):
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(AugustusAgent, self).end()


class AugustusTool(Tool):
    def __init__(self, config):
        super(AugustusTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.augustus_out_gff = ""
        self.augustus_out_faa = ""

    def run_augustus(self):
        auto_aug = os.path.join(self.work_dir, "autoAug")
        if os.path.exists(auto_aug):
            shutil.rmtree(auto_aug)
        path_env = "{0}/program/perl-5.24.0/bin:{0}/bioinfo/Genomic/Sofware/augustus/bin".format(self.config.SOFTWARE_DIR)
        perl_lib = "{0}/program/perl-5.24.0/lib/site_perl/5.24.0:{0}/bioinfo/Genomic/Sofware/maker-2.31.9/maker/lib:{0}/bioinfo/Genomic/Sofware/maker-2.31.9/maker/perl/lib:{0}/app/program/perl-5.24.0/lib/site_perl/5.24.0/x86_64-linux-thread-multi".format(self.config.SOFTWARE_DIR)
        self.set_environ(PATH = path_env)
        self.set_environ(PERL5LIB = perl_lib)
        self.set_environ(AUGUSTUS_CONFIG_PATH=self.config.SOFTWARE_DIR+"/bioinfo/Genomic/Sofware/augustus/config")
        self.set_environ(LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + '/program/miniconda2/pkgs/boost-cpp-1.64.0-1/lib')  #add for ningbo zouguanqing 20181029
        self.augustus_path = "bioinfo/Genomic/Sofware/augustus/scripts/autoAug.pl"
        self.has_value_common_options = ["trainingset","hints","estali","cdna"]
        self.no_value_common_options = ["pasa","pasapolyAhints"]
        self.run_options = ""
        option_str = ""
        for ops in self.has_value_common_options:
            if self.option(ops):
                option_str += " --{} {} ".format(ops,self.option(ops))
        for ops in self.no_value_common_options:
            if self.option(ops) != "0":
                option_str += " --{} ".format(ops)
        option_str += self.option("augustus_other_options")
        self.run_options = "--genome={} --species={} {}".format(self.genome_fasta, self.option("species"), option_str)

        cmd = '{} {}'.format(self.augustus_path, self.run_options)
        command = self.add_command("run_augustus", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("augustus运行完成")
            self.augustus_out_gff = os.path.join(self.work_dir,"autoAug/results/predictions/augustus.abinitio.gff")
            self.augustus_out_faa = os.path.join(self.work_dir,"autoAug/results/predictions/augustus.abinitio.aa")
            if os.path.exists(self.augustus_out_gff):
                self.logger.info("运行生成{}".format(self.augustus_out_gff))
            else:
                self.logger.info("运行没有生成{}".format(self.augustus_out_gff))
            if os.path.exists(self.augustus_out_faa):
                self.logger.info("运行生成{}".format(self.augustus_out_faa))
            else:
                self.logger.info("运行没有生成{}".format(self.augustus_out_faa))

        else:
            self.set_error("augustus运行出错!", code="33300101")


    def set_output(self):
        gff_new_path = os.path.join(self.output_dir,"augustus.abinitio.gff")
        faa_new_path = os.path.join(self.output_dir,"augustus.abinitio.aa")
        if os.path.exists(gff_new_path):
            os.remove(gff_new_path)
        if os.path.exists(faa_new_path):
            os.remove(faa_new_path)
        os.link(self.augustus_out_gff,gff_new_path)
        os.link(self.augustus_out_faa,faa_new_path)




    def run(self):
        super(AugustusTool, self).run()
        self.run_augustus()
        #self.set_output()
        self.end()


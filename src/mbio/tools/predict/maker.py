# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'
# version 1.0
# last_modify: 2018.06.04

import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import shutil


class MakerAgent(Agent):
    """
    maker 进行基因预测
    """

    def __init__(self, parent):
        super(MakerAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"}, # 组装拼接好的scaffold文件
            {"name": "pyu_hmm", "type": "string", "default": ""},
            {"name": "gmhmm_mod", "type": "string", "default": ""},
            {"name": "base", "type": "string","default": "maker1"},
            {"name": "ref_protein", "type": "infile", "format": "sequence.fasta"},
            {"name": "species", "type": "string","default":""},
            {"name": "get_fasta","type": "string","default":"T"}  # 最终想获得faa和ffn序列，选T ，否则选F
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome", code="33300601")
        if not os.path.exists(self.option("pyu_hmm")):
            raise OptionError("pyu.hmm文件路径不存在", code="33300602")
        if not os.path.exists(self.option("gmhmm_mod")):
            raise OptionError("gmhmm.mod文件路径不存在", code="33300603")
        if not self.option("ref_protein").is_set and not self.option("species"):
            raise OptionError("ref_protein和species参数，必须设置一个", code="33300604")
        if self.option("ref_protein").is_set and  self.option("species"):
            raise OptionError("ref_protein和species参数，必须设置一个", code="33300605")


    def set_resource(self):
        self._cpu = 16
        self._memory = '40G'

    def end(self):
        super(MakerAgent, self).end()


class MakerTool(Tool):
    def __init__(self, config):
        super(MakerTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.sample_name = (self.genome_fasta).split('/')[-1].rsplit('.')[0]
        self.maker_out_log = ""
        self.maker_out_gff = ""
        self.maker_out_faa = ""
        self.maker_out_ffn = ""

        perl_lib="{0}/program/perl-5.24.0/lib/site_perl/5.24.0:{0}/bioinfo/Genomic/Sofware/maker-2.31.9/maker/lib:{0}/bioinfo/Genomic/Sofware/maker-2.31.9/maker/perl/lib:{0}/program/perl-5.24.0/lib/site_perl/5.24.0/x86_64-linux-thread-multi".format(self.config.SOFTWARE_DIR)
        self.set_environ(PERL5LIB=perl_lib)
        self.set_environ(PATH=self.config.SOFTWARE_DIR+"/program/perl-5.24.0/bin")
        self.set_environ(AUGUSTUS_CONFIG_PATH=self.config.SOFTWARE_DIR+"/bioinfo/Genomic/Sofware/augustus/config")
        self.set_environ(ZOE=self.config.SOFTWARE_DIR+"/bioinfo/Genomic/Sofware/SNAP-master/SNAP-master/Zoe")

    def run_maker(self):
        maker_out = os.path.join(self.work_dir,"{0}.maker.output".format(self.option("base")))
        if os.path.exists(maker_out):
            shutil.rmtree(maker_out)

        os.system("cp ~/.gm_key ./")
        self.maker_path = "/bioinfo/Genomic/Sofware/maker/bin/maker"
        cp_ctl = "cp {}/bioinfo/Genomic/Sofware/maker/ctl/*.ctl ./".format(self.config.SOFTWARE_DIR)
        cp_pyu_hmm = "cp {} ./".format(self.option("pyu_hmm"))
        cp_gmhmm_mod = "cp {} ./".format(self.option("gmhmm_mod"))
        os.system(cp_ctl)
        os.system(cp_pyu_hmm)
        os.system(cp_gmhmm_mod)

        if self.option("ref_protein").is_set:
            self.ref_protein = self.option("ref_protein").prop['path']
        else:
            self.ref_protein = ""
        sed_ref = "sed -i 's#^protein=#protein={}#' maker_opts.ctl".format(self.ref_protein)
        sed_genome = "sed -i 's#^genome=#genome={}#' maker_opts.ctl".format(self.genome_fasta)
        sed_spe = "sed -i 's/^augustus_species=/augustus_species={}/' maker_opts.ctl".format(self.option("species"))
        sed_exe = "sed -i 's#/mnt/ilustre/users/sanger-dev/app#{}#' maker_exe.ctl".format(self.config.SOFTWARE_DIR)
        os.system(sed_genome)
        os.system(sed_ref)
        os.system(sed_spe)
        os.system(sed_exe)

        cmd = '{} -g {}  -R -base {} -TMP /tmp'.format(self.maker_path, self.genome_fasta, self.option("base"))
        command = self.add_command("run_maker", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("maker运行完成")
            self.maker_out_log = os.path.join(self.work_dir,"{0}.maker.output/{0}_master_datastore_index.log".format(self.option("base")))
            if os.path.exists(self.maker_out_log):
                self.logger.info("maker运行生成{}".format(self.maker_out_log))
            else:
                self.logger.info("maker运行没有生成{}".format(self.maker_out_log))
        else:
            self.set_error("maker运行出错!", code="33300601")

    def run_gff3_merge(self):
        if not os.path.exists(self.maker_out_log):
            self.set_error("缺失文件%s", variables=(self.maker_out_log), code="33300602")
        self.gff3_merge = "bioinfo/Genomic/Sofware/maker/bin/gff3_merge"
        cmd = "{} -d {}".format(self.gff3_merge,self.maker_out_log)
        command = self.add_command("run_gff3_merge",cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gff3_merge运行成功")
            self.maker_out_gff = os.path.join(self.work_dir,"%s.all.gff"%(self.option("base")))
            if os.path.exists(self.maker_out_gff):
                self.logger.info("gff3_merge 生成文件%s"%(self.maker_out_gff))
            else:
                self.logger.info("gff3_merge 没有生成文件%s"%(self.maker_out_gff))
        else:
            self.logger.info("gff3_merge运行失败")

    def run_fasta_merge(self):
        if not os.path.exists(self.maker_out_log):
            self.set_error("缺失文件%s", variables=(self.maker_out_log), code="33300603")
        self.fasta_merge = "/bioinfo/Genomic/Sofware/maker/bin/fasta_merge"
        cmd = "{} -d {}".format(self.fasta_merge,self.maker_out_log)
        command = self.add_command("run_fasta_merge",cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("fasta_merge运行完成")
            self.maker_out_faa = os.path.join(self.work_dir,"%s.all.maker.proteins.fasta"%(self.option("base")))
            self.maker_out_ffn = os.path.join(self.work_dir,"%s.all.maker.transcripts.fasta"%(self.option("base")))
            if os.path.exists(self.maker_out_faa):
                self.logger.info("maker运行生成{}".format(self.maker_out_faa))
            else:
                self.logger.info("maker运行没有生成{}".format(self.maker_out_faa))

            if os.path.exists(self.maker_out_ffn):
                self.logger.info("maker运行生成{}".format(self.maker_out_ffn))
            else:
                self.logger.info("maker运行没有生成{}".format(self.maker_out_ffn))
        else:
            self.logger.info("fasta_merge运行失败")

    def set_output(self):
        new_gff = os.path.join(self.output_dir,"{}.all.gff".format(self.option("base")))
        if os.path.exists(new_gff):
            os.remove(new_gff)
        os.link(self.maker_out_gff,new_gff)
        if self.option("get_fasta") == "T":
            new_faa = os.path.join(self.output_dir,"{}.proteins.fasta".format(self.option("base")))
            new_ffn = os.path.join(self.output_dir,"{}.transcripts.fasta".format(self.option("base")))
            for i in [new_faa,new_ffn]:
                if os.path.exists(i):
                    os.remove(i)
            os.link(self.maker_out_faa,new_faa)
            os.link(self.maker_out_ffn,new_ffn)



    def run(self):
        super(MakerTool, self).run()
        self.run_maker()
        self.run_gff3_merge()
        if self.option("get_fasta") == "T":
            self.run_fasta_merge()
        self.set_output()
        self.end()


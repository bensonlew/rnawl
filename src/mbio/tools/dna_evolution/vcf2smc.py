# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.09.08

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import shutil


class Vcf2smcAgent(Agent):
    """
    工具：根据vcf文件，计算每个分群的每个chr的smc
    """
    def __init__(self, parent):
        super(Vcf2smcAgent, self).__init__(parent)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf", "required": True},  # pop.recode.vcf或者比较分析的结果
            {"name": "group_list", "type": "infile", "format": "dna_evolution.group_table", "required": True},  # 分组文件,输出文件以分组文件名称.分割命名
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 9
        self._memory = "30G"

    def end(self):
        super(Vcf2smcAgent, self).end()


class Vcf2smcTool(Tool):
    def __init__(self, config):
        super(Vcf2smcTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/7.2.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/7.2.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/Python35/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/program/Python35/lib')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/gcc_7.2/gmp-6.1.0/lib')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/gcc_7.2/gsl23/lib')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/gcc_7.2/mpfr-3.1.4/lib')
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.bgzip_path = "bioinfo/dna_evolution/bgzip"
        self.tabix_path = "bioinfo/dna_evolution/tabix"
        self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.smc_path = self.config.SOFTWARE_DIR + "/program/Python35/bin/smc++"

    def get_chr_list(self):
        """
        得到chr_list，以逗号分隔
        """
        chr_list, sca_list = [], []
        with open(self.option("vcf_file").prop["path"], "r") as f:
            for line in f:
                if line.startswith("##contig=<ID="):
                    chr = line.split(",")[0].split("##contig=<ID=")[1]
                    if chr.startswith("chr") or chr.startswith("Chr"):
                        chr_list.append(chr)
                    else:
                        sca_list.append(chr)
                if line.startswith("#CHROM"):
                    break
        if chr_list:
            self.chr_list = chr_list
        else:
            self.chr_list = sca_list

    def run_bgzip(self):
        """
        bgzip
        """
        base_name = os.path.basename(self.option("vcf_file").prop["path"])
        vcf_path = os.path.join(self.work_dir, base_name)
        if os.path.exists(vcf_path + ".gz"):
            os.remove(vcf_path + ".gz")
        if os.path.exists(vcf_path):
            os.remove(vcf_path)
        os.link(self.option("vcf_file").prop["path"], vcf_path)
        cmd = "{} {}".format(self.bgzip_path, vcf_path)
        command = self.add_command("bgzip", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bgzip运行成功")
        else:
            self.set_error("bgzip运行失败")

    def run_tabix(self):
        """
        tabix
        """
        base_name = os.path.basename(self.option("vcf_file").prop["path"])
        cmd = "{} -p vcf {}".format(self.tabix_path, os.path.join(self.work_dir, base_name + ".gz"))
        command = self.add_command("tabix", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("tabix运行成功")
        else:
            self.set_error("tabix运行失败")

    def run_smc(self):
        """
        smc++
        /mnt/ilustre/users/sanger-dev/app/program/Python35/bin/smc++ vcf2smc
        /mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_evolution/file/smc/1.vcf.gz
        chr1.smc.gz chr1 Pop1:Native_1,R01
        """
        cmd_list = []
        group_info = self.option("group_list").prop["group_info"]
        for group in group_info.keys():
            pop_info = group + ":" + ",".join(group_info[group])
            pop_dir = os.path.join(self.output_dir, group)
            if os.path.exists(pop_dir):
                shutil.rmtree(pop_dir)
            os.mkdir(pop_dir)
            for chr in self.chr_list:
                chr_smc = os.path.join(pop_dir, chr + ".smc.gz")
                base_name = os.path.basename(self.option("vcf_file").prop["path"])
                cmd = "{} vcf2smc {}".format(self.smc_path, os.path.join(self.work_dir, base_name + ".gz"))
                cmd += " {} {} {}".format(chr_smc, chr, pop_info)
                cmd_list.append(cmd)
        cmd_file = os.path.join(self.work_dir, "smc_pop_chr.list")
        wrong_cmd = os.path.join(self.work_dir, "failed_smc_pop_chr.txt")
        with open(cmd_file, "w") as f:
            for cmd in cmd_list:
                f.write(cmd + "\n")
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.parafly, cmd_file, 8, wrong_cmd)
        command = self.add_command("smc_pop_chr", cmd_more).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format("smc_pop_chr"))
        else:
            self.set_error("{}运行失败".format("smc_pop_chr"))

    def run(self):
        super(Vcf2smcTool, self).run()
        self.get_chr_list()
        self.run_bgzip()
        self.run_tabix()
        self.run_smc()
        self.end()

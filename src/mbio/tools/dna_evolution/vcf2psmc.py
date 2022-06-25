# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.09.08

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import shutil


class Vcf2psmcAgent(Agent):
    """
    工具：根据vcf文件，得到每个组的fa文件
    """
    def __init__(self, parent):
        super(Vcf2psmcAgent, self).__init__(parent)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf", "required": True},  # pop.recode.vcf或者比较分析的结果
            {"name": "group_list", "type": "infile", "format": "dna_evolution.group_table", "required": True},  # 分组文件,输出文件以分组文件名称.分割命名
            {"name": "psmcfa_list", "type": "outfile", "format": "dna_evolution.group_table"},  # psmcfa.list
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self.logger.info(self.option("group_list").prop["group_info"])
        group_num = len(self.option("group_list").prop["group_info"].keys())
        self._cpu = group_num + 1
        self._memory = "15G"

    def end(self):
        super(Vcf2psmcAgent, self).end()


class Vcf2psmcTool(Tool):
    def __init__(self, config):
        super(Vcf2psmcTool, self).__init__(config)
        self.perl_path = "miniconda2/bin/perl"
        # self.vcf2fa_path = self.config.PACKAGE_DIR + "/dna_evolution/vcf2fa.pl"
        self.vcf2fa_path = self.config.PACKAGE_DIR + "/dna_evolution/vcf2psmc.pl"
        self.fq2psmcfa_path = self.config.SOFTWARE_DIR + "/bioinfo/dna_evolution/fq2psmcfa"
        self.fq2psmcfa_sh = self.config.PACKAGE_DIR + "/dna_evolution/fq2psmcfa.sh"
        self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"

    def run_vcf2fa(self):
        """
        vcf2fa.pl
        """
        self.fa_dir = os.path.join(self.work_dir, "fasta_dir")
        if os.path.exists(self.fa_dir):
            shutil.rmtree(self.fa_dir)
        os.mkdir(self.fa_dir)
        # cmd = "{} {} -g {}".format(self.perl_path, self.vcf2fa_path, self.option("group_list").prop["path"])
        # cmd += " -i {} -o {}".format(self.option("vcf_file").prop["path"], self.fa_dir)
        cmd = "{} {} -i {}".format(self.perl_path, self.vcf2fa_path, self.option("vcf_file").prop["path"])
        cmd += " -p {} -o {}".format(self.option("group_list").prop["path"], self.fa_dir)
        command = self.add_command("vcf2fa", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("vcf2fa运行成功")
        else:
            self.set_error("vcf2fa运行失败")
        # self.fa_list = os.path.join(self.fa_dir, "fa.list")
        # if not os.path.exists(self.fa_list):
        #     self.set_error("vcf2fa没有生成fa.list文件，请检查")

    def get_group_fasta(self):
        """
        将同组的fasta文件cat到一起
        """
        self.fa_info = {}
        with open(self.option("group_list").prop["path"], "r") as f:
            for line in f:
                item = line.strip().split("\t")
                if item[1] not in self.fa_info.keys():
                    self.fa_info[item[1]] = []
                self.fa_info[item[1]].append(item[0])
        with open(self.output_dir + "/psmcfa.list", "w") as w:
            for group_name in self.fa_info.keys():
                fa_list = []
                for sample_name in self.fa_info[group_name]:
                    fa_list.append(os.path.join(self.fa_dir, sample_name + ".fasta"))
                group_fa = os.path.join(self.output_dir, group_name + ".fasta")
                os.system("cat {} > {}".format(" ".join(fa_list), group_fa))
                w.write(group_name + "\t" + group_fa + "\n")
        self.option("psmcfa_list", self.output_dir + "/psmcfa.list")

    def run_fq2psmcfa(self):
        """
        fq2psmcfa
        """
        cmd_list = []
        self.fa_info = {}
        with open(self.fa_list, "r") as f, open(self.output_dir + "/psmcfa.list", "w") as w:
            for line in f:
                item = line.strip().split("\t")
                self.fa_info[item[0]] = item[1]
                w.write(item[0] + "\t" + self.output_dir + "/" + item[0] + ".psmcfa" + "\n")
        for base_name in self.fa_info.keys():
            fa = self.fa_info[base_name]
            cmd = "{} {} -q20 {} {}".format(self.fq2psmcfa_sh, self.fq2psmcfa_path, fa, os.path.join(self.output_dir, base_name + ".psmcfa"))
            cmd_list.append(cmd)
        cmd_file = os.path.join(self.work_dir, "fq2psmcfa_cmd.list")
        wrong_cmd = os.path.join(self.work_dir, "failed_fq2psmcfa_cmd.txt")
        with open(cmd_file, "w") as f:
            for cmd in cmd_list:
                f.write(cmd + "\n")
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.parafly, cmd_file, len(cmd_list), wrong_cmd)
        command = self.add_command("fq2fsmcfa_more", cmd_more).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format("fq2fsmcfa_more"))
        else:
            self.set_error("{}运行失败".format("fq2fsmcfa_more"))
            raise Exception("{}运行失败".format("fq2fsmcfa_more"))
        self.option("psmcfa_list", self.output_dir + "/psmcfa.list")

    def run(self):
        super(Vcf2psmcTool, self).run()
        self.run_vcf2fa()
        # self.run_fq2psmcfa()
        self.get_group_fasta()
        self.end()

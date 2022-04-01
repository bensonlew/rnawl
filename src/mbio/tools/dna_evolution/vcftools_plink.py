    # -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180823

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class VcftoolsPlinkAgent(Agent):
    """
    两步命令：
        第一步：vcftool --vcf ./pop.recode.vcf --plink --out ./pop
        第二步：plink --file ./pop --make-bed --out ./pop --allow-extra-chr --chr-set 35
    flag是否运行第二步plink,单倍体图谱不需要，群体结构需要第二步
    module:hapmap调取该tool
    """
    def __init__(self, parent):
        super(VcftoolsPlinkAgent, self).__init__(parent)
        options = [
            {"name": "recode_vcf_path", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "flag", "type": "bool", "default": True},  # 默认做第二步
            {"name": "make_bed", "type": "bool", "default": False},   # --make-bed
            {"name": "allow_extra_chr", "type": "bool", "default": False},   # --allow-extra-chr
            {"name": "chr_set", "type": "int", "default": 35},   # --chr-set 35
            {"name": "blocks", "type": "string"},   # --blocks no-pheno-req
            {"name": "r2", "type": "bool", "default": False},  # --r2
            {"name": "ld_window", "type": "int", "default": 99999},  # --ld-window
            {"name": "ld_window_kb", "type": "int", "default": 500},  # --ld-window-kb, 固定大小
            {"name": "ld_window_r2", "type": "float", "default": 0.8},  # --ld-window-r2, R2
            {"name": "ldsnp_list", "type": "infile", "format": "dna_evolution.region_snp"},  # pop1_pop2.1.pi_tajimaD_fst.select.snp
            {"name": "compare_ldsnp_list", "type": "infile", "format": "dna_evolution.region_snp"},  # pop1_pop2.1.pi_tajimaD_fst.select.snp
        ]
        self.add_option(options)
        self.step.add_steps('VcftoolsPlink')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.VcftoolsPlink.start()
        self.step.update()

    def step_end(self):
        self.step.VcftoolsPlink.finish()
        self.step.update()

    def check_options(self):
        if not self.option("recode_vcf_path").is_set:
            raise OptionError("请设置recode_vcf_path")
        if not isinstance(self.option('make_bed'), bool):
            raise OptionError("make_bed参数类型错误")
        if not isinstance(self.option('allow_extra_chr'), bool):
            raise OptionError("allow_extra_chr参数类型错误")
        if self.option("ld_window_r2") > 1 or self.option("ld_window_r2") < 0:
            raise OptionError("R2的大小范围是0-1，而不是:{}".format(self.option("ld_window_r2")))
        if self.option("compare_ldsnp_list").is_set:
            if not self.option("ldsnp_list").is_set:
                raise OptionError("设置了compare_ldsnp_list文件的时候请同时设置ldsnp_list文件")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 8
        size = os.path.getsize(self.option("recode_vcf_path").prop['path'])
        size = size / 1024 / 1024 / 1024
        if size < 50:
            self._memory = "40G"
        elif size < 80:
            self._memory = "70G"
        elif size < 100:
            self._memory = "90G"
        else:
            self._memory = "120G"

    def end(self):
        super(VcftoolsPlinkAgent, self).end()


class VcftoolsPlinkTool(Tool):
    def __init__(self, config):
        super(VcftoolsPlinkTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.vcftools_path = "bioinfo/dna_evolution/vcftools"
        self.plink_path = "bioinfo/dna_evolution/plink/plink"
        self.parafly = "program/parafly-r2013-01-21/src/ParaFly"
        self.plink = self.config.SOFTWARE_DIR + "/bioinfo/dna_evolution/plink/plink"

    def run_vcftools_plink(self):
        """
        这一步慢一些
        """
        cmd = "{} --vcf {}".format(self.vcftools_path, self.option('recode_vcf_path').prop['path'])
        cmd += " --plink --out {}".format(self.output_dir + "/pop")
        self.logger.info(cmd)
        self.logger.info("开始进行VcftoolsPlink")
        command = self.add_command("vcftoolsplink", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("VcftoolsPlink完成！")
        else:
            self.set_error("VcftoolsPlink出错！")

    def run_plink(self):
        """
        较快，不耗内存
        :return:
        """
        cmd = "{}".format(self.plink_path)
        cmd += " --file {} ".format(self.output_dir + "/pop")
        if self.option('make_bed') is True:
            cmd += " --make-bed"
        cmd += " --out {}".format(self.output_dir + '/pop')
        if self.option('allow_extra_chr') is True:
            cmd += " --allow-extra-chr"
        if self.option('chr_set'):
            cmd += " --chr-set {}".format(int(self.option('chr_set')))
        self.logger.info(cmd)
        self.logger.info("开始进行Plink")
        command = self.add_command("plink", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Plink完成！")
        else:
            self.set_error("Plink出错！")

    def run_plink_snp_single(self):
        """
        对snp的结果进行plink
        """
        cmd1 = "{} --file {} --ld-window {}".format(self.plink_path, self.output_dir + "/pop", self.option("ld_window"))
        cmd1 += " --ld-snp-list {} --ld-window-kb {}".format(self.option("ldsnp_list").prop["path"], self.option("ld_window_kb"))
        cmd1 += " --ld-window-r2 {} --chr-set {}".format(self.option("ld_window_r2"), self.option("chr_set"))
        if self.option("r2"):
            cmd1 += " --r2"
        if self.option("allow_extra_chr"):
            cmd1 += " --allow-extra-chr"
        cmd1 += " --out {}".format(os.path.join(self.output_dir, os.path.basename(self.option("ldsnp_list").prop["path"])))
        command = self.add_command("plink_ld_snp_single", cmd1).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format("plink_ld_snp_single"))
        else:
            self.set_error("{}运行失败".format("plink_ld_snp_single"))

    def run_plink_snp_double(self):
        """
        对snp的结果进行plink
        """
        cmd1 = "{} --file {} --ld-window {}".format(self.plink, self.output_dir + "/pop", self.option("ld_window"))
        cmd1 += " --ld-snp-list {} --ld-window-kb {}".format(self.option("ldsnp_list").prop["path"], self.option("ld_window_kb"))
        cmd1 += " --ld-window-r2 {} --chr-set {}".format(self.option("ld_window_r2"), self.option("chr_set"))
        if self.option("r2"):
            cmd1 += " --r2"
        if self.option("allow_extra_chr"):
            cmd1 += " --allow-extra-chr"
        cmd1 += " --out {}".format(os.path.join(self.output_dir, os.path.basename(self.option("ldsnp_list").prop["path"])))
        cmd2 = " {} --file {} --ld-window {}".format(self.plink, self.output_dir + "/pop", self.option("ld_window"))
        cmd2 += " --ld-snp-list {} --ld-window-kb {}".format(self.option("compare_ldsnp_list").prop["path"], self.option("ld_window_kb"))
        cmd2 += " --ld-window-r2 {} --chr-set {}".format(self.option("ld_window_r2"), self.option("chr_set"))
        if self.option("r2"):
            cmd2 += " --r2"
        if self.option("allow_extra_chr"):
            cmd2 += " --allow-extra-chr"
        cmd2 += " --out {}".format(os.path.join(self.output_dir, os.path.basename(self.option("compare_ldsnp_list").prop["path"])))
        cmd_file = os.path.join(self.work_dir, "plink_ld_snp_cmd.list")
        wrong_cmd = os.path.join(self.work_dir, "plink_ld_snp_cmd.txt")
        with open(cmd_file, "w") as f:
            f.write(cmd1 + "\n")
            f.write(cmd2 + "\n")
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.parafly, cmd_file, 2, wrong_cmd)
        command = self.add_command("plink_ld_snp_double", cmd_more).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format("plink_ld_snp_double"))
        else:
            self.set_error("{}运行失败".format("plink_ld_snp_double"))

    def run(self):
        super(VcftoolsPlinkTool, self).run()
        self.run_vcftools_plink()
        if self.option("compare_ldsnp_list").is_set:
            self.run_plink_snp_double()
        elif self.option("ldsnp_list").is_set:
            self.run_plink_snp_single()
        elif self.option('flag') is True:
            self.run_plink()
        self.end()

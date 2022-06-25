# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.09.18

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class RegionSnpAgent(Agent):
    """
    工具： region_snp.pl
    通过region_snp.pl筛选出snp信息
    """
    def __init__(self, parent):
        super(RegionSnpAgent, self).__init__(parent)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf", "required": True},  # pop.recode.vcf或者比较分析的结果
            {"name": "select_file", "type": "infile", "format": "dna_evolution.region_snp", "required": True},  # pop1_pop2.1.pi_tajimaD_fst.select
            {"name": "compare_select_file", "type": "infile", "format": "dna_evolution.pi_tajimad_fst"},  # pop1_pop2.2.pi_tajimaD_fst.select
            {"name": "region_snp", "type": "outfile", "format": "dna_evolution.region_snp"},  # pop1_pop2.1.pi_tajimaD_fst.select.snp
            {"name": "compare_region_snp", "type": "outfile", "format": "dna_evolution.pi_tajimad_fst"},  # pop1_pop2.1.pi_tajimaD_fst.select.snp
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 3
        self._memory = "15G"

    def end(self):
        super(RegionSnpAgent, self).end()


class RegionSnpTool(Tool):
    def __init__(self, config):
        super(RegionSnpTool, self).__init__(config)
        self.perl = "miniconda2/bin/perl"
        self.perl_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/perl"
        self.parafly = "program/parafly-r2013-01-21/src/ParaFly"
        self.region_snp_path = self.config.PACKAGE_DIR + "/dna_evolution/region_snp.pl"

    def run_region_snp_single(self):
        """
        region_snp.pl
        """
        region_snp1 = self.output_dir + "/" + os.path.basename(self.option("select_file").prop["path"]) + ".snp"
        cmd = "{} {} -select {}".format(self.perl, self.region_snp_path, self.option("select_file").prop["path"])
        cmd += " -input {} -output {}".format(self.option("vcf_file").prop["path"], region_snp1)
        command = self.add_command("region_snp_single", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("region_snp_single运行成功")
        else:
            self.set_error("region_snp_single运行失败")
        self.option("region_snp", region_snp1)

    def run_region_snp_double(self):
        """
        region_snp.pl
        """
        region_snp1 = self.output_dir + "/" + os.path.basename(self.option("select_file").prop["path"]) + ".snp"
        region_snp2 = self.output_dir + "/" + os.path.basename(self.option("compare_select_file").prop["path"]) + ".snp"
        cmd1 = "{} {} -select {}".format(self.perl_path, self.region_snp_path, self.option("select_file").prop["path"])
        cmd1 += " -input {} -output {}".format(self.option("vcf_file").prop["path"], region_snp1)
        cmd2 = "{} {} -select {}".format(self.perl_path, self.region_snp_path, self.option("compare_select_file").prop["path"])
        cmd2 += " -input {} -output {}".format(self.option("vcf_file").prop["path"], region_snp2)
        cmd_file = os.path.join(self.work_dir, "region_snp_cmd.list")
        wrong_cmd = os.path.join(self.work_dir, "region_snp_cmd.txt")
        with open(cmd_file, "w") as f:
            f.write(cmd1 + "\n")
            f.write(cmd2 + "\n")
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.parafly, cmd_file, 2, wrong_cmd)
        command = self.add_command("region_snp_double", cmd_more).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format("region_snp_double"))
        else:
            self.set_error("{}运行失败".format("region_snp_double"))
        self.option("region_snp", region_snp1)
        self.option("compare_region_snp", region_snp2)

    def run(self):
        super(RegionSnpTool, self).run()
        if self.option("compare_select_file").is_set:
            self.run_region_snp_double()
        else:
            self.run_region_snp_single()
        self.end()

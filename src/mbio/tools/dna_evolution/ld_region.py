# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.08.22

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class LdRegionAgent(Agent):
    """
    工具： ld_region.pl
    通过ld_region.pl得到关联区域
    """
    def __init__(self, parent):
        super(LdRegionAgent, self).__init__(parent)
        options = [
            {"name": "snp_ld", "type": "infile", "format": "dna_evolution.test", "required": True},  # pop1_pop2.1.pi_tajimaD_fst.select.snp.ld
            {"name": "compare_snp_ld", "type": "infile", "format": "dna_evolution.test"},  # pop1_pop2.2.pi_tajimaD_fst.select.snp.ld
            # {"name": "region_file", "type": "outfile", "format": "dna_evolution.test"},  # pop1_pop2.1.pi_tajimaD_fst.select.region
            # {"name": "compare_region_file", "type": "outfile", "format": "dna_evolution.test"},  # pop1_pop2.1.pi_tajimaD_fst.select.region
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 3
        self._memory = "15G"

    def end(self):
        super(LdRegionAgent, self).end()


class LdRegionTool(Tool):
    def __init__(self, config):
        super(LdRegionTool, self).__init__(config)
        self.perl = "miniconda2/bin/perl"
        self.perl_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/perl"
        self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.ld_region = self.config.PACKAGE_DIR + "/dna_evolution/sweep_ld_region.pl"

    def run_ld_region_single(self):
        """
        ld_region.pl
        """
        ld_region1 = self.output_dir + "/" + os.path.basename(self.option("snp_ld").prop["path"]) + ".snp"
        cmd = "{} {} -input {} -output {}".format(self.perl, self.ld_region, self.option("snp_ld").prop["path"], ld_region1)
        command = self.add_command("ld_region_single", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ld_region_single运行成功")
        else:
            self.set_error("ld_region_single运行失败")
        # self.option("region_file", ld_region1)

    def run_ld_region_double(self):
        """
        ld_region.pl
        """
        ld_region1 = self.output_dir + "/" + os.path.basename(self.option("snp_ld").prop["path"]).split(".snp.ld")[0] + ".region "
        ld_region2 = self.output_dir + "/" + os.path.basename(self.option("compare_snp_ld").prop["path"]).split(".snp.ld")[0] + ".region "
        cmd1 = "{} {} -input {} -output {}".format(self.perl_path, self.ld_region, self.option("snp_ld").prop["path"], ld_region1)
        cmd2 = "{} {} -input {} -output {}".format(self.perl_path, self.ld_region, self.option("compare_snp_ld").prop["path"], ld_region2)
        cmd_file = os.path.join(self.work_dir, "ld_region_cmd.list")
        wrong_cmd = os.path.join(self.work_dir, "ld_region_cmd.txt")
        with open(cmd_file, "w") as f:
            f.write(cmd1 + "\n")
            f.write(cmd2 + "\n")
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.parafly, cmd_file, 2, wrong_cmd)
        command = self.add_command("ld_region_double", cmd_more).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format("ld_region_double"))
        else:
            self.set_error("{}运行失败".format("ld_region_double"))
        # self.option("region_file", ld_region1)
        # self.option("compare_region_file", ld_region2)

    def run(self):
        super(LdRegionTool, self).run()
        if self.option("compare_snp_ld").is_set:
            self.run_ld_region_double()
        else:
            self.run_ld_region_single()
        self.end()

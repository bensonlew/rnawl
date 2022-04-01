# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class RegionVcfAgent(Agent):
    """
    关联区域变异位点分析，需要region-vcf.pl
    perl ~/sg-users/xuanhongdong/BSA/script/region-vcf.pl -i ~/sg-users/xuanhongdong/BSA/Demo/data/pop.final.vcf
    -o ~/sg-users/xuanhongdong/BSA/Demo/region_vcf/region.threshold.vcf
    -r ~/sg-users/xuanhongdong/BSA/Demo/slidingwin/sliding-win.threshold.select
    version 1.0
    author: HONGDONG
    last_modified:2018.02.22
    laste modified by zengjing 增加参数wp、mp、wb、mb、step
    """

    def __init__(self, parent):
        super(RegionVcfAgent, self).__init__(parent)
        options = [
            {"name": "p_f_vcf", "type": "string"},  # pop.final.vcf
            {"name": "s_w_select", "type": "string"},  # sliding-win.threshold.select
            {"name": "wp", "type": "string"},  # 野生型亲本名称
            {"name": "mp", "type": "string"},  # 突变型亲本名称
            {"name": "wb", "type": "string"},  # 野生型混池名称
            {"name": "mb", "type": "string"}  # 突变型混池名称
            # {"name": "step", "type": "int"}  # 滑窗策略较小的数值
        ]
        self.add_option(options)
        self.step.add_steps('RegionVcf')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.RegionVcf.start()
        self.step.update()

    def step_end(self):
        self.step.RegionVcf.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('p_f_vcf'):
            raise OptionError('必须提供index_calc_result结果表', code="31500901")
        if not self.option('s_w_select'):
            raise OptionError('必须提供sliding-win.threshold.select结果表', code="31500902")
        if not self.option("mb"):
            raise OptionError('必须提供突变型混池mb', code="31500903")
        # if not self.option("step"):
        #     raise OptionError('必须提供滑窗策略较小值')
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(RegionVcfAgent, self).end()


class RegionVcfTool(Tool):

    def __init__(self, config):
        super(RegionVcfTool, self).__init__(config)
        self._version = "1.0.1"
        self.script_path = self.config.PACKAGE_DIR + '/bsa/region-vcf.pl'
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '

    def run_regionvcf(self):
        one_cmd = self.perl_path + self.script_path
        one_cmd += " -vcf {} -o {} -r {}".format(self.option("p_f_vcf"), self.output_dir + '/region.threshold.vcf',
                                                 self.option("s_w_select"))
        one_cmd += " -mb {}".format(self.option("mb"))
        if self.option("mp"):
            one_cmd += " -mp {}".format(self.option("mp"))
        if self.option("wb"):
            one_cmd += " -wb {}".format(self.option("wb"))
        if self.option("wp"):
            one_cmd += " -wp {}".format(self.option("wp"))
        self.logger.info(one_cmd)
        self.logger.info("开始运行one_cmd")
        cmd = self.add_command("one_cmd", one_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行one_cmd成功" )
        else:
            self.set_error("运行one_cmd出错！", code="31500901")
            self.set_error("运行one_cmd出错！", code="31500904")

    def run(self):
        """
        运行
        """
        super(RegionVcfTool, self).run()
        self.run_regionvcf()
        self.end()

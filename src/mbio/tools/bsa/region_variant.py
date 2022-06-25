# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class RegionVariantAgent(Agent):
    """
    关联区域变异位点分析，需要region-variant.pl
    version 1.0
    author: HONGDONG
    last_modified:2018.02.22
    """
    
    def __init__(self, parent):
        super(RegionVariantAgent, self).__init__(parent)
        options = [
            {"name": "i_c_result", "type": "string"},  # index_calc_result
            {"name": "s_w_select", "type": "string"}  # sliding-win.threshold.select
        ]
        self.add_option(options)
        self.step.add_steps('RegionVariant')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.RegionVariant.start()
        self.step.update()
        
    def step_end(self):
        self.step.RegionVariant.finish()
        self.step.update()
        
    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('i_c_result'):
            raise OptionError('必须提供index_calc_result结果表', code="31500801")
        if not self.option('s_w_select'):
            raise OptionError('必须提供sliding-win.threshold.select结果表', code="31500802")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '10G'
        
    def end(self):
        super(RegionVariantAgent, self).end()


class RegionVariantTool(Tool):

    def __init__(self, config):
        super(RegionVariantTool, self).__init__(config)
        self._version = "1.0.1"
        self.script_path = self.config.PACKAGE_DIR + '/bsa/region-variant.pl'
        self.perl_path = 'miniconda2/bin/perl '

    def run_regionvariant(self):
        one_cmd = self.perl_path + self.script_path + \
                  " -i {} -o {} -r {}".format(self.option("i_c_result"), self.output_dir + '/region.threshold.variant',
                                              self.option("s_w_select"))
        self.logger.info(one_cmd)
        self.logger.info("开始运行one_cmd")
        cmd = self.add_command("one_cmd", one_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行one_cmd成功")
        else:
            self.set_error("运行one_cmd出错！", code="31500801")
            self.set_error("运行one_cmd出错！", code="31500804")

    def run(self):
        """
        运行
        """
        super(RegionVariantTool, self).run()
        self.run_regionvariant()
        self.end()



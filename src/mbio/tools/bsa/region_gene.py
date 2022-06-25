# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class RegionGeneAgent(Agent):
    """
    关联区域变异位点分析，需要region-gene.pl
    perl ~/sg-users/xuanhongdong/BSA/script/region-gene.pl -a ~/sg-users/xuanhongdong/BSA/Demo/data/pop.summary
    -o ~/sg-users/xuanhongdong/BSA/Demo/region_gene/region.threshold.gene
    -i ~/sg-users/xuanhongdong/BSA/Demo/slidingwin/sliding-win.threshold.select
    version 1.0
    author: HONGDONG
    last_modified:2018.02.22
    """
    
    def __init__(self, parent):
        super(RegionGeneAgent, self).__init__(parent)
        options = [
            {"name": "pop_summary", "type": "string"},  # pop.summary
            {"name": "s_w_select", "type": "string"}  # sliding-win.threshold.select
        ]
        self.add_option(options)
        self.step.add_steps('RegionGene')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.RegionGene.start()
        self.step.update()
        
    def step_end(self):
        self.step.RegionGene.finish()
        self.step.update()
        
    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('pop_summary'):
            raise OptionError('必须提供pop_summary结果表', code="31500701")
        if not self.option('s_w_select'):
            raise OptionError('必须提供sliding-win.threshold.select结果表', code="31500702")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '10G'
        
    def end(self):
        super(RegionGeneAgent, self).end()


class RegionGeneTool(Tool):

    def __init__(self, config):
        super(RegionGeneTool, self).__init__(config)
        self._version = "1.0.1"
        self.script_path = self.config.PACKAGE_DIR + '/bsa/region-gene.pl'
        self.perl_path = 'miniconda2/bin/perl '

    def run_regiongene(self):
        one_cmd = self.perl_path + self.script_path + \
                  " -a {} -o {} -i {}".format(self.option("pop_summary"), self.output_dir + '/region.threshold.gene',
                                              self.option("s_w_select"))
        self.logger.info(one_cmd)
        self.logger.info("开始运行one_cmd")
        cmd = self.add_command("one_cmd", one_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行one_cmd成功")
        else:
            self.set_error("运行one_cmd出错！", code="31500701")
            self.set_error("运行one_cmd出错！", code="31500704")

    def run(self):
        """
        运行
        """
        super(RegionGeneTool, self).run()
        self.run_regiongene()
        self.end()



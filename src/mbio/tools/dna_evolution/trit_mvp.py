# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180823

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class TritMvpAgent(Agent):
    """
    GWAS中MVP的R程序
    cmd: MVP.single.R
    Tool: vcf2hapmap 的cmd:cmd: vcf2hapmap.pl输出pop_hapmap
    /mnt/ilustre/users/qingmei.cui/newmdt/sanger/7.pop/GWAS_20180730/step01.vcf-filter/pop.hapmap
    /mnt/ilustre/users/qingmei.cui/newmdt/sanger/7.pop/GWAS_20180730/step01.vcf-filter/trit//Chlorogenic_acid.trt
    """
    def __init__(self, parent):
        super(TritMvpAgent, self).__init__(parent)
        options = [
            {"name": "trait_name", "type": "string"},   # trt的list中第一列的name
            {"name": "trait_trt", "type": "string"},    # trt的list中第二列的path文件路径
            {"name": "pop_hapmap", "type": "string"},   # Tool: vcf2hapmap的输出
        ]
        self.add_option(options)
        self.step.add_steps('TritMvp')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.TritMvp.start()
        self.step.update()

    def step_end(self):
        self.step.TritMvp.finish()
        self.step.update()

    def check_options(self):
        if not self.option("trait_name"):
            raise OptionError("请设置trait_name")
        if not self.option("trait_trt"):
            raise OptionError("请设置trait_trt")
        if not self.option("pop_hapmap"):
            raise OptionError("请设置pop_hapmap参数")

    def set_resource(self):
        """
        运行所需资源
        hl这是8cpu 100g,实际上小一些
        """
        self._cpu = 8
        self._memory = '20G'

    def end(self):
        super(TritMvpAgent, self).end()


class TritMvpTool(Tool):
    def __init__(self, config):
        super(TritMvpTool, self).__init__(config)
        self.rscript_path = "program/R-3.3.3/bin/Rscript"
        self.perl_path = 'miniconda2/bin/perl '     # test
        # self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.mvp_path = self.config.PACKAGE_DIR + '/dna_evolution/MVP.single.R'

    def TritMvp(self):
        """
        要重新写下！！！
        :return:
        """
        if not os.path.exists(self.output_dir + "/" + self.option("trait_name")):
            os.mkdir(self.output_dir + "/" + self.option("trait_name"))
        cmd = "{} {} --trait  {} --hapmap {} --output {}"\
            .format(self.rscript_path, self.mvp_path, self.option("trait_trt"), self.option("pop_hapmap"),
                    self.output_dir + "/" + self.option("trait_name"))
        # test
        # scri = "/mnt/ilustre/users/sanger-dev/app/database/dna_wgs_geneome/check_dirfile.pl"
        # d = "/mnt/ilustre/users/sanger-dev/app/database/dna_wgs_geneome/Vitis_vinifera/CRIBI_University_of_Padua/v2.1/2015.04.16"
        # cmd = "{} {} -d {} -o {}".format(self.perl_path, scri, d, self.output_dir + "/" + self.option("trait_name"))
        # test
        self.logger.info(cmd)
        self.logger.info("开始进行TritMvp")
        command = self.add_command("tritmvp", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("TritMvp完成！")
        else:
            self.set_error("TritMvp出错！")
            raise Exception("TritMvp出错！")

    def run(self):
        super(TritMvpTool, self).run()
        self.TritMvp()
        self.end()

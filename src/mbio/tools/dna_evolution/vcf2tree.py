# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180830

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class Vcf2treeAgent(Agent):
    """
    群体结构子模块，将vcf转为tree格式的文件
    """
    def __init__(self, parent):
        super(Vcf2treeAgent, self).__init__(parent)
        options = [
            {"name": "recode_vcf_path", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "phylip_path", "type": "outfile", "format": "dna_evolution.phylip"},
            {"name": "pop_fasta", "type": "outfile", "format": "sequence.fasta"}
        ]
        self.add_option(options)
        self.step.add_steps('script')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.script.start()
        self.step.update()

    def step_end(self):
        self.step.script.finish()
        self.step.update()

    def check_options(self):
        if not self.option("recode_vcf_path").is_set:
            raise OptionError("缺少%s, 请添加%s!", variables=("recode_vcf_path", "vcf_file"), code="11111")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '15G'

    def end(self):
        super(Vcf2treeAgent, self).end()


class Vcf2treeTool(Tool):
    def __init__(self, config):
        super(Vcf2treeTool, self).__init__(config)
        self.perl_path = 'miniconda2/bin/perl'
        self.script_path = self.config.PACKAGE_DIR + "/dna_evolution/vcf2tree.pl"

    def vcf2tree(self):
        """
        perl /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/dna_evolution/vcf2tree.pl
         -i /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/evolution/jinhua/pop.recode.vcf -o ./pop
        :return:
        """
        cmd = "{} {} -i {} -o {}"\
            .format(self.perl_path, self.script_path, self.option("recode_vcf_path").prop['path'],
                    self.output_dir + "/pop")
        self.logger.info(cmd)
        self.logger.info("开始进行script")
        command = self.add_command("vcf2tree", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("script完成！")
            if os.path.exists(self.output_dir + "/pop.phylip"):
                self.option("phylip_path", self.output_dir + "/pop.phylip")
            if os.path.exists(self.output_dir + "/pop.fasta"):
                self.option("pop_fasta", self.output_dir + "/pop.fasta")
        else:
            self.set_error("script出错！")

    def run(self):
        super(Vcf2treeTool, self).run()
        self.vcf2tree()
        self.end()

# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'


from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class CalPpmAgent(Agent):
    """
    宏基因组ppm丰度计算
    """

    def __init__(self, parent):
        super(CalPpmAgent, self).__init__(parent)
        options = [
            {"name": "clean_stat", "type": "infile", "format": "meta.profile", "required":True},
            {"name": "geneset_table", "type": "infile", "format": "sequence.profile_table", "required":True}, ## reads_number.xls
            {"name": "out_table", "type": "outfile", "format": "meta.otu.otu_table,sequence.profile_table"},
        ]
        self.add_option(options)

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(CalPpmAgent, self).end()


class CalPpmTool(Tool):
    def __init__(self, config):
        super(CalPpmTool, self).__init__(config)
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        self.ppm_script = self.config.PACKAGE_DIR + '/metagenomic/scripts/ppm.pl'

    def run(self):
        """
        运行
        :return:
        """
        super(CalPpmTool, self).run()
        self.get_ppm_table()
        self.set_output()
        self.end()

    def get_ppm_table(self):
        reads_table = self.option("geneset_table").prop['path']
        clean_stat = self.option("clean_stat").prop['path']
        self.ppm = os.path.join(self.output_dir, "PPM.xls")
        cmd = '{} {} -i {} -stat {} -o {}'.format(self.perl_path , self.ppm_script, reads_table, clean_stat, self.ppm)
        command = self.add_command('calculate_ppm', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("calculate_ppm succeed")
        elif command.return_code in [-9]:  # add memory limit by guhaidong @ 20180417
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("calculate_ppm failed")

    def set_output(self):
        self.logger.info('开始设置输出结果文件')
        self.option("out_table", self.ppm)


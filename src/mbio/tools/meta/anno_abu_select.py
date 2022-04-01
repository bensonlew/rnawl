# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
# last modify by shaohua.yuan
# last modify date: 20171012

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class AnnoAbuSelectAgent(Agent):
    """
    宏基因组注释丰度筛选
    """

    def __init__(self, parent):
        super(AnnoAbuSelectAgent, self).__init__(parent)
        options = [
            {"name": "lowest_level_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "abu_num", "type": "string"},
            {"name": "abu_proportion", "type": "string"},
            {"name": "database", "type": "string", "default": "others"},
            {"name": "abu_select_level", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("lowest_level_profile").is_set:
            raise OptionError("必须设置最低层级注释丰度表")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        self.option("abu_select_level", self.output_dir + "/abu_select_gene.xls")
        super(AnnoAbuSelectAgent, self).end()


class AnnoAbuSelectTool(Tool):
    def __init__(self, config):
        super(AnnoAbuSelectTool, self).__init__(config)
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        self.script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/abu_select.pl'

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoAbuSelectTool, self).run()
        self.run_select_abu()
        #self.set_output()
        self.end()

    def select_table(self):
        lowest_level_profile = self.option("lowest_level_profile").prop['path']
        self.lowest_level_profile = os.path.join(self.output_dir, "select_lowest_gene.xls")
        cmd = '{} {} -i {} -sam {} -o {}'.format(self.python_path , self.python_script, lowest_level_profile, sams, self.lowest_level_profile)
        command = self.add_command('select_lowest_table', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("select_lowest_table succeed")
        else:
            self.set_error("select_lowest_table failed")
            raise Exception("select_lowest_table failed")

    def run_select_abu(self):
        self.logger.info("start annotation abudance select")
        self.logger.info(self.option("lowest_level_profile"))
        self.logger.info(self.option("database"))
        self.logger.info(self.option("abu_num"))
        self.lowest_level_profile = self.option("lowest_level_profile").prop['path']
        abu_num = self.option("abu_num")
        abu_proportion = self.option("abu_proportion")
        outfile = os.path.join(self.output_dir, "abu_select_gene.xls")
        cmd = self.perl_path + ' {} -i {} -num {} -pro {} -o {}'. \
            format(self.script, self.lowest_level_profile, abu_num,
                   abu_proportion, outfile)
        if self.option("database") == "nr":
            cmd += ' -database nr'
        command = self.add_command('abudance select', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("annotation abudance select succeed")
        else:
            self.set_error("annotation abudance select failed")
            raise Exception("annotation abudance select failed")


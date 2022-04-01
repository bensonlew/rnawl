# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __last_modify__ = '20191226'


from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
from mbio.packages.toolapps.common import link_dir

class FqToFaAgent(Agent):
    """
    将fq处理成fa格式
    """

    def __init__(self, parent):
        super(FqToFaAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq"},
            {"name": "fasta", "type": "outfile", "format": "sequence.fasta"},
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("fastq").is_set:
            raise OptionError("请输入read1文件！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '20G'


class FqToFaTool(Tool):
    def __init__(self, config):
        super(FqToFaTool, self).__init__(config)
        self.perl = "/program/perl/perls/perl-5.24.0/bin/perl"
        self.fatofa = self.config.PACKAGE_DIR + "/metagenomic/fq_to_fa.pl"
        self.file = self.option("fastq").prop['path'].split("/")[-1]
        self.out = self.work_dir + "/" + self.file + ".fa"

    def run_fqtofa(self):
        """
        description
        :return:
        """
        cmd = "{} {} -i {} -o {}".format(self.perl, self.fatofa, self.option("fastq").prop['path'], self.out)
        command = self.add_command("run_fqtofa", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_fqtofa运行完成！")
        else:
            self.set_error("run_fqtofa运行完成运行出错!")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        if os.path.exists(self.output_dir + "/" + self.file + ".fa"):
            os.remove(self.output_dir + "/" + self.file + ".fa")
        os.link(self.out, self.output_dir + "/" + self.file + ".fa")
        self.option("fasta", self.output_dir + "/" + self.file + ".fa")

    def run(self):
        super(FqToFaTool, self).run()
        self.run_fqtofa()
        self.set_output()
        self.end()
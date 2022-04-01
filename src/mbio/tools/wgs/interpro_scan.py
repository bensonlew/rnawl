# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20190225

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class InterproScanAgent(Agent):
    """
    参考基因组的pfam数据库注释
    """
    def __init__(self, parent):
        super(InterproScanAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入拆开后的fasta
            {"name": "cpu", "type": "int", 'default': 8}
        ]
        self.add_option(options)
        self.step.add_steps('analysis')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.analysis.start()
        self.step.update()

    def step_end(self):
        self.step.analysis.finish()
        self.step.update()

    def check_options(self):
        if not self.option("query"):
            raise OptionError("please input fasta file!")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = self.option('cpu') + 1
        self._memory = '30G'

    def end(self):
        super(InterproScanAgent, self).end()


class InterproScanTool(Tool):
    def __init__(self, config):
        super(InterproScanTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/7.2.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/7.2.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/bin")
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/wgs_v2")
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/program/Python35/bin")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/program/Python35/lib")
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/bin")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/library_gcc_7.2/lib")
        self.interpro_scan_path = "database/dna_wgs_geneome/Anno_database/InterProScan/" \
                                  "interproscan-5.30-69.0/interproscan.sh"

    def script_run(self):
        """
        interproscan-5.30-69.0/interproscan.sh -i sub.1.fa -o sub.1.fa.interPro.interproscan -f TSV -dp -goterms
         -pa -iprlookup -cpu 8  -t n
        :return:
        """
        outfile = os.path.join(self.output_dir,
                               os.path.basename(self.option('query').prop['path']) + '.interPro.interproscan')
        cmd = "{} -i {} -o {} -f TSV -dp -goterms -pa -iprlookup -cpu {} -t n"\
            .format(self.interpro_scan_path, self.option('query').prop['path'], outfile, self.option('cpu'))
        self.logger.info(cmd)
        self.logger.info("start interpro_scan")
        command = self.add_command("interpro_scan", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("interpro_scan script finish！")
        else:
            self.set_error("interpro_scan script failed！")

    def run(self):
        super(InterproScanTool, self).run()
        self.script_run()
        self.end()

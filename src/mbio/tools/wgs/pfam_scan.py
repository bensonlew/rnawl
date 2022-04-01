# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20190222

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class PfamScanAgent(Agent):
    """
    参考基因组的pfam数据库注释
    """
    def __init__(self, parent):
        super(PfamScanAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入拆开后的fasta
            {"name": "database", "type": "string"},  # pfam的db文件
            {"name": "cpu", "type": "int", 'default': 8},
            {"name": "anno_method", "type": "string", "default": "interpro_scan"}  # interpro_scan or pfam_scan
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
        if not self.option("database"):
            raise OptionError("please set pfam database!")
        if self.option("anno_method") not in ['pfam_scan', 'interpro_scan']:
            raise OptionError("anno method is not right, must be pfam_scan or interpro_scan")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = self.option('cpu') + 1
        self._memory = '30G'

    def end(self):
        super(PfamScanAgent, self).end()


class PfamScanTool(Tool):
    def __init__(self, config):
        super(PfamScanTool, self).__init__(config)
        if self.option("anno_method") == 'pfam_scan':
            self.set_environ(PERL5LIB=self.config.SOFTWARE_DIR + "/bioinfo/wgs_v2/PfamScan")
            self.set_environ(PATH=self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/bin")
            self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/wgs_v2")
            self.pfam_scan_path = "bioinfo/wgs_v2/PfamScan/pfam_scan.pl"
        else:
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
        pfam_scan.pl -fasta sub.1.fa -dir pfam/2018-10-9  -outfile sub.1.fa.pfam.anno -cpu 8 -as -translate all
        :return:
        """
        outfile = os.path.join(self.output_dir, os.path.basename(self.option('query').prop['path']) + '.pfam.anno')
        cmd = "{} -fasta {} -dir {} -outfile {} -cpu {} -as -translate all"\
            .format(self.pfam_scan_path, self.option('query').prop['path'], self.option("database"), outfile,
                    self.option('cpu'))
        self.logger.info(cmd)
        self.logger.info("start pfam_scan")
        command = self.add_command("pfam_scan", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("pfam_scan script finish！")
        else:
            if not self.check_error:
                self.set_error("pfam_scan出错！")
            else:
                self.logger.info("pfam_scan，正常退出！")

    def interpro_script_run(self):
        """
        interproscan-5.30-69.0/interproscan.sh -i sub.1.fa -o sub.1.fa.interPro.interproscan -f TSV -dp -goterms
         -pa -iprlookup -cpu 8  -t n

         print SH "$interproscan -i $_ -o $out/$file.Pfam.interproscan -appl Pfam -f TSV -dp -iprlookup -cpu 8";
        if ($type eq 'nuc'){
            print SH " -t n ";
        }
        :return:
        """
        outfile = os.path.join(self.output_dir,
                               os.path.basename(self.option('query').prop['path']) + '.Pfam.interproscan')
        cmd = "{} -i {} -o {} -appl Pfam -f TSV -dp -iprlookup -cpu {} -t n" \
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
        super(PfamScanTool, self).run()
        if self.option("anno_method") == 'pfam_scan':
            self.script_run()
        else:
            self.interpro_script_run()
        self.end()

    def check_error(self, file_path):
        with open(file_path, 'r') as r:
            for line in r:
                if re.match('Unrecognised character \[\*\] in.*', line):
                    return True
        return False

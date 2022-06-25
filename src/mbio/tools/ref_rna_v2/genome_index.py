# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'liubinxu'


class GenomeIndexAgent(Agent):
    """
    build genome index
    """
    def __init__(self, parent):
        super(GenomeIndexAgent, self).__init__(parent)
        options = [
            {'default': '', 'type': 'string', 'name': 'fasta'},
            {'default': '', 'type': 'string', 'name': 'gtf'},
            {'default': '', 'type': 'string', 'name': 'index'},
            {'default': '', 'type': 'string', 'name': 'transcript'},
            {'default': '', 'type': 'string', 'name': 'stat'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 13
        self._memory = "{}G".format('30')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        super(GenomeIndexAgent, self).end()


class GenomeIndexTool(Tool):
    """
    build genome index
    """
    def __init__(self, config):
        super(GenomeIndexTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.bowtie2_index = '/bioinfo/align/bowtie2-2.2.9/bowtie2-build'
        self.hisat2_index = '/bioinfo/align/hisat2/hisat2-2.1.0/hisat2-build'
        self.samtools_index = '/miniconda2/bin/samtools'
        self.gff_read = '/bioinfo/rna/cufflinks-2.2.1/gffread'
        self.stat = self.config.PACKAGE_DIR + "/ref_rna_v2/gtf2genome_stat.sh"
        self.seqkit = software_dir + '/bioinfo/seq/seqkit'
        self.tabadd = self.config.PACKAGE_DIR + "/ref_rna_v2/tabletools_add.pl"

        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run_build(self):
        cmd = '{} {} {}'.format(self.bowtie2_index, self.option('fasta'), self.option('index'))
        cmd_name = 'bowtie2_index'
        command1 = self.add_command(cmd_name, cmd)
        command1.run()

        cmd = '{} -p 10 {} {}'.format(self.hisat2_index, self.option('fasta'), self.option('index'))
        cmd_name = 'hisat2_index'
        command2 = self.add_command(cmd_name, cmd)
        command2.run()

        cmd = '{} faidx {}'.format(self.samtools_index, self.option('fasta'))
        cmd_name = 'samtools_index'
        command3 = self.add_command(cmd_name, cmd)
        command3.run()

        cmd = '{} {} -g {} -w {}'.format(self.gff_read, self.option("gtf"), self.option("fasta"), self.option("transcript"))
        cmd_name = 'gffread'
        command4 = self.add_command(cmd_name, cmd)
        command4.run()

        cmd = 'sh {} {} {} {} {} {}'.format(self.stat, self.option('fasta'), self.option('gtf'), self.option('stat'), self.seqkit, self.tabadd)
        cmd_name = 'stat'
        command5 = self.add_command(cmd_name, cmd)
        command5.software_dir = "/bin/"
        command5.run()

        self.wait()
        for command in [command1, command2, command3, command4]:
            if command.return_code == 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            elif command.return_code is None:
                self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
                command.rerun()
                self.wait()
                if command.return_code is 0:
                    self.logger.info("{} Finished successfully".format(cmd_name))
                else:
                    self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
            else:
                self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))

    def set_output(self):
        pass

    def run(self):
        super(GenomeIndexTool, self).run()
        self.run_build()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation_v2/'
        data = {
            "id": "genome_index" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v2.genome_index",
            "instant": False,
            "options": dict(
                fasta= test_dir + "Saccharomyces_cerevisiae.dna.toplevel.fa",
                gtf= test_dir + "Saccharomyces_cerevisiae.R64-1-1.39.gtf",
                index=test_dir + "Saccharomyces_cerevisiae.dna.toplevel_index",
                transcript=test_dir + "transcript.fa",
                stat=test_dir + "genome_stat.xls"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

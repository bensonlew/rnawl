# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'liubinxu'


class FastaUniqAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(FastaUniqAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'config', 'format': 'small_rna.ini'},
            {'default': 'unique.fasta', 'type': 'string', 'name': 'u'},
            {'default': 'uniq.fasta', 'type': 'string', 'name': 'uniq'},
            {'default': 'table.xls', 'type': 'string', 'name': 'table'},
            {'default': 'listfile', 'type': 'string', 'name': 'list'},
            {'default': True, 'type': 'bool', 'name': 'uniq_all'},
        ]
        self.add_option(options)
        self._memory_increase_step = 50

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('25')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
            ])
        """
        super(FastaUniqAgent, self).end()


class FastaUniqTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(FastaUniqTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.perl_path = 'program/perl-5.24.0/bin/perl'

        self.uniqueFastas = self.config.PACKAGE_DIR + "/small_rna/uniqueFastas.pl"

    def run_uniqueFastas(self):
        cmd = '{} {} '.format(self.perl_path, self.uniqueFastas)
        cmd += '-{} {} '.format("i", self.option("config").prop['path'])
        cmd += '-{} {} '.format("u", self.option("u"))
        cmd += '-{} {} '.format("uniq", self.option("uniq"))
        cmd += '-{} {} '.format("table", self.option("table"))
        cmd += '-{} {} '.format("list", self.option("list"))
        if self.option('uniq_all'):
            cmd += '-uniq_all'
        cmd_name = 'fasta_uniq'
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code in [-9, -7]:  # 加入return_code检测，在sanger超出内存的返回值为-9
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by shicaiping @ 20190128
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
        all_files = [self.option("u"), self.option("uniq"), self.option("table"), self.option("list")]
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(FastaUniqTool, self).run()
        self.run_uniqueFastas()
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
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA/data0'
        data = {
            "id": "FastaUniq" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "small_rna.fasta_uniq",
            "instant": False,
            "options": dict(
                config=test_dir + "/" + "qc_file.config",
                u="unique.fasta",
                uniq="uniq.fasta",
                table="table.xls",
                list="listfile",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

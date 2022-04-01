# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import ConfigParser
import unittest
__author__ = 'fengyitong'


class IsomiridAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(IsomiridAgent, self).__init__(parent)
        options = [
            {'name': 'config', 'type': 'string'},
            {'name': 'ref', 'type': 'string'},
            {'name': 'pre', 'type': 'string'},
            {'name': 'mature', 'type': 'string'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 6
        self._memory = "{}G".format('60')

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
        super(IsomiridAgent, self).end()

class IsomiridTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(IsomiridTool, self).__init__(config)
        self.python_path = 'program/Python/bin/python'
        self.isomiRID = self.config.PACKAGE_DIR + "/small_rna/isomiRID.py"
        self.bowtie = self.config.SOFTWARE_DIR + "/bioinfo/align/bowtie-1.1.2/bowtie"
        self.bowtie_build = self.config.SOFTWARE_DIR + "/bioinfo/align/bowtie-1.1.2/bowtie-build"

    def run_isomiRID(self):
        cmd = '{} {} {} '.format(self.python_path, self.isomiRID, self.creat_config())
        cmd_name = 'isomirs'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
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

    def creat_config(self):
        config_txt = ''
        cfg = ConfigParser.ConfigParser()
        cfg.read(self.option('config'))
        libs = cfg.options('FASTA')
        config_txt += 'bowtie_path: ' + self.bowtie + '\n'
        config_txt += 'bowtie-build_path: ' + self.bowtie_build + '\n'
        for n, lib in enumerate(libs):
            config_txt += 'lib: ' + cfg.get('FASTA', lib) + ' ' + 'Sample' + str(n+1) + ' ' + 'fa' + '\n'
        config_txt += 'main_ref: ' + self.option('pre') + ' ' + 'yes' + '\n'
        config_txt += 'filter_ref: ' + self.option('ref') + ' ' + 'yes' + '\n'
        config_txt += 'known_miRNAs: ' + self.option('mature') + ' ' + 'yes' + '\n'
        config_txt += '''
M3: yes 4
M5: yes 3
RangeSize: 18 26
cutoff: 20
cutoff: 50
cutoff: 100
        '''
        with open('Config.txt', 'w') as config_w:
            config_w.write(config_txt)
        return os.getcwd() + '/Config.txt'

    def set_output(self):
        all_files = os.listdir(self.work_dir)
        for each in all_files:
            if each.endswith('.xls'):
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)

    def run(self):
        super(IsomiridTool, self).run()
        self.run_isomiRID()
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
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/smallrna/family'
        data = {
            "id": "isomiRID" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "small_rna.isomiRID",
            "instant": False,
            "options": dict(
                config='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/smallrna/bias/Uniq.cfg.ini',
                pre='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/smallrna/isomiRID/mmu.hairpin.fa',
                mature='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/smallrna/isomiRID/mmu.mature.fa',
                ref='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/smallrna/isomiRID/ref.fa'
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
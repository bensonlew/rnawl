# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'

import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest

class Gzfastq2fastqAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(Gzfastq2fastqAgent, self).__init__(parent)
        options = [
            {'name': 'fastq_path', 'type': 'string', 'default': ''},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = '60G'

    def end(self):
        super(Gzfastq2fastqAgent, self).end()

class Gzfastq2fastqTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(Gzfastq2fastqTool, self).__init__(config)
        self.parafly = self.config.SOFTWARE_DIR + '/program/parafly-r2013-01-21/src/ParaFly'

    def gzfastq2fastq(self):
        list_info = ''
        with open(self.option("fastq_path") + '/list.txt', 'r') as list_r, open('ungz.bash', 'w') as ungz_w:
            list_info = list_r.read()
            list_info_new = ''
            for line in list_info.split('\n'):
                line = line.strip().split('\t')
                special_character = ['&', '>', '<']
                for character in special_character:
                    if character in line[0]:
                        self.set_error("{}中存在特殊字符{}，无法解压，请重命名后再分析！".format(line[0], character))
                if line[0].split(".")[-1] in ["gz"]:
                    gzed = ".".join(line[0].split(".")[:-2]) + ".fastq"
                    list_info_new += gzed + '\t' + '\t'.join(line[1:]) + '\n'
                    ungz_w.write('gunzip -c ' + self.option('fastq_path') + '/' + line[0] + ' > ' + self.option('fastq_path') + '/' + gzed +'\n')
                else:
                    if line:
                        list_info_new += '\t'.join(line) + '\n'
        os.remove(self.option("fastq_path") + '/list.txt')
        with open(self.option("fastq_path") + '/list.txt', 'w') as list_w:
            list_w.write(list_info_new.strip('\n') + '\n')
        cmd = self.parafly + ' -c ungz.bash -CPU 10'
        cmd_name = 'gzfastq2fastq'
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.software_dir = ""
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
                self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35003701")
        else:
            self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35003702")

    def run(self):
        super(Gzfastq2fastqTool, self).run()
        self.gzfastq2fastq()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/pipline/ref'
        data = {
            "id": "gzfastq2fastq_gff" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "prok_rna.gzfastq2fastq",
            "instant": False,
            "options": dict(
                fastq_path = '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/data1'
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()

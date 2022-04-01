# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool


class Pigzfastq2fastqAgent(Agent):
    def __init__(self, parent):
        super(Pigzfastq2fastqAgent, self).__init__(parent)
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
        super(Pigzfastq2fastqAgent, self).end()


class Pigzfastq2fastqTool(Tool):
    def __init__(self, config):
        super(Pigzfastq2fastqTool, self).__init__(config)
        self.parafly = self.config.SOFTWARE_DIR + '/program/parafly-r2013-01-21/src/ParaFly'
        self.pigz = "/bioinfo/seq/pigz-2.4/"
        self.set_environ(PATH=self.config.SOFTWARE_DIR + self.pigz)

    def gzfastq2fastq(self):
        with open(self.option("fastq_path") + '/list.txt', 'r') as list_r, open('ungz.bash', 'w') as ungz_w:
            list_info = list_r.read()
            list_info_new = ''
            for line in list_info.split('\n'):
                line = line.strip().split('\t')
                if line[0].split(".")[-1] in ["gz"]:
                    gzed = ".".join(line[0].split(".")[:-2]) + ".fastq"
                    list_info_new += gzed + '\t' + '\t'.join(line[1:]) + '\n'
                    ungz_w.write('pigz -c -p 4 -k -d ' + self.option('fastq_path') + '/' + line[0] + ' > ' + self.output_dir+ '/' + gzed + '\n')
                else:
                    if line:
                        list_info_new += '\t'.join(line) + '\n'

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
                self.set_error("%s Failed. >>> %s", variables=(cmd_name, cmd), code="35003701")
        else:
            self.set_error("%s Failed. >>> %s", variables=(cmd_name, cmd), code="35003702")
        # os.remove(self.option("fastq_path") + '/list.txt')
        # with open(self.option("fastq_path") + '/list.txt', 'w') as list_w:
        #     list_w.write(list_info_new.strip('\n') + '\n')

    def gzfastq2fastq_new(self):
        with open(self.option("fastq_path") + '/list.txt', 'r') as list_r, open('ungz.bash', 'w') as ungz_w:
            list_info = list_r.read()
            list_info_new = ''
            list_info_new = '#!/bin/bash' + "\n"
            for line in list_info.split('\n'):
                line = line.strip().split('\t')
                if line[0].split(".")[-1] in ["gz"]:
                    gzed = ".".join(line[0].split(".")[:-2]) + ".fastq"
                    list_info_new += gzed + '\t' + '\t'.join(line[1:]) + '\n'
                    ungz_w.write('{}pigz -c -p 10 -k -d '.format(self.config.SOFTWARE_DIR + self.pigz) + self.option('fastq_path') + '/' + line[0] + ' > ' + self.output_dir+ '/' + gzed + '\n')
                else:
                    if line:
                        list_info_new += '\t'.join(line) + '\n'

        # cmd = self.parafly + ' -c ungz.bash -CPU 10'
        code = os.system('chmod +x {}'.format('ungz.bash'))
        cmd = '{}/ungz.bash '.format(self.work_dir)
        cmd_name = 'gzfastq2fastq'
        command = self.add_command(cmd_name, cmd, shell=True,ignore_error=True)
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
                self.set_error("%s Failed. >>> %s", variables=(cmd_name, cmd), code="35003701")
        else:
            self.set_error("%s Failed. >>> %s", variables=(cmd_name, cmd), code="35003702")

    def gzfastq2fastq_single(self):
        with open(self.option("fastq_path") + '/list.txt', 'r') as list_r, open('ungz.bash', 'w') as ungz_w:
            list_info = list_r.read()
            list_info_new = ''
            cmd_list =[]
            for line in list_info.split('\n'):
                line = line.strip().split('\t')
                if line[0].split(".")[-1] in ["gz"]:
                    gzed = ".".join(line[0].split(".")[:-2]) + ".fastq"
                    list_info_new += gzed + '\t' + '\t'.join(line[1:]) + '\n'
                    cmd = '{}/pigz -c -p 20 -k -d '.format(self.config.SOFTWARE_DIR + self.pigz) + self.option('fastq_path') + '/' + line[0] + ' > ' + self.output_dir+ '/' + gzed
                    cmd_list.append(cmd)
                    ungz_w.write('pigz -c -p 20 -k -d ' + self.option('fastq_path') + '/' + line[0] + ' > ' + self.output_dir+ '/' + gzed + '\n')
                else:
                    if line:
                        list_info_new += '\t'.join(line) + '\n'

        cmd = cmd
        cmd_name = 'gzfastq2fastq'
        command = self.add_command(cmd_name, cmd,shell=True, ignore_error=True)
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
                self.set_error("%s Failed. >>> %s", variables=(cmd_name, cmd), code="35003701")
        else:
            self.set_error("%s Failed. >>> %s", variables=(cmd_name, cmd), code="35003702")

    def run(self):
        super(Pigzfastq2fastqTool, self).run()
        self.gzfastq2fastq()
        # self.gzfastq2fastq_new()
        # self.gzfastq2fastq_single()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        data = {
            "id": "zhijie_gzfastq2fastq_single" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "ref_rna_v2.pigzfastq2fastq",
            "instant": False,
            "options": dict(
                fastq_path='/mnt/ilustre/users/isanger/workspace/20210301/Denovorna_majorbio_324951/remote_input/fastq_dir/Rawdata_FX2021022200180'
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

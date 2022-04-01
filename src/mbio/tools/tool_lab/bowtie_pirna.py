# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import unittest
from Bio import SeqIO


class BowtiePirnaAgent(Agent):
    '''
    last_modify: 2019.04.09
    '''
    def __init__(self, parent):
        super(BowtiePirnaAgent, self).__init__(parent)
        options = [
            {'name': 'fa', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},

        ]
        self.add_option(options)


    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(BowtiePirnaAgent, self).end()

class BowtiePirnaTool(Tool):
    def __init__(self, config):
        super(BowtiePirnaTool, self).__init__(config)
        self.cmd_path = "bioinfo/align/bowtie-1.2.3-linux-x86_64/bowtie-build"

    def run(self):
        super(BowtiePirnaTool, self).run()
        self.bowtie()
        self.set_output()
        self.end()

    def bowtie(self):
        cmd = '{} {} {}'.format(self.cmd_path, self.option('fa').path, self.option('fa').path)
        cmd_name = 'bowtie'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd), code="33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd), code="33704307")


    def len_split(self):
        fasta = self.option('fa').path
        name = os.path.basename(fasta).split('.')[0] + '_' + 'len'
        fasta_len_dir = os.path.join(os.path.dirname(fasta), name)
        os.mkdir(fasta_len_dir)
        pirna = SeqIO.parse(fasta, 'fasta')
        for record in pirna:
            length = len(str(record.seq))
            with open('{}/{}_{}.fa'.format(fasta_len_dir, name, length), 'a+') as f:
                f.write('>{}'.format(record.id) + '\n' + str(record.seq) + '\n')



    def run_makedb(self):
        cmd = os.path.join(self.cmd_path, "makeblastdb")
        db_name = os.path.join(self.option('db_name'), os.path.basename(self.option('db_name')))
        cmd += ' -dbtype {} -in {} -parse_seqids -out {}'.format('nucl', self.option('fa').path, db_name)
        cmd_name = 'run_makedb'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")

    def set_output(self):
        pass

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bowtie__{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'tool_lab.fasta_len_split',
            'instant': False,
            'options': {
                'fa': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/smallrna/piRNA/ssc/ssc.fa'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

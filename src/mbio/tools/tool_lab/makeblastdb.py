# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import unittest

class MakeblastdbAgent(Agent):
    '''
    last_modify: 2019.04.09
    '''
    def __init__(self, parent):
        super(MakeblastdbAgent, self).__init__(parent)
        options = [
            {'name': 'fa', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'db_name', 'type': 'string'}
        ]
        self.add_option(options)


    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(MakeblastdbAgent, self).end()

class MakeblastdbTool(Tool):
    def __init__(self, config):
        super(MakeblastdbTool, self).__init__(config)
        self.cmd_path = "bioinfo/align/ncbi-blast-2.3.0+/bin"

    def run(self):
        super(MakeblastdbTool, self).run()
        self.run_makedb()
        self.set_output()
        self.end()

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
    def test_pe(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'kallisto_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.expression.kallisto',
            'instant': False,
            'options': {
                'fq_type': 'PE',
                'fasta_file': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/all_transcripts.fa',
                'fastq_l_file': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/sickle_dir/Ctr_Liver_1_sickle_l.fastq',
                'fastq_r_file': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/sickle_dir/Ctr_Liver_1_sickle_r.fastq',
                'sample_name': 'Ctr_Liver_1'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_se(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'kallisto_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.expression.kallisto',
            'instant': False,
            'options': {
                'fq_type': 'SE',
                'fasta_file': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/all_transcripts.fa',
                'fastq_file': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/sickle_dir/Ctr_Liver_1_sickle_l.fastq',
                'sample_name': 'Ctr_Liver_1'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_se')])
    unittest.TextTestRunner(verbosity=2).run(suite)

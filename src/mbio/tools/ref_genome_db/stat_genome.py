# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
import os
from biocluster.tool import Tool
import unittest

class StatGenomeAgent(Agent):
    def __init__(self, parent):
        super(StatGenomeAgent, self).__init__(parent)
        options = [
            {'name': 'genome', 'type': 'infile', 'format': 'ref_genome_db.fasta'},
            # ensembl: gtf, ncbi: gff
            {'name': 'annotation', 'type': 'infile', 'format': 'ref_genome_db.gtf'},
            # unclass, ensembl, ncbi
            {'name': 'source', 'type': 'string', 'default': 'ensembl'},
            {'name': 'result_xls', 'type': 'outfile', 'format': 'ref_genome_db.common'}
        ]
        self.add_option(options)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} - {}'.format('genome', self.option('genome').prop['path']))
        self.logger.debug('{} - {}'.format('annotation', self.option('annotation').prop['path']))
        self.logger.debug('{} - {}'.format('source', self.option('source')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 2
        self._memory = '{}G'.format(os.path.getsize(self.option('genome').prop['path'])/1024**3+10)

    def end(self):
        super(StatGenomeAgent, self).end()

class StatGenomeTool(Tool):
    def __init__(self, config):
        super(StatGenomeTool, self).__init__(config)
        self.seqkit = 'bioinfo/seq/seqkit'
        self.python = 'miniconda2/bin/python'
        self.stat_annotation = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db/stat_annotation.py')
        self.stat_merge = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db/stat_merge.py')
        self.annotation = self.option('annotation').prop['path']
        self.tabular = os.path.join(self.work_dir, '{}.tabular'.format(os.path.basename(self.annotation)))
        self.dna = self.option('genome').prop['path']
        self.fx2tab = os.path.join(self.work_dir, '{}.fx2tab'.format(os.path.basename(self.dna)))
        self.result_xls = os.path.join(self.work_dir, '{}.genome_stat.xls'.format(os.path.basename(self.annotation)))

    def run(self):
        super(StatGenomeTool, self).run()
        self.run_stat_annotation()
        self.run_seqkit_fx2tab()
        self.run_stat_merge()
        self.set_output()
        self.end()

    def run_stat_annotation(self):
        cmd = '{} {} '.format(self.python, self.stat_annotation)
        cmd += '-i {} '.format(self.annotation)
        cmd += '-s {} '.format(self.option('source'))
        cmd += '-o {}'.format(self.tabular)
        cmd_name = 'stat_annotation'
        self.run_code(cmd_name, cmd)

    def run_seqkit_fx2tab(self):
        cmd = '{} fx2tab -g -l -n -i '.format(self.seqkit)
        cmd += '{} '.format(self.dna)
        cmd += '-o {}'.format(self.fx2tab)
        cmd_name = 'seqkit_fx2tab'
        self.run_code(cmd_name, cmd)

    def run_stat_merge(self):
        cmd = '{} {} '.format(self.python, self.stat_merge)
        cmd += '-t {} '.format(self.tabular)
        cmd += '-f {} '.format(self.fx2tab)
        cmd += '-o {}'.format(self.result_xls)
        cmd_name = 'stat_merge'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd):
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, try again'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('succeed in rerunning {}'.format(cmd_name))
            else:
                self.set_error('fail to rerun {}, abord'.format(cmd_name))
        else:
            self.set_error('fail to run {}, abord'.format(cmd_name))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        source = self.result_xls
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('result_xls', link_name)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_ensembl(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'stat_genome_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db.stat_genome',
            'instant': False,
            'options': {
                'genome': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_genome_db/ensembl/Acyrthosiphon_pisum.GCA_000142985.2.dna_rm.toplevel.fa',
                'annotation': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_genome_db/ensembl/Acyrthosiphon_pisum.GCA_000142985.2.36.gtf',
                'source': 'ensembl',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_ncbi(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'stat_genome_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db.stat_genome',
            'instant': False,
            'options': {
                'genome': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_genome_db/ncbi/GCF_001659605.1_Manihot_esculenta_v6_genomic.fna',
                'annotation': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_genome_db/ncbi/GCF_001659605.1_Manihot_esculenta_v6_genomic.gff',
                'source': 'ncbi',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_unclass_gtf(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'stat_genome_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db.stat_genome',
            'instant': False,
            'options': {
                'genome': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_genome_db/unclass/haemonchus_contortus.PRJEB506.WBPS11.genomic.fa',
                'annotation': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_genome_db/unclass/haemonchus_contortus.PRJEB506.WBPS11.gtf',
                'source': 'unclass',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_unclass_gff(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'stat_genome_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db.stat_genome',
            'instant': False,
            'options': {
                'genome': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_genome_db/ncbi/GCF_001659605.1_Manihot_esculenta_v6_genomic.fna',
                'annotation': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_genome_db/unclass/haemonchus_contortus.PRJEB506.WBPS11.annotations.gff3',
                'source': 'unclass',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
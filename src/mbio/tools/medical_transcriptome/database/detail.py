# -*- coding: utf-8 -*-
# __author__ = 'gudeqing,qinjincheng'

from biocluster.agent import Agent
from mbio.packages.ref_rna_v3.functions import toolfuncdeco, runcmd
from biocluster.tool import Tool
import os
import unittest

class DetailAgent(Agent):
    '''
    last_modify: 2019.06.24
    '''
    def __init__(self, parent):
        super(DetailAgent, self).__init__(parent)
        options = [
            {'name': 'gene_fa', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'txpt_fa', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'biomart_file', 'type': 'infile', 'format': 'ref_rna_v3.common'},
            {'name': 'biomart_type', 'type': 'string', 'format': None},
            {'name': 't2g', 'type': 'infile', 'format': 'ref_rna_v3.common'},
            {'name': 'known_cds', 'type': 'infile', 'format': 'ref_rna_v3.common'},
            {'name': 'known_pep', 'type': 'infile', 'format': 'ref_rna_v3.common'},
            {'name': 'novel_cds', 'type': 'infile', 'format': 'ref_rna_v3.common'},
            {'name': 'novel_pep', 'type': 'infile', 'format': 'ref_rna_v3.common'},
            {'name': 'database', 'type': 'outfile', 'format': 'ref_rna_v3.common'}
        ]
        self.add_option(options)
        self.step.add_steps('detail')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def step_start(self):
        self.step.detail.start()
        self.step.update()

    def step_end(self):
        self.step.detail.finish()
        self.step.update()

    def set_resource(self):
        self._cpu = 1
        size = os.path.getsize(self.option('gene_fa').path) + os.path.getsize(self.option('txpt_fa').path)
        self._memory = '{}G'.format(int(size / 1024.0 ** 3 * 10 + 10))

    @toolfuncdeco
    def end(self):
        super(DetailAgent, self).end()

class DetailTool(Tool):
    def __init__(self, config):
        super(DetailTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'detail': os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/database/detail.py')
        }
        self.file = {
            'database': os.path.join(self.output_dir, 'refrna_seqs.db'),
            'output_detail':os.path.join(self.output_dir, "detail")
        }

    def run(self):
        super(DetailTool, self).run()
        self.run_detail()
        self.set_output()
        self.end()

    @toolfuncdeco
    def run_detail(self):
        cmd = '{} {}'.format(self.program['python'], self.script['detail'])
        cmd += ' --gf {}'.format(self.option('gene_fa').path)
        cmd += ' --tf {}'.format(self.option('txpt_fa').path)
        cmd += ' --bf {}'.format(self.option('biomart_file').path)
        cmd += ' --bt {}'.format(self.option('biomart_type'))
        cmd += ' --t2g {}'.format(self.option('t2g').path)
        cmd += ' --kcf {}'.format(self.option('known_cds').path)
        cmd += ' --kpf {}'.format(self.option('known_pep').path)
        if self.option('novel_cds').is_set and self.option('novel_pep').is_set:
            cmd += ' --ncf {}'.format(self.option('novel_cds').path)
            cmd += ' --npf {}'.format(self.option('novel_pep').path)
        cmd += ' --output {}'.format(self.file['database'])
        cmd += ' --output_detail {}'.format(self.file['output_detail'])
        cmd_name = 'run_detail'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def set_output(self):
        self.option('database').set_path(self.file['database'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'medical_transcriptome.database.detail',
            'instant': False,
            'options': {
                'gene_fa': '/mnt/ilustre/users/sanger-dev/workspace/20201127/MedicalTranscriptome_6sieeiseu0ubhjnsm1o6302l5r/GeneFa/output/gene.fasta',
                'biomart_file':'/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/biomart/biomart.txt',
                'biomart_type':'type1',
                'known_cds':'/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/cds/Homo_sapiens.GRCh38.cds.all.fa',
                'known_pep':'/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/cds/Homo_sapiens.GRCh38.pep.all.fa',
                'txpt_fa':'/mnt/ilustre/users/sanger-dev/workspace/20201127/MedicalTranscriptome_6sieeiseu0ubhjnsm1o6302l5r/RefrnaAssemble/output/NewTranscripts/all_transcripts.fa',
                'novel_cds':'/mnt/ilustre/users/sanger-dev/workspace/20201127/MedicalTranscriptome_6sieeiseu0ubhjnsm1o6302l5r/AnnotOrfpfam/output/new_transcripts.fa.transdecoder.cds',
                'novel_pep':"/mnt/ilustre/users/sanger-dev/workspace/20201127/MedicalTranscriptome_6sieeiseu0ubhjnsm1o6302l5r/AnnotOrfpfam/output/new_transcripts.fa.transdecoder.pep",
                't2g':'/mnt/ilustre/users/sanger-dev/workspace/20201127/MedicalTranscriptome_6sieeiseu0ubhjnsm1o6302l5r/RefrnaAssemble/output/NewTranscripts/trans2gene',
                # 'jcpklist': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/count/jc.pk.list',
                # 'jcecpklist': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/count/jcec.pk.list'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

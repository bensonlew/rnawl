# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.module import Module
import sys
import os
import glob
import unittest

class ExpressionModule(Module):
    '''
    last_modify: 2019.04.11
    '''
    def __init__(self, work_id):
        super(ExpressionModule, self).__init__(work_id)
        options = [
            {'name': 'fq_type', 'type': 'string', 'default': None}, # PE SE
            {'name': 'fq_list', 'type': 'infile', 'format': 'ref_rna_v2.fastq_list'},
            {'name': 'fa_input', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'txpt2gene', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'method', 'type': 'string', 'default': None}, # rsem salmon kallisto
            {'name': 'bowtie2-mismatch-rate', 'type': 'float', 'default': 0.1},
            {'name': 'kmer', 'type': 'int', 'default': 31},
            {'name': 'strand_specific', 'type': 'bool', 'default': False},
            {'name': 'lib_type', 'type': 'string', 'default': None}, # rf (workflow forward) fr (workflow reverse)
            {'name': 'phred_quals', 'type': 'int', 'default': None}, # 33 64
        ]
        self.add_option(options)
        self.calc_tools = list()
        self.stat_tools = list()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.option('fq_list').prepare()
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(ExpressionModule, self).run()
        self.run_prepare()

    def run_prepare(self):
        self.step.add_steps('prepare')
        self.prepare = self.add_tool('tool_lab.expression.prepare')
        options = {
            'raw_fasta': self.option('fa_input'),
            'txpt2gene': self.option('txpt2gene')
        }
        self.prepare.set_options(options)
        self.prepare.on('start', self.set_step, {'start': self.step.prepare})
        self.prepare.on('end', self.set_step, {'end': self.step.prepare})
        self.prepare.on('end', self.run_tool)
        self.prepare.run()

    def run_tool(self):
        if self.option('method') == 'rsem':
            self.set_rsem()
        elif self.option('method') == 'salmon':
            self.set_salmon()
        elif self.option('method') == 'kallisto':
            self.set_kallisto()
        if len(self.calc_tools) == 1:
            self.calc_tools[0].on('end', self.run_statistics)
        else:
            self.on_rely(self.calc_tools, self.run_statistics)
        for tool in self.calc_tools:
            tool.run()

    def set_rsem(self):
        for n, sample in enumerate(self.option('fq_list').samples):
            self.step.add_steps('rsem_{}'.format(n))
            rsem = self.add_tool('tool_lab.expression.rsem')
            options = {
                'fq_type': self.option('fq_type'),
                'g2t': self.prepare.option('g2t'),
                'reference_fasta': self.prepare.option('clean_fasta'),
                'phred_quals': self.option('phred_quals'),
                'sample_name': sample,
                'bowtie2-mismatch-rate': self.option('bowtie2-mismatch-rate')
            }
            if self.option('fq_type') == 'PE':
                options.update({
                    'upstream_read': self.option('fq_list').pe1_fastq[sample],
                    'downstream_read': self.option('fq_list').pe2_fastq[sample],
                })
            elif self.option('fq_type') == 'SE':
                options.update({'upstream_read': self.option('fq_list').se_fastq[sample]})
            if self.option('strand_specific'):
                if self.option('lib_type') == 'rf':
                    options.update({'forward_prob': 0})
                elif self.option('lib_type') == 'fr':
                    options.update({'forward_prob': 1})
            rsem.set_options(options)
            rsem.on('start', self.set_step, {'start': getattr(self.step, 'rsem_{}'.format(n))})
            rsem.on('end', self.set_step, {'end': getattr(self.step, 'rsem_{}'.format(n))})
            self.calc_tools.append(rsem)

    def set_salmon(self):
        for n, sample in enumerate(self.option('fq_list').samples):
            self.step.add_steps('salmon_{}'.format(n))
            salmon = self.add_tool('tool_lab.expression.salmon')
            options = {
                'fq_type': self.option('fq_type'),
                't2g': self.prepare.option('t2g'),
                'transcripts': self.prepare.option('clean_fasta'),
                'sample_name': sample,
                'kmer': self.option('kmer')
            }
            if self.option('fq_type') == 'PE':
                options.update({
                    'mates1': self.option('fq_list').pe1_fastq[sample],
                    'mates2': self.option('fq_list').pe2_fastq[sample],
                })
            elif self.option('fq_type') == 'SE':
                options.update({'unmated_reads': self.option('fq_list').se_fastq[sample]})
            salmon.set_options(options)
            salmon.on('start', self.set_step, {'start': getattr(self.step, 'salmon_{}'.format(n))})
            salmon.on('end', self.set_step, {'end': getattr(self.step, 'salmon_{}'.format(n))})
            self.calc_tools.append(salmon)

    def set_kallisto(self):
        for n, sample in enumerate(self.option('fq_list').samples):
            self.step.add_steps('kallisto_{}'.format(n))
            kallisto = self.add_tool('tool_lab.expression.kallisto')
            options = {
                'fq_type': self.option('fq_type'),
                'fasta_file': self.prepare.option('clean_fasta'),
                'sample_name': sample,
                'kmer': self.option('kmer')
            }
            if self.option('fq_type') == 'PE':
                options.update({
                    'fastq_l_file': self.option('fq_list').pe1_fastq[sample],
                    'fastq_r_file': self.option('fq_list').pe2_fastq[sample],
                    'strand_specific': self.option('strand_specific'),
                    'stranded': self.option('lib_type')
                })
            elif self.option('fq_type') == 'SE':
                options.update({'fastq_file': self.option('fq_list').se_fastq[sample]})
            kallisto.set_options(options)
            kallisto.on('start', self.set_step, {'start': getattr(self.step, 'kallisto_{}'.format(n))})
            kallisto.on('end', self.set_step, {'end': getattr(self.step, 'kallisto_{}'.format(n))})
            self.calc_tools.append(kallisto)

    def run_statistics(self):
        def get_loc2name(exp_type, method, tools, dir_out):
            dct = {
                'T': {'rsem': 'isoforms_results', 'salmon': 'quant_sf', 'kallisto': 'abundance_tsv'},
                'G': {'rsem': 'genes_results', 'salmon': 'quant_genes_sf', 'kallisto': 'abundance_tsv'}
            }
            l2n_file = os.path.join(dir_out, 'location2name.{}.{}.tsv'.format(method, exp_type))
            open(l2n_file, 'w').writelines(['{}\t{}\n'.format(
                tool.option(dct[exp_type][method]).path, tool.option('sample_name')
            ) for tool in tools])
            return l2n_file
        for exp_type in ['T', 'G']:
            self.step.add_steps('statistics_{}'.format(exp_type))
            statistics = self.add_tool('tool_lab.expression.statistics')
            options = {
                'exp_type': exp_type,
                'method': self.option('method'),
                'loc2name': get_loc2name(exp_type, self.option('method'), self.calc_tools, self.work_dir),
                't2g': self.prepare.option('t2g')
            }
            statistics.set_options(options)
            statistics.on('start', self.set_step, {'start': getattr(self.step, 'statistics_{}'.format(exp_type))})
            statistics.on('end', self.set_step, {'end': getattr(self.step, 'statistics_{}'.format(exp_type))})
            self.stat_tools.append(statistics)
        if len(self.stat_tools) == 1:
            self.stat_tools[0].on('end', self.set_output)
        else:
            self.on_rely(self.stat_tools, self.set_output)
        for tool in self.stat_tools:
            tool.run()

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        for tool in self.stat_tools:
            for source in glob.glob(os.path.join(tool.output_dir, '*')):
                link_name = os.path.join(self.output_dir, os.path.basename(source))
                if os.path.exists(link_name):
                    os.remove(link_name)
                os.link(source, link_name)
                self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        super(ExpressionModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_rsem_pe(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'expression_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'ref_rna_v2.expression',
            'instant': False,
            'options': {
                'fq_type': 'PE',
                'fq_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/clean_data/fq.list',
                'fa_input': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/transcript.fa',
                'txpt2gene': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/method-stringtie/output/NewTranscripts/trans2gene',
                'method': 'rsem',
                'strand_specific': False,
                'lib_type': None,
                'phred_quals': 33
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_salmon_pe(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'expression_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'ref_rna_v2.expression',
            'instant': False,
            'options': {
                'fq_type': 'PE',
                'fq_list': '/mnt/lustre/users/sanger/workspace/20190617/Refrna_i-sanger_185630/HiseqQc/output/sickle_dir/fq_list.txt',
                'fa_input': '/mnt/lustre/users/sanger/workspace/20190617/Refrna_i-sanger_185630/RefrnaAssemble/output/NewTranscripts/all_transcripts.fa',
                'txpt2gene': '/mnt/lustre/users/sanger/workspace/20190617/Refrna_i-sanger_185630/RefrnaAssemble/output/NewTranscripts/trans2gene',
                'method': 'salmon',
                'strand_specific': True,
                'lib_type': None,
                'phred_quals': 33
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


    def test_kallisto_pe(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'expression_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'ref_rna_v2.expression',
            'instant': False,
            'options': {
                'fq_type': 'PE',
                'fq_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/clean_data/fq.list',
                'fa_input': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/transcript.fa',
                'txpt2gene': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/method-stringtie/output/NewTranscripts/trans2gene',
                'method': 'kallisto',
                'strand_specific': False,
                'lib_type': None,
                'phred_quals': 33
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    import sys
    suite = unittest.TestSuite()
    suite.addTests([TestFunction(sys.argv[1])])
    unittest.TextTestRunner(verbosity=2).run(suite)

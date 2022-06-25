# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, shicaiping, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class NewTranscriptsAgent(Agent):
    '''
    last_modify: 2019.01.24
    export 10 files
    '''
    def __init__(self, parent):
        super(NewTranscriptsAgent, self).__init__(parent)
        options = [
            {'name': 'merged_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'}, # from stringtie_merge
            {'name': 'gffcmp_tmap', 'type': 'infile', 'format': 'lnc_rna.tmap'}, # from gffcompare
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'add_code_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'}, # 0
            {'name': 'change_id_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'}, # 1
            {'name': 'old_genes_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'}, # 2
            {'name': 'old_trans_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'}, # 3
            {'name': 'new_genes_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'}, # 4
            {'name': 'new_trans_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'}, # 5
            {'name': 'new_trans_fa', 'type': 'outfile', 'format': 'lnc_rna.fasta'}, # 6
            {'name': 'ref_and_new_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},  # 7
            {'name': 'trans2gene', 'type': 'outfile', 'format': 'lnc_rna.common'}, # 8
            {'name': 'all_trans_fa', 'type': 'outfile', 'format': 'lnc_rna.fasta'}, # 9
        ]
        self.add_option(options)
        self.step.add_steps('new_transcripts')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.new_transcripts.start()
        self.step.update()

    def stepfinish(self):
        self.step.new_transcripts.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} - {}'.format('merged_gtf', self.option('merged_gtf').prop['path']))
        if not self.option('merged_gtf').is_set:
            raise OptionError('gffcompare GTF file must be provided')
        self.logger.debug('{} - {}'.format('gffcmp_tmap', self.option('gffcmp_tmap').prop['path']))
        if not self.option('gffcmp_tmap').is_set:
            raise OptionError('gffcompare TMAP file must be provided')
        self.logger.debug('{} - {}'.format('ref_gtf', self.option('ref_gtf').prop['path']))
        if not self.option('ref_gtf').is_set:
            raise OptionError('reference annotation GTF must be provided')
        self.logger.debug('{} - {}'.format('ref_fa', self.option('ref_fa').prop['path']))
        if not self.option('ref_fa').is_set:
            raise OptionError('reference genome FASTA must be provided')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(os.path.getsize(self.option('ref_fa').prop['path']) / 1024 ** 3 + 8)

    def end(self):
        super(NewTranscriptsAgent, self).end()

class NewTranscriptsTool(Tool):
    def __init__(self, config):
        super(NewTranscriptsTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.perl = '/miniconda2/bin/perl'
        self.gffread = 'bioinfo/rna/cufflinks-2.2.1/gffread'
        self.add_class_code_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/add_class_code.py')
        self.gtfmerge_pl = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/gtfmerge.pl')
        self.assembly_stat_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/assembly_stat.py')
        self.operate_gtf = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/operate_gtf.py')
        self.add_code_gtf = os.path.join(self.work_dir, 'add_code_merged.gtf') # 0
        self.change_id_gtf = os.path.join(self.work_dir, 'change_id_merged.gtf') # 1
        self.new_transcripts_fa = os.path.join(self.work_dir, 'new_transcripts.fa') # 6
        self.ref_and_new_gtf = os.path.join(self.work_dir, 'ref_and_new.gtf') # 7
        self.trans2gene = os.path.join(self.work_dir, 'trans2gene') # 8
        self.all_transcripts_fa = os.path.join(self.work_dir, 'all_transcripts.fa') # 9

    def run(self):
        super(NewTranscriptsTool, self).run()
        self.run_add_code()
        self.run_change_id()
        self.run_generate_merged_gtf()
        self.run_generate_changed_gtf()
        self.run_gffread_new()
        self.run_operate_gtf()
        self.run_gffread_all()
        self.set_output()
        self.end()

    def run_add_code(self):
        cmd = '{} {} '.format(self.python, self.add_class_code_py)
        cmd += '-i {} '.format(self.option('merged_gtf').prop['path'])
        cmd += '-t {} '.format(self.option('gffcmp_tmap').prop['path'])
        cmd += '-o {}'.format(self.add_code_gtf)
        cmd_name = 'run_add_code'
        self.run_code(cmd_name, cmd)

    def run_change_id(self):
        cmd = '{} {} '.format(self.perl, self.gtfmerge_pl)
        cmd += '-i {} '.format(self.add_code_gtf)
        cmd += '-compare {} '.format(self.option('gffcmp_tmap').prop['path'])
        cmd += '-ref {} '.format(self.option('ref_gtf').prop['path'])
        cmd += '-o {} '.format(self.change_id_gtf)
        cmd_name = 'run_change_id'
        self.run_code(cmd_name, cmd)

    def run_generate_merged_gtf(self):
        self.merged_dir = os.path.join(self.work_dir, 'merged')
        if not os.path.exists(self.merged_dir):
            os.mkdir(self.merged_dir)
        cmd = '{} {}'.format(self.python, self.assembly_stat_py)
        cmd += ' -transcript_file {}'.format(self.add_code_gtf)
        cmd += ' -out_new_trans {}'.format(os.path.join(self.merged_dir, 'new_transcripts.gtf'))
        cmd += ' -out_new_genes {}'.format(os.path.join(self.merged_dir, 'new_genes.gtf'))
        cmd += ' -out_old_trans {}'.format(os.path.join(self.merged_dir, 'old_transcripts.gtf'))
        cmd += ' -out_old_genes {}'.format(os.path.join(self.merged_dir, 'old_genes.gtf'))
        cmd += ' -tmapfile {}'.format(self.option('gffcmp_tmap').prop['path'])
        cmd_name = 'run_generate_merged_gtf'
        self.run_code(cmd_name, cmd)

    def run_generate_changed_gtf(self):
        self.changed_dir = os.path.join(self.work_dir, 'changed')
        if not os.path.exists(self.changed_dir):
            os.mkdir(self.changed_dir)
        cmd = '{} {} '.format(self.python, self.assembly_stat_py)
        cmd += '-transcript_file {} '.format(self.change_id_gtf)
        cmd += '-out_new_trans {} '.format(os.path.join(self.changed_dir, 'new_transcripts.gtf'))
        cmd += '-out_new_genes {} '.format(os.path.join(self.changed_dir, 'new_genes.gtf'))
        cmd += '-out_old_trans {} '.format(os.path.join(self.changed_dir, 'old_transcripts.gtf'))
        cmd += '-out_old_genes {} '.format(os.path.join(self.changed_dir, 'old_genes.gtf'))
        cmd += '-tmapfile {}'.format(self.option('gffcmp_tmap').prop['path'])
        cmd_name = 'run_generate_changed_gtf'
        self.run_code(cmd_name, cmd)

    def run_gffread_new(self):
        cmd = '{} {} '.format(self.gffread, os.path.join(self.changed_dir, 'new_transcripts.gtf'))
        cmd += '-g {} '.format(self.option('ref_fa').prop['path'])
        cmd += '-w {}'.format(self.new_transcripts_fa)
        cmd_name = 'run_gffread_new'
        self.run_code(cmd_name, cmd)

    def run_operate_gtf(self):
        cmd = '{} {} '.format(self.python, self.operate_gtf)
        cmd += '--ref {} '.format(self.option('ref_gtf').prop['path'])
        cmd += '--new {} '.format(os.path.join(self.changed_dir, 'new_transcripts.gtf'))
        cmd += '--all {} '.format(self.ref_and_new_gtf)
        cmd += '--t2g {}'.format(self.trans2gene)
        cmd_name = 'run_operate_gtf'
        self.run_code(cmd_name, cmd)

    def run_gffread_all(self):
        cmd = '{} {} '.format(self.gffread, self.ref_and_new_gtf)
        cmd += '-g {} '.format(self.option('ref_fa').prop['path'])
        cmd += '-w {}'.format(self.all_transcripts_fa)
        cmd_name = 'run_gffread_all'
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
        # 0
        add_code_gtf = os.path.join(self.output_dir, 'add_code_merged.gtf')
        if os.path.exists(add_code_gtf):
            os.remove(add_code_gtf)
        os.link(self.add_code_gtf, add_code_gtf)
        self.logger.info('succeed in linking {} to {}'.format(self.add_code_gtf, add_code_gtf))
        self.option('add_code_gtf').set_path(add_code_gtf)
        # 1
        change_id_gtf = os.path.join(self.output_dir, 'change_id_merged.gtf')
        if os.path.exists(change_id_gtf):
            os.remove(change_id_gtf)
        os.link(self.change_id_gtf, change_id_gtf)
        self.logger.info('succeed in linking {} to {}'.format(self.change_id_gtf, change_id_gtf))
        self.option('change_id_gtf').set_path(change_id_gtf)
        # 2
        old_genes_gtf = os.path.join(self.output_dir, 'old_genes.gtf')
        if os.path.exists(old_genes_gtf):
            os.remove(old_genes_gtf)
        os.link(os.path.join(self.work_dir, 'changed/old_genes.gtf'), old_genes_gtf)
        self.logger.info('succeed in linking {} to {}'.format(
            os.path.join(self.work_dir, 'changed/old_genes.gtf'),
            old_genes_gtf
        ))
        self.option('old_genes_gtf').set_path(old_genes_gtf)
        # 3
        old_trans_gtf = os.path.join(self.output_dir, 'old_transcripts.gtf')
        if os.path.exists(old_trans_gtf):
            os.remove(old_trans_gtf)
        os.link(os.path.join(self.work_dir, 'changed/old_transcripts.gtf'), old_trans_gtf)
        self.logger.info('succeed in linking {} to {}'.format(
            os.path.join(self.work_dir, 'changed/old_transcripts.gtf'),
            old_trans_gtf
        ))
        self.option('old_trans_gtf').set_path(old_trans_gtf)
        # 4
        new_genes_gtf = os.path.join(self.output_dir, 'new_genes.gtf')
        if os.path.exists(new_genes_gtf):
            os.remove(new_genes_gtf)
        os.link(os.path.join(self.work_dir, 'changed/new_genes.gtf'), new_genes_gtf)
        self.logger.info('succeed in linking {} to {}'.format(
            os.path.join(self.work_dir, 'changed/new_genes.gtf'),
            new_genes_gtf
        ))
        self.option('new_genes_gtf').set_path(new_genes_gtf)
        # 5
        new_trans_gtf = os.path.join(self.output_dir, 'new_transcripts.gtf')
        if os.path.exists(new_trans_gtf):
            os.remove(new_trans_gtf)
        os.link(os.path.join(self.work_dir, 'changed/new_transcripts.gtf'), new_trans_gtf)
        self.logger.info('succeed in linking {} to {}'.format(
            os.path.join(self.work_dir, 'changed/new_transcripts.gtf'),
            new_trans_gtf
        ))
        self.option('new_trans_gtf').set_path(new_trans_gtf)
        # 6
        ref_and_new_gtf = os.path.join(self.output_dir, 'ref_and_new.gtf')
        if os.path.exists(ref_and_new_gtf):
            os.remove(ref_and_new_gtf)
        os.link(self.ref_and_new_gtf, ref_and_new_gtf)
        self.logger.info('succeed in linking {} to {}'.format(self.ref_and_new_gtf, ref_and_new_gtf))
        self.option('ref_and_new_gtf').set_path(ref_and_new_gtf)
        # 7
        new_trans_fa = os.path.join(self.output_dir, 'new_transcripts.fa')
        if os.path.exists(new_trans_fa):
            os.remove(new_trans_fa)
        os.link(self.new_transcripts_fa, new_trans_fa)
        self.logger.info('succeed in linking {} to {}'.format(self.new_transcripts_fa, new_trans_fa))
        self.option('new_trans_fa').set_path(new_trans_fa)
        # 8
        trans2gene = os.path.join(self.output_dir, 'trans2gene')
        if os.path.exists(trans2gene):
            os.remove(trans2gene)
        os.link(self.trans2gene, trans2gene)
        self.logger.info('succeed in linking {} to {}'.format(self.trans2gene, trans2gene))
        self.option('trans2gene').set_path(trans2gene)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        # 9
        all_trans_fa = os.path.join(self.output_dir, 'all_transcripts.fa')
        if os.path.exists(all_trans_fa):
            os.remove(all_trans_fa)
        os.link(self.all_transcripts_fa, all_trans_fa)
        self.logger.info('succeed in linking {} to {}'.format(self.all_transcripts_fa, all_trans_fa))
        self.option('all_trans_fa').set_path(all_trans_fa)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_hsa(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'new_transcripts_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.assemble.new_transcripts',
            'instant': False,
            'options': {
                'merged_gtf': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/assemble/merged.gtf',
                'gffcmp_tmap': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/assemble/gffcompare.tmap',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
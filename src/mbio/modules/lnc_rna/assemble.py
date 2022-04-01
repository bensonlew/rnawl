# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, shicaiping, qinjincheng'

from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import os
import sys
import shutil
import unittest

class AssembleModule(Module):
    '''
    last_modify: 2019.04.23
    '''
    def __init__(self, work_id):
        super(AssembleModule, self).__init__(work_id)
        options = [
            {'name': 'assemble_method', 'type': 'string', 'default': ''},
            {'name': 'fr_stranded', 'type': 'string', 'default': 'fr-unstranded'},
            {'name': 'strand_direct', 'type': 'string', 'default': 'none'},
            {'name': 'bam_loc', 'type': 'infile', 'format': 'lnc_rna.loc2name'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            # stringtie arguments
            {'name': 'min_coverage', 'type': 'int', 'default': 3},
            {'name': 'min_read', 'type': 'int', 'default': 5},
            {'name': 'min_cov', 'type': 'int', 'default': 5},
            {'name': 'min_tpm', 'type': 'int', 'default': 1},
            {'name': 'min_iso', 'type': 'float', 'default': 0.4},
            # cufflinks arguments
            {'name': 'F', 'type': 'float', 'default': 0.1},
            {'name': 'fpkm_cut', 'type': 'int', 'default': 1},
            {'name': 'min_isoform_fraction', 'type': 'float', 'default': 0.4},
            # output
            {'name': 'merged_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'gffcmp_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'change_id_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'new_gene_gtf', 'type': 'outfile', 'format': 'ref_rna_v2.gtf'},
            {'name': 'new_transcripts_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'new_transcripts_fa', 'type': 'outfile', 'format': 'lnc_rna.fasta'},
            {'name': 'ref_and_new_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'all_transcripts_fa', 'type': 'outfile', 'format': 'lnc_rna.fasta'},
            {'name': 'trans2gene', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps(
            'stringtie', 'stringtie_merge',
            'cufflinks', 'cuffmerge',
            'gffcompare', 'new_transcripts', 'assemble_stat'
        )
        self.tools = list()
        self.sum_tools = list()

    def check_options(self):
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run(self):
        super(AssembleModule, self).run()
        if self.option('assemble_method').lower() == '':
            sys.exit(1)
        elif self.option('assemble_method').lower() == 'stringtie':
            self.run_stringtie()
        elif self.option('assemble_method').lower() == 'cufflinks':
            self.run_cufflinks()

    def run_stringtie(self):
        self.step.stringtie.start()
        self.step.update()
        name2loc_dict = self.option('bam_loc').prop['locations']
        for n, sample_name in enumerate(self.option('bam_loc').prop['names']):
            self.step.add_steps('stringtie_{}'.format(n))
            stringtie = self.add_tool('lnc_rna.assemble.stringtie')
            stringtie.set_options({
                'fr_stranded': self.option('fr_stranded'),
                'strand_direct': self.option('strand_direct'),
                'sample_name': sample_name,
                'sample_bam': name2loc_dict[sample_name],
                'ref_gtf': self.option('ref_gtf'),
                'min_coverage': self.option('min_coverage'),
                "min_read": self.option("min_read"),
                'ref_fa': self.option('ref_fa'),
            })
            step = getattr(self.step, 'stringtie_{}'.format(n))
            step.start()
            stringtie.on('end', self.finish_update, 'stringtie_{}'.format(n))
            self.tools.append(stringtie)
            self.sum_tools.append(stringtie)
        if len(self.tools) == 1:
            self.tools[0].on('end', self.run_stringtie_merge)
        else:
            self.on_rely(self.tools, self.run_stringtie_merge)
        for tool in self.tools:
            tool.run()
        else:
            self.step.stringtie.finish()
            self.step.update()

    def run_stringtie_merge(self):
        self.step.stringtie_merge.start()
        self.step.update()
        gtf_list = os.path.join(self.work_dir, 'gtf.list')
        with open(gtf_list, 'w') as f:
            for tool in self.tools:
                sample_gtf = tool.option('sample_gtf').prop['path']
                if os.path.getsize(sample_gtf) > 0:
                    f.write('{}\n'.format(sample_gtf))
        self.stringtie_merge = self.add_tool('lnc_rna.assemble.stringtie_merge')
        self.stringtie_merge.set_options({
            'gtf_list': gtf_list,
            'ref_gtf': self.option('ref_gtf'),
            'min_cov': self.option('min_cov'),
            'min_tpm': self.option('min_tpm'),
            'min_iso': self.option('min_iso'),
            'ref_fa': self.option('ref_fa'),
        })
        self.stringtie_merge.on('end', self.run_gffcompare, 'stringtie')
        self.stringtie_merge.run()
        self.sum_tools.append(self.stringtie_merge)
        self.step.stringtie_merge.finish()
        self.step.update()

    def run_cufflinks(self):
        self.step.cufflinks.start()
        self.step.update()
        name2loc_dict = self.option('bam_loc').prop['locations']
        for n, sample_name in enumerate(self.option('bam_loc').prop['names']):
            self.step.add_steps('cufflinks_{}'.format(n))
            cufflinks = self.add_tool('lnc_rna.assemble.cufflinks')
            cufflinks.set_options({
                'fr_stranded': self.option('fr_stranded'),
                'strand_direct': self.option('strand_direct'),
                'sample_name': sample_name,
                'sample_bam': name2loc_dict[sample_name],
                'ref_gtf': self.option('ref_gtf'),
                'F': self.option('F'),
                'fpkm_cut': self.option('fpkm_cut'),
                'ref_fa': self.option('ref_fa'),
            })
            step = getattr(self.step, 'cufflinks_{}'.format(n))
            step.start()
            cufflinks.on('end', self.finish_update, 'cufflinks_{}'.format(n))
            self.tools.append(cufflinks)
            self.sum_tools.append(cufflinks)
        if len(self.tools) == 1:
            self.tools[0].on('end', self.run_cuffmerge)
        else:
            self.on_rely(self.tools, self.run_cuffmerge)
        for tool in self.tools:
            tool.run()
        else:
            self.step.stringtie.finish()
            self.step.update()

    def run_cuffmerge(self):
        self.step.cuffmerge.start()
        self.step.update()
        gtf_list = os.path.join(self.work_dir, 'gtf.list')
        with open(gtf_list, 'w') as f:
            for tool in self.tools:
                sample_gtf = tool.option('sample_gtf').prop['path']
                if os.path.getsize(sample_gtf) > 0:
                    f.write('{}\n'.format(sample_gtf))
        self.cuffmerge = self.add_tool('lnc_rna.assemble.cuffmerge')
        self.cuffmerge.set_options({
            'gtf_list': gtf_list,
            'ref_gtf': self.option('ref_gtf'),
            'ref_fa': self.option('ref_fa'),
            'min_isoform_fraction': self.option('min_isoform_fraction'),
        })
        self.cuffmerge.on('end', self.run_gffcompare, 'cufflinks')
        self.cuffmerge.run()
        self.sum_tools.append(self.cuffmerge)
        self.step.cuffmerge.finish()
        self.step.update()

    def run_gffcompare(self, event):
        self.step.gffcompare.start()
        self.step.update()
        if not event['data']:
            sys.exit(2)
        elif event['data'] == 'stringtie':
            self.merged_gtf = self.stringtie_merge.option('merged_gtf').prop['path']
        elif event['data'] == 'cufflinks':
            self.merged_gtf = self.cuffmerge.option('merged_gtf').prop['path']
        self.gffcompare = self.add_tool('lnc_rna.assemble.gffcompare')
        self.gffcompare.set_options({
            'merged_gtf': self.merged_gtf,
            'ref_gtf': self.option('ref_gtf'),
        })
        self.gffcompare.on('end', self.run_new_transcripts)
        self.gffcompare.run()
        self.sum_tools.append(self.gffcompare)
        self.step.gffcompare.finish()
        self.step.update()

    def run_new_transcripts(self):
        self.step.new_transcripts.start()
        self.step.update()
        if os.path.getsize(self.gffcompare.option('gffcmp_tmap').prop['path']) > 0:
            self.gffcmp_tmap = self.gffcompare.option('gffcmp_tmap').prop['path']
        self.new_transcripts = self.add_tool('lnc_rna.assemble.new_transcripts')
        self.new_transcripts.set_options({
            'merged_gtf': self.merged_gtf,
            'gffcmp_tmap': self.gffcmp_tmap,
            'ref_fa': self.option('ref_fa'),
            'ref_gtf': self.option('ref_gtf'),
        })
        self.new_transcripts.on('end', self.run_assemble_stat)
        self.new_transcripts.run()
        self.sum_tools.append(self.new_transcripts)
        self.step.new_transcripts.finish()
        self.step.update()

    def run_assemble_stat(self):
        self.step.assemble_stat.start()
        self.step.update()
        self.all_files_dir = os.path.join(self.work_dir, 'AllFilesDir')
        if os.path.exists(self.all_files_dir):
            shutil.rmtree(self.all_files_dir)
        os.mkdir(self.all_files_dir)
        for tool in self.sum_tools:
            self.linkdir(tool.output_dir, self.all_files_dir)
        self.assemble_stat = self.add_tool('lnc_rna.assemble.assemble_stat')
        self.assemble_stat.set_options({
            'all_files_dir': self.all_files_dir,
            'assemble_method': self.option('assemble_method'),
        })
        self.assemble_stat.on('end', self.set_output)
        self.assemble_stat.run()
        self.step.assemble_stat.finish()
        self.step.update()

    def linkdir(self, dir_in, dir_out):
        allfiles = os.listdir(dir_in)
        if not os.path.exists(dir_out):
            os.mkdir(dir_out)
        oldfiles = [os.path.join(dir_in, i) for i in allfiles]
        newfiles = [os.path.join(dir_out, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                elif os.path.isdir(newfile):
                    shutil.rmtree(newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                if os.path.exists(newfiles[i]):
                    os.remove(newfiles[i])
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                if os.path.exists(dir_out):
                    os.remove(dir_out)
                os.link(oldfiles[i], dir_out)

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        dir_list = ['Gffcompare', 'NewTranscripts', 'Statistics']
        if self.option('assemble_method').lower() == '':
            sys.exit(3)
        elif self.option('assemble_method').lower() == 'stringtie':
            single_dir, merged_dir = 'Stringtie', 'StringtieMerge'
        elif self.option('assemble_method').lower() == 'cufflinks':
            single_dir, merged_dir = 'Cufflinks', 'Cuffmerge'
        dir_list.extend([single_dir, merged_dir])
        for subdir in dir_list:
            os.mkdir(os.path.join(self.output_dir, subdir))
        for f in os.listdir(self.all_files_dir):
            source = os.path.join(self.all_files_dir, f)
            if f.endswith('_out.gtf') or f.endswith('_out.fa'):
                link_name = os.path.join(self.output_dir, single_dir, f)
            elif f.endswith('merged.gtf') or f.endswith('merged.fa'):
                link_name = os.path.join(self.output_dir, merged_dir, f)
            elif f.startswith('gffcmp'):
                link_name = os.path.join(self.output_dir, 'Gffcompare', f)
            elif f.startswith('old_') or f.startswith('new_'):
                link_name = os.path.join(self.output_dir, 'NewTranscripts', f)
            elif f == 'ref_and_new.gtf' or f == 'all_transcripts.fa' or f == 'trans2gene':
                link_name = os.path.join(self.output_dir, 'NewTranscripts', f)
            else:
                continue
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        for f in os.listdir(self.assemble_stat.output_dir):
            source = os.path.join(self.assemble_stat.output_dir, f)
            link_name = os.path.join(self.output_dir, 'Statistics', f)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('merged_gtf').set_path(os.path.join(self.output_dir, merged_dir, 'merged.gtf'))
        self.option('gffcmp_gtf').set_path(os.path.join(self.output_dir, 'Gffcompare', 'gffcmp.gtf'))
        self.option('change_id_gtf').set_path(os.path.join(self.output_dir, merged_dir, 'change_id_merged.gtf'))
        self.option('new_gene_gtf').set_path(os.path.join(self.output_dir, 'NewTranscripts', 'new_genes.gtf'))
        self.option('new_transcripts_gtf').set_path(os.path.join(self.output_dir, 'NewTranscripts', 'new_transcripts.gtf'))
        self.option('new_transcripts_fa').set_path(os.path.join(self.output_dir, 'NewTranscripts', 'new_transcripts.fa'))
        self.option('ref_and_new_gtf').set_path(os.path.join(self.output_dir, 'NewTranscripts', 'ref_and_new.gtf'))
        self.option('all_transcripts_fa').set_path(os.path.join(self.output_dir, 'NewTranscripts', 'all_transcripts.fa'))
        self.option('trans2gene').set_path(os.path.join(self.output_dir, 'NewTranscripts', 'trans2gene'))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        super(AssembleModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run script to do test.
    '''

    def test_stringtie(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'assemble_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'lnc_rna.assemble',
            'instant': False,
            'options': {
                'assemble_method': 'stringtie',
                'bam_loc': '/mnt/ilustre/users/sanger-dev/workspace/20190318/Single_RnaseqMapping_star_test6499/RnaseqMapping/output/bam_loc.txt',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_cufflinks(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'assemble_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'lnc_rna.assemble',
            'instant': False,
            'options': {
                'assemble_method': 'cufflinks',
                'bam_loc': '/mnt/ilustre/users/sanger-dev/workspace/20190318/Single_RnaseqMapping_star_test6499/RnaseqMapping/output/bam_loc.txt',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_stringtie')])
    unittest.TextTestRunner(verbosity=2).run(suite)
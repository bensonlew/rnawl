# -*- coding: utf-8 -*-
# __author__ = 'liubinxu, qinjincheng'

from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import os
import shutil
import unittest

class AnnotClassModule(Module):
    '''
    last_modify: 2019.02.28b
    '''
    def __init__(self, work_id):
        super(AnnotClassModule, self).__init__(work_id)
        options = [
            # essential options for determining the process
            {'name': 'db', 'type': 'string', 'default': 'nr,swissprot,kegg,eggnog,pfam,go'},
            # ['ref', 'new']
            {'name': 'type', 'type': 'string', 'default': ''},
            # input files and related options
            {'name': 'blast_nr_xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'}, # nr
            {'name': 'blast_swissprot_xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'}, # swissprot
            {'name': 'blast_eggnog_xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'}, # eggnog
            {'name': 'blast_kegg_xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'}, # kegg
            {'name': 'known_ko', 'type': 'infile', 'format': 'lnc_rna.common'},
            # ['All', 'Animals', 'Plants', 'Protists']
            {'name': 'taxonomy', 'type': 'string', 'default': ''},
            # {'ref': 'yellow', 'new': 'green'}
            {'name': 'link_bgcolor', 'type': 'string', 'default': ''},
            # {'ref': '#FFFF00', 'new': '#00CD00'}
            {'name': 'png_bgcolor', 'type': 'string', 'default': ''},
            {'name': 'blast2go_annot', 'type': 'infile', 'format': 'lnc_rna.common'}, # go
            {'name': 'pfam_domain', 'type': 'infile', 'format': 'lnc_rna.common'}, # pfam
            # essential options for doing statistics on annotation
            {'name': 'gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'g2t2p', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'gene2trans', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'des', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'des_type', 'type': 'string', 'default': ''},
            {'name': 'enterz', 'type': 'infile', 'format': 'lnc_rna.common'},
            # output files
            {'name': 'gene_go_list', 'type': 'outfile', 'format': 'ref_rna_v2.go_list'},
            {'name': 'gene_kegg_table', 'type': 'outfile', 'format': 'ref_rna_v2.kegg_table'},
            {'name': 'kegg_table', 'type': 'outfile', 'format': 'ref_rna_v2.kegg_table'},
            {'name': 'gene_go_level_2', 'type': 'outfile', 'format': 'ref_rna_v2.level2'},
            {"name": "nr_version", "type": "string", "default": ""},
            {"name": "swissprot_version", "type": "string", "default": ""},
            {"name": "eggnog_version", "type": "string", "default": ""},
            {"name": "string_version", "type": "string", "default": ""},
            {"name": "cog_version", "type": "string", "default": ""},
            {"name": "kegg_version", "type": "string", "default": "202003"},
            {"name": "kegg_subtax1", "type": "string", "default": None},
            {"name": "kegg_species", "type": "string", "default": "2019"},
            {"name": "kegg_subtax2", "type": "string", "default": None},
            {"name": "pir_version", "type": "string", "default":"2019"},
            {"name": "ncbi_taxonmy_version", "type": "string", "default":"2019"},
            {"name": "go_version", "type": "string", "default": "2019"},
            {"name": "pfam_version", "type": "string", "default": "32"},

        ]
        self.add_option(options)
        self.tools = list()
        self.step.add_steps(
            'nr_annot',
            'swissprot_annot',
            'eggnog_annot',
            'kegg_annot',
            'pfam_annot',
            'go_annot',
            'annot_file',
            'annot_stat',
            'annot_query'
        )

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} = {}'.format('db', self.option('db')))
        if self.option('db') == '':
            raise OptionError('required database must be specified')
        self.logger.debug('{} = {}'.format('type', self.option('type')))
        if self.option('type') == '':
            raise OptionError('type of classification of annotation must be specified')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_step(self, event):
        if 'start' in event['data']:
            event['data']['start'].start()
        if 'end' in event['data']:
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        if 'nr' in self.option('db'):
            self.set_nr_annot()
            self.tools.append(self.nr_annot)
        if 'swissprot' in self.option('db'):
            self.set_swissprot_annot()
            self.tools.append(self.swissprot_annot)
        if 'eggnog' in self.option('db'):
            self.set_eggnog_annot()
            self.tools.append(self.eggnog_annot)
        if 'kegg' in self.option('db'):
            self.set_kegg_annot()
            self.tools.append(self.kegg_annot)
        if 'go' in self.option('db'):
            self.set_go_annot()
            self.tools.append(self.go_annot)
        if 'pfam' in self.option('db'):
            self.set_pfam_annot()
            self.tools.append(self.pfam_annot)
        super(AnnotClassModule, self).run()
        if len(self.tools) == 1:
            self.tools[0].on('end', self.run_annot_file)
        else:
            self.on_rely(self.tools, self.run_annot_file)
        for tool in self.tools:
            tool.run()

    def set_nr_annot(self):
        self.nr_annot = self.add_tool('lnc_rna.annotation.nr_annot')
        options = {
            'blast_xml': self.option('blast_nr_xml')
        }
        self.nr_annot.set_options(options)
        self.nr_annot.on('start', self.set_step, {'start': self.step.nr_annot})
        self.nr_annot.on('end', self.set_step, {'end': self.step.nr_annot})
        self.nr_annot.on('end', self.set_output, 'nr_annot')

    def set_swissprot_annot(self):
        self.swissprot_annot = self.add_tool('lnc_rna.annotation.swissprot_annot')
        options = {
            'blast_xml': self.option('blast_swissprot_xml')
        }
        self.swissprot_annot.set_options(options)
        self.swissprot_annot.on('start', self.set_step, {'start': self.step.swissprot_annot})
        self.swissprot_annot.on('end', self.set_step, {'end': self.step.swissprot_annot})
        self.swissprot_annot.on('end', self.set_output, 'swissprot_annot')

    def set_eggnog_annot(self):
        self.eggnog_annot = self.add_tool('lnc_rna.annotation.eggnog_annot')
        options = {
            'blast_xml': self.option('blast_eggnog_xml'),
            'eggnog_version': self.option('eggnog_version')
        }
        self.eggnog_annot.set_options(options)
        self.eggnog_annot.on('start', self.set_step, {'start': self.step.eggnog_annot})
        self.eggnog_annot.on('end', self.set_step, {'end': self.step.eggnog_annot})
        self.eggnog_annot.on('end', self.set_output, 'eggnog_annot')

    def set_kegg_annot(self):
        self.kegg_annot = self.add_tool('lnc_rna.annotation.kegg_annot')
        options = {
            'blast_xml': self.option('blast_kegg_xml'),
            'known_ko': self.option('known_ko'),
            'taxonomy': self.option('taxonomy'),
            'link_bgcolor': self.option('link_bgcolor'),
            'png_bgcolor': self.option('png_bgcolor'),
            'kegg_version': self.option('kegg_version'),
            'kegg_subtax1': self.option('kegg_subtax1'),
            'kegg_subtax2': self.option('kegg_subtax2'),
            'kegg_species': self.option('kegg_species')
        }
        self.kegg_annot.set_options(options)
        self.kegg_annot.on('start', self.set_step, {'start': self.step.kegg_annot})
        self.kegg_annot.on('end', self.set_step, {'end': self.step.kegg_annot})
        self.kegg_annot.on('end', self.set_output, 'kegg_annot')

    def set_go_annot(self):
        self.go_annot = self.add_tool('lnc_rna.annotation.go_annot')
        options = {
            'blast2go_annot': self.option('blast2go_annot'),
            'pir_version': self.option('pir_version'),
            'go_version': self.option('go_version')
        }
        self.go_annot.set_options(options)
        self.go_annot.on('start', self.set_step, {'start': self.step.go_annot})
        self.go_annot.on('end', self.set_step, {'end': self.step.go_annot})
        self.go_annot.on('end', self.set_output, 'go_annot')

    def set_pfam_annot(self):
        self.pfam_annot = self.add_tool('lnc_rna.annotation.pfam_annot')
        options = {
            'file_in': self.option('pfam_domain'),
            'pfam_version': self.option('pfam_version')
        }
        self.pfam_annot.set_options(options)
        self.pfam_annot.on('start', self.set_step, {'start': self.step.pfam_annot})
        self.pfam_annot.on('end', self.set_step, {'end': self.step.pfam_annot})
        self.pfam_annot.on('end', self.set_output, 'pfam_annot')

    def run_annot_file(self):
        self.annot_file = self.add_tool('lnc_rna.annotation.annot_file')
        if self.option('type') == 'ref':
            options = {
                'gtf': self.option('gtf'),
                'g2t2p': self.option('g2t2p')
            }
        elif self.option('type') == 'new':
            options = {
                'gene2trans': self.option('gene2trans')
            }
        self.annot_file.set_options(options)
        self.annot_file.on('start', self.set_step, {'start': self.step.annot_file})
        self.annot_file.on('end', self.set_step, {'end': self.step.annot_file})
        self.annot_file.on('end', self.set_output, 'annot_file')
        self.annot_file.on('end', self.run_annot_stat)
        self.annot_file.run()

    def run_annot_stat(self):
        self.annot_stat = self.add_tool('lnc_rna.annotation.annot_stat')
        options = {
            'database': self.option('db'),
            'gene2trans': self.annot_file.option('i2u')
        }
        if 'nr' in self.option('db'):
            options['nr_xml'] = self.option('blast_nr_xml')
        if 'swissprot' in self.option('db'):
            options['swissprot_xml'] = self.option('blast_swissprot_xml')
        if 'eggnog' in self.option('db'):
            options['eggnog_cog_table'] = self.eggnog_annot.option('cog_table')
        if 'kegg' in self.option('db'):
            options['kegg_xml'] = self.option('blast_kegg_xml')
            options['kegg_anno_table'] = self.kegg_annot.option('kegg_table')
            options['link_bgcolor'] = self.option('link_bgcolor')
            options['png_bgcolor'] = self.option('png_bgcolor')
            options['taxonomy'] = self.option('taxonomy')
            options['kegg_version'] = self.option('kegg_version')
            if self.option('type') == 'ref':
                options['known_ko'] = self.option('known_ko')
        if 'go' in self.option('db'):
            options['blast2go_annot'] = self.option('blast2go_annot')
            options['gos_list'] = self.go_annot.option('go_list')
        if 'pfam' in self.option('db'):
            options['pfam_domain'] = self.option('pfam_domain')
        self.annot_stat.set_options(options)
        self.annot_stat.on('start', self.set_step, {'start': self.step.annot_stat})
        self.annot_stat.on('end', self.set_step, {'end': self.step.annot_stat})
        self.annot_stat.on('end', self.set_output, 'annot_stat')
        self.annot_stat.on('end', self.run_annot_query)
        self.annot_stat.run()

    def run_annot_query(self):
        self.annot_query = self.add_tool('lnc_rna.annotation.annot_query')
        options = {
            'gene2trans': self.annot_file.option('i2u'),
            'des': self.option('des'),
            'des_type': self.option('des_type'),
            'enterz': self.option('enterz'),
            'kegg_version': self.option('kegg_version')
        }
        if 'nr' in self.option('db'):
            options['blast_nr_table'] = self.nr_annot.option('blast_table')
        if 'swissprot' in self.option('db'):
            options['blast_swissprot_table'] = self.swissprot_annot.option('blast_table')
        if 'eggnog' in self.option('db'):
            options['cog_list'] = self.eggnog_annot.option('cog_table')
        if 'kegg' in self.option('db'):
            options['kegg_table'] = self.kegg_annot.option('kegg_table')
        if 'pfam' in self.option('db'):
            options['pfam_domain'] = self.pfam_annot.option('file_out')
        if 'go' in self.option('db'):
            options['gos_list'] = self.go_annot.option('go_list')
        self.annot_query.set_options(options)
        self.annot_query.on('start', self.set_step, {'start': self.step.annot_query})
        self.annot_query.on('end', self.set_step, {'end': self.step.annot_query})
        self.annot_query.on('end', self.set_output, 'annot_query')
        self.annot_query.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'annot_stat':
            self.linkdir(obj.output_dir, 'anno_stat')
            if 'eggnog' in self.option('db'):
                os.link(
                    os.path.join(obj.work_dir, 'cog_summary.xls'), os.path.join(self.output_dir, 'cog/cog_summary.xls')
                )
            if 'kegg' in self.option('db'):
                self.option('gene_kegg_table').set_path(
                    os.path.join(self.output_dir, 'anno_stat/kegg_stat/gene_kegg_table.xls')
                )
            if 'go' in self.option('db'):
                self.option('gene_go_list').set_path(
                    os.path.join(self.output_dir, 'anno_stat/go_stat/gene_gos.list'))
                self.option('gene_go_level_2').set_path(
                    os.path.join(self.output_dir, 'anno_stat/go_stat/gene_go12level_statistics.xls')
                )
        elif event['data'] == 'nr_annot':
            self.linkdir(obj.output_dir, 'nr')
        elif event['data'] == 'swissprot_annot':
            self.linkdir(obj.output_dir, 'swissprot')
        elif event['data'] == 'eggnog_annot':
            self.linkdir(obj.output_dir, 'cog')
        elif event['data'] == 'kegg_annot':
            self.linkdir(obj.output_dir, 'kegg')
        elif event['data'] == 'pfam_annot':
            self.linkdir(obj.output_dir, 'pfam')
            if os.path.exists(os.path.join(self.output_dir, 'pfam_domain')):
                os.remove(os.path.join(self.output_dir, 'pfam_domain'))
            os.link(os.path.join(obj.output_dir, 'pfam.filter.xls'), os.path.join(self.output_dir, 'pfam_domain'))
        elif event['data'] == 'go_annot':
            self.linkdir(obj.output_dir, 'go')
        elif event['data'] == 'annot_file':
            if os.path.exists(os.path.join(self.output_dir, 'tran2gene.txt')):
                os.remove(os.path.join(self.output_dir, 'tran2gene.txt'))
            os.link(os.path.join(obj.output_dir, 'i2u'), os.path.join(self.output_dir, 'tran2gene.txt'))
        elif event['data'] == 'annot_query':
            if os.path.exists(os.path.join(self.output_dir, 'anno_stat/all_anno_detail.xls')):
                os.remove(os.path.join(self.output_dir, 'anno_stat/all_anno_detail.xls'))
            os.link(
                os.path.join(self.annot_query.output_dir, 'all_anno_detail.xls'),
                os.path.join(self.output_dir, 'anno_stat/all_anno_detail.xls')
            )
            if os.path.exists(os.path.join(self.output_dir, 'anno_stat/gene_anno_detail.xls')):
                os.remove(os.path.join(self.output_dir, 'anno_stat/gene_anno_detail.xls'))
            self.end()

    def linkdir(self, olddir, newname):
        newdir = os.path.join(self.output_dir, newname)
        if os.path.exists(newdir):
            shutil.rmtree(newdir)
        os.mkdir(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            else:
                new1 = os.path.join(newdir, os.path.basename(oldfiles[i]))
                os.system('cp -r {} {}'.format(oldfiles[i], new1))

    def end(self):
        super(AnnotClassModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annot_class_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'lnc_rna.annot_class',
            'instant': False,
            'options': {
                'db': 'nr,swissprot,kegg,eggnog,pfam,go',
                'taxonomy': 'Animals',
                'blast_nr_xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/nr/blast.filter.xml',
                'blast_swissprot_xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/swissprot/blast.filter.xml',
                'blast_eggnog_xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/eggnog/blast.filter.xml',
                'blast_kegg_xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/kegg/blast.filter.xml',
                'pfam_domain': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/pfam/pfam.filter.xls',
                'blast2go_annot': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/go/blast2go.filter.xls',
                'g2t2p': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/Annotation_v2/g2t2p',
                'gtf': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'des': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/biomart/Homo_sapiens.GRCh38.biomart_gene.txt',
                'des_type': 'type1',
                'enterz': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/NCBI/Homo_sapiens.GRCh38.biomart_enterz.txt',
                'type': 'ref',
                'known_ko': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/KEGG/Homo_sapiens.GRCh38.pathway',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()

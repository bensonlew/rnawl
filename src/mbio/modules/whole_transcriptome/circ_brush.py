# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun,qinjincheng'

import os
import shutil
import unittest

from biocluster.module import Module


class CircBrushModule(Module):
    def __init__(self, work_id):
        super(CircBrushModule, self).__init__(work_id)
        CIRC_METHOD = ('ciri2,find_circ', 'ciri2', 'find_circ', 'circ_finder', 'circexplorer2')
        options = [
            {'name': 'circ_method', 'type': 'string', 'default': CIRC_METHOD[0]},
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'annotate', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            {'name': 'organism_name', 'type': 'string', 'default': None},
            {'name': 'junction_reads', 'type': 'int', 'default': 2},
            {'name': 'circrna_length', 'type': 'int', 'default': 100000}
        ]
        self.add_option(options)
        self.modules = list()

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(CircBrushModule, self).run()
        if self.option('circ_method') in ('ciri2,find_circ', 'ciri2', 'find_circ', 'circexplorer2'):
            self.run_bwaindex()
        if self.option('circ_method') in ('circ_finder'):
            self.run_circ_finder()

    def run_bwaindex(self):
        if self.option('circ_method') == 'ciri2,find_circ':
            self.bwaindex = self.add_tool('whole_transcriptome.circrna.bwaindex')
            self.bwaindex.set_options({'genome': self.option('genome')})
            self.bwaindex.on('end', self.run_ciriaddfindcirc)
            self.bwaindex.run()
        if self.option('circ_method') == 'ciri2':
            self.bwaindex = self.add_tool('whole_transcriptome.circrna.bwaindex')
            self.bwaindex.set_options({'genome': self.option('genome')})
            self.bwaindex.on('end', self.run_ciri2)
            self.bwaindex.run()
        if self.option('circ_method') == 'find_circ':
            self.bwaindex = self.add_tool('whole_transcriptome.circrna.bwaindex')
            self.bwaindex.set_options({'genome': self.option('genome')})
            self.bwaindex.on('end', self.run_find_circ)
            self.bwaindex.run()
        if self.option('circ_method') == 'circexplorer2':
            self.bwaindex = self.add_tool('whole_transcriptome.circrna.bwaindex')
            self.bwaindex.set_options({'genome': self.option('genome')})
            self.bwaindex.on('end', self.run_circexplorer)
            self.bwaindex.run()

    def run_circ_finder(self):
        is_se, sp2fqs_dict = read_fastq_dir(self.option('fastq_dir').path)
        for sample, fastq_list in sp2fqs_dict.items():
            circfinder = self.add_module('whole_transcriptome.circrna.circfinder')
            opts = {'genome': self.option('genome').path,
                    'annotate': self.option('annotate').path,
                    'sample': sample,
                    'junction_reads': self.option('junction_reads'),
                    'circrna_length': self.option('circrna_length')
                    }
            if is_se:
                opts.update({'fq1': fastq_list[0]})
            else:
                opts.update({'fq1': fastq_list[0], 'fq2': fastq_list[1]})
            circfinder.set_options(opts)
            self.modules.append(circfinder)
        else:
            self.on_rely(self.modules, self.run_count)
        for module in self.modules:
            module.run()

    def run_circexplorer(self):
        is_se, sp2fqs_dict = read_fastq_dir(self.option('fastq_dir').path)
        for sample, fastq_list in sp2fqs_dict.items():
            circexplorer = self.add_module('whole_transcriptome.circrna.circexplorer')
            opts = {'genome': self.option('genome').path,
                    'annotate': self.option('annotate').path,
                    'sample': sample,
                    'junction_reads': self.option('junction_reads'),
                    'circrna_length': self.option('circrna_length')
                    }
            if is_se:
                opts.update({'fq1': fastq_list[0]})
            else:
                opts.update({'fq1': fastq_list[0], 'fq2': fastq_list[1]})
            circexplorer.set_options(opts)
            self.modules.append(circexplorer)
        else:
            self.on_rely(self.modules, self.run_count)
        for module in self.modules:
            module.run()

    def run_ciri2(self):
        is_se, sp2fqs_dict = read_fastq_dir(self.option('fastq_dir').path)
        for sample, fastq_list in sp2fqs_dict.items():
            ciri2 = self.add_module('whole_transcriptome.circrna.ciri')
            opts = {'genome': self.option('genome').path,
                    'annotate': self.option('annotate').path,
                    'sample': sample,
                    'junction_reads': self.option('junction_reads'),
                    'circrna_length': self.option('circrna_length')
                    }
            if is_se:
                opts.update({'fq1': fastq_list[0]})
            else:
                opts.update({'fq1': fastq_list[0], 'fq2': fastq_list[1]})
            ciri2.set_options(opts)
            self.modules.append(ciri2)
        else:
            self.on_rely(self.modules, self.run_count)
        for module in self.modules:
            module.run()

    def run_find_circ(self):
        is_se, sp2fqs_dict = read_fastq_dir(self.option('fastq_dir').path)
        for sample, fastq_list in sp2fqs_dict.items():
            findcircbwa = self.add_module('whole_transcriptome.circrna.findcirc')
            opts = {'genome': self.option('genome').path,
                    'annotate': self.option('annotate').path,
                    'sample': sample,
                    'junction_reads': self.option('junction_reads'),
                    'circrna_length': self.option('circrna_length')
                    }
            if is_se:
                opts.update({'fq1': fastq_list[0]})
            else:
                opts.update({'fq1': fastq_list[0], 'fq2': fastq_list[1]})
            findcircbwa.set_options(opts)
            self.modules.append(findcircbwa)
        else:
            self.on_rely(self.modules, self.run_count)
        for module in self.modules:
            module.run()

    def run_ciriaddfindcirc(self):
        is_se, sp2fqs_dict = read_fastq_dir(self.option('fastq_dir').path)
        for sample, fastq_list in sp2fqs_dict.items():
            ciriaddfindcirc = self.add_module('whole_transcriptome.circrna.ciriaddfindcirc')
            opts = {
                'genome': self.option('genome').path,
                'annotate': self.option('annotate').path,
                'sample': sample,
                'junction_reads': self.option('junction_reads'),
                'circrna_length': self.option('circrna_length')
            }
            if is_se:
                opts.update({'fq1': fastq_list[0]})
            else:
                opts.update({'fq1': fastq_list[0], 'fq2': fastq_list[1]})
            ciriaddfindcirc.set_options(opts)
            self.modules.append(ciriaddfindcirc)
        else:
            self.on_rely(self.modules, self.run_count)
        for module in self.modules:
            module.run()

    def run_count(self):
        self.count = self.add_tool('whole_transcriptome.circrna.circcount')
        circ_list = os.path.join(self.work_dir, 'circ_id.list')
        with open(circ_list, 'w') as handle:
            for module in self.modules:
                handle.write('{}\t{}\n'.format(module.option('sample'), module.rpm.option('RPM').path))
        opts = {
            'list_id': circ_list,

        }

        self.count.set_options(opts)
        self.count.on('end', self.run_circrpm)
        self.count.run()

    def run_circrpm(self):
        self.circrpm = self.add_tool('whole_transcriptome.circrna.circrpm')
        circ_list = os.path.join(self.work_dir, 'circ_id.list')
        opts = {
            'list_id': circ_list
        }

        self.circrpm.set_options(opts)
        self.circrpm.on('end', self.run_detail)
        self.circrpm.run()

    def run_detail(self):
        self.circdetail = self.add_tool('whole_transcriptome.circrna.circdetail')
        circ_list = os.path.join(self.work_dir, 'circ_id.list')
        opts = {
            'list_id': circ_list
        }
        self.circdetail.set_options(opts)
        self.circdetail.on('end', self.run_getfasta)
        self.circdetail.run()

    def run_getfasta(self):
        self.getfasta = self.add_tool('whole_transcriptome.circrna.getfasta')
        opts = {
            'genome': self.option('genome').path,
            'details': self.circdetail.option('details'),
        }
        self.getfasta.set_options(opts)
        organism = ['Arabidopsis_thaliana', 'Camellia_sinensis', 'Glycine_max', 'Gossypium_arboreum',
                    'Gossypium_hirsutum', 'Homo_sapiens', 'Hordeum_vulgare', 'Mus_musculus', 'Nicotiana_benthamiana',
                    'Oryza_sativa', 'Pyrus_betulifolia', 'Solanum_lycopersicum', 'Solanum_tuberosum',
                    'Triticum_aestivum', 'Zea_mays', "Gossypium_raimondii", "Poncirus_trifoliata", 'Caenorhabditis_elegans', 'Latimeria_chalumnae']
        if self.option('organism_name') in organism:
            self.getfasta.on('end', self.run_ifdatabase)
        else:
            self.getfasta.on('end', self.run_refnew)
        self.getfasta.run()

    def run_refnew(self):
        self.refnew = self.add_tool('whole_transcriptome.circrna.refnew')
        opts = {
            'circfasta': self.getfasta.option('circfasta')
        }
        self.refnew.set_options(opts)
        self.refnew.on('end', self.set_output)
        self.refnew.run()

    def run_ifdatabase(self):
        self.ifdatabase = self.add_tool('whole_transcriptome.circrna.ifdatabase')
        opts = {
            'details': self.circdetail.option('details'),
            'organism_name': self.option('organism_name'),
            'circfasta': self.getfasta.option('circfasta'),
        }
        self.ifdatabase.set_options(opts)
        self.ifdatabase.on('end', self.run_cutfasta)
        self.ifdatabase.run()

    def run_cutfasta(self):
        self.cutfasta = self.add_tool('whole_transcriptome.circrna.cutfasta')
        opts = {
            'details': self.ifdatabase.option('detail_circbase'),
            'circfasta': self.getfasta.option('circfasta')
        }
        self.cutfasta.set_options(opts)
        self.cutfasta.on('end', self.set_output)
        self.cutfasta.run()

    def set_output(self):
        p = self.count.option('counts').path
        link_names1 = os.path.join(self.output_dir, os.path.basename(p))
        if os.path.exists(link_names1):
            os.remove(link_names1)
        shutil.copy(p, link_names1)

        # os.link(p, link_names1)
        self.count.option('counts').set_path(link_names1)

        q = self.circrpm.option('rpms').path
        link_names2 = os.path.join(self.output_dir, os.path.basename(q))
        if os.path.exists(link_names2):
            os.remove(link_names2)
        shutil.copy(q, link_names2)
        # os.link(q, link_names2)
        self.circrpm.option('rpms').set_path(link_names2)

        c = self.getfasta.option('circfasta').path
        link_names5 = os.path.join(self.output_dir, os.path.basename(c))
        if os.path.exists(link_names5):
            os.remove(link_names5)
        shutil.copy(c, link_names5)
        # os.link(q, link_names2)
        self.getfasta.option('circfasta').set_path(link_names5)

        organism = ['Arabidopsis_thaliana', 'Camellia_sinensis', 'Glycine_max', 'Gossypium_arboreum',
                    'Gossypium_hirsutum', 'Homo_sapiens', 'Hordeum_vulgare', 'Mus_musculus', 'Nicotiana_benthamiana',
                    'Oryza_sativa', 'Pyrus_betulifolia', 'Solanum_lycopersicum', 'Solanum_tuberosum',
                    'Triticum_aestivum', 'Zea_mays']
        if self.option('organism_name') in organism:
            d = self.ifdatabase.option('detail_circbase').path
            link_names3 = os.path.join(self.output_dir, os.path.basename(d))
            if os.path.exists(link_names3):
                os.remove(link_names3)
            shutil.copy(d, link_names3)
            # os.link(d, link_names3)
            self.ifdatabase.option('detail_circbase').set_path(link_names3)

            j = self.cutfasta.option('ref').path
            k = self.cutfasta.option('new').path
            link_names8 = os.path.join(self.output_dir, os.path.basename(j))
            if os.path.exists(link_names8):
                os.remove(link_names8)
            shutil.copy(j, link_names8)
            link_names9 = os.path.join(self.output_dir, os.path.basename(k))
            if os.path.exists(link_names9):
                os.remove(link_names9)
            shutil.copy(k, link_names9)
        else:
            w = self.circdetail.option('details').path
            link_names4 = os.path.join(self.output_dir, os.path.basename(w))
            if os.path.exists(link_names4):
                os.remove(link_names4)
            shutil.copy(w, link_names4)
            self.circdetail.option('details').set_path(link_names4)

            x = self.refnew.option('ref_no').path
            y = self.refnew.option('new_no').path
            link_names6 = os.path.join(self.output_dir, os.path.basename(x))
            if os.path.exists(link_names6):
                os.remove(link_names6)
            shutil.copy(x, link_names6)
            link_names7 = os.path.join(self.output_dir, os.path.basename(y))
            if os.path.exists(link_names7):
                os.remove(link_names7)
            shutil.copy(y, link_names7)

        self.end()

    def end(self):
        super(CircBrushModule, self).end()


def read_fastq_dir(fastq_dir):
    is_se = False
    sp2fqs_dict = dict()
    for line in open(os.path.join(fastq_dir, 'list.txt')):
        eles = line.strip().split('\t')
        fastq = os.path.join(fastq_dir, eles[0])
        sample, mate_type = eles[1:]
        if sample in sp2fqs_dict:
            if mate_type == 'l':
                sp2fqs_dict[sample].insert(0, fastq)
            elif mate_type == 'r':
                sp2fqs_dict[sample].append(fastq)
        else:
            sp2fqs_dict[sample] = [fastq]
    else:
        pe_sample_count = len(filter(lambda item: len(item[1]) > 1, sp2fqs_dict.items()))
        if pe_sample_count:
            if pe_sample_count == len(sp2fqs_dict):
                is_se = False
            else:
                raise Exception('mix mate type found in {} -> {}'.format(fastq_dir, sp2fqs_dict))
        else:
            is_se = True
        return is_se, sp2fqs_dict


class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test_ath(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'circ_brush_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.circ_brush',
            'instant': False,
            'options': {
                'circ_method': 'ciri2,find_circ',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/workspace/20191024/Longrna_workflow_9211_6466/FastpRna/output/fastq',
                'annotate': '/mnt/ilustre/users/sanger-dev/workspace/20191024/Longrna_workflow_9211_6466/LargeGush/output/filter_by_express/filtered_file/all_mrna.gtf',
                'genome': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/TAIR10_ensembl_44/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa',
                'organism_name': 'Arabidopsis_thaliana',
                'junction_reads': 2,
                'circrna_length': 100000
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'circ_brush_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.circ_brush',
            'instant': False,
            'options': {
                'circ_method': 'ciri2',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/workspace/20191018/Longrna_workflow_3076_7227/FastpRna/output/fastq',
                'annotate': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/gtf/Homo_sapiens.GRCh38.96.gtf',
                'genome': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa',
                'organism_name': 'Homo_sapiens',
                'junction_reads': 2,
                'circrna_length': 100000
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_ath')])
    unittest.TextTestRunner(verbosity=2).run(suite)

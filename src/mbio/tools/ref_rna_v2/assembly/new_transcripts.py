# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import glob
import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.ref_rna.trans_step import merged_add_code, gtf_to_trans2gene


class NewTranscriptsAgent(Agent):
    def __init__(self, parent):
        super(NewTranscriptsAgent, self).__init__(parent)
        options = [
            {"name": "tmap", "type": "infile", "format": "assembly.tmap"},
            {"name": "ref_fa", "type": "infile", "format": "ref_rna_v2.fasta"},
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},
            {"name": "merged_gtf", "type": "infile", "format": "gene_structure.gtf"},

            {"name": "change_id_gtf", "type": "outfile", "format": "ref_rna_v2.gtf"},
            {"name": "change_id_fa", "type": "outfile", "format": "ref_rna_v2.fasta"},
            {"name": "ref_and_new_gtf", "type": "outfile", "format": "ref_rna_v2.gtf"},
            {"name": "trans2gene", "type": "outfile", "format": "ref_rna_v2.common"},

            {'name': 'all_transcripts', 'type': 'outfile', 'format': 'ref_rna_v2.fasta'},
            {"name": "add_code_merged", "type": "outfile", "format": "ref_rna_v2.gtf"},
            {"name": "new_transcripts", "type": "outfile", "format": "ref_rna_v2.gtf"},
            {"name": "new_genes", "type": "outfile", "format": "ref_rna_v2.gtf"},
            {"name": "old_transcripts", "type": "outfile", "format": "ref_rna_v2.gtf"},
            {"name": "old_genes", "type": "outfile", "format": "ref_rna_v2.gtf"},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "8G"

    def end(self):
        super(NewTranscriptsAgent, self).end()


class NewTranscriptsTool(Tool):
    def __init__(self, config):
        super(NewTranscriptsTool, self).__init__(config)
        self.program = {
            'perl': 'miniconda2/bin/perl',
            'python': 'miniconda2/bin/python',
            'gffread': 'bioinfo/rna/cufflinks-2.2.1/gffread'
        }
        self.script = {
            'gtfmerge': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/scripts/gtfmerge.pl'),
            'assembly_stat': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/scripts/assembly_stat.py')
        }
        self.file = {
            'add_code_merged_gtf': os.path.join(self.work_dir, 'add_code_merged.gtf'),
            'change_id_merged_gtf': os.path.join(self.work_dir, 'change_id_merged.gtf'),
            'ref_and_new_gtf': os.path.join(self.work_dir, 'ref_and_new.gtf'),
            'trans2gene': os.path.join(self.output_dir, 'trans2gene'),
            'new_transcripts_fa': os.path.join(self.output_dir, 'new_transcripts.fa'),
            'all_transcripts_fa': os.path.join(self.output_dir, 'all_transcripts.fa'),
        }
        self.dir = {
            'merged': os.path.join(self.work_dir, 'merged'),
            'changed': os.path.join(self.work_dir, 'changed')
        }

    def run(self):
        super(NewTranscriptsTool, self).run()
        self.cal_merged_add_code()
        self.run_gtfmerge()
        self.run_assembly_stat_merged()
        self.run_assembly_stat_changed()
        self.concat_ref_and_new()
        self.cal_gtf_to_trans2gene()
        self.run_gffread_new()
        self.run_gffread_all()
        self.set_output()
        self.end()

    def cal_merged_add_code(self):
        merged_add_code(trans_file=self.option('merged_gtf').path, tmap_file=self.option('tmap').path,
                        new_trans=self.file['add_code_merged_gtf'])

    def run_gtfmerge(self):
        cmd = '{} {}'.format(self.program['perl'], self.script['gtfmerge'])
        cmd += ' -i {}'.format(self.file['add_code_merged_gtf'])
        cmd += ' -compare {}'.format(self.option('tmap').path)
        cmd += ' -ref {}'.format(self.option('ref_gtf').path)
        cmd += ' -o {}'.format(self.file['change_id_merged_gtf'])
        command = self.add_command("run_gtfmerge", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            return
        elif command.return_code == 1:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state("memory_limit", "memory is low!")
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in running {}".format(self.script['gtfmerge']))
                return
        else:
            self.set_error("fail to run {}".format(self.script['gtfmerge']))

    def run_assembly_stat_merged(self):
        if not os.path.isdir(self.dir['merged']):
            os.mkdir(self.dir['merged'])
        cmd = '{} {}'.format(self.program['python'], self.script['assembly_stat'])
        cmd += ' -tmapfile {}'.format(self.option('tmap').path)
        cmd += ' -transcript_file {}'.format(self.file['add_code_merged_gtf'])
        cmd += ' -out_new_trans {}'.format(os.path.join(self.dir['merged'], 'new_trans.gtf'))
        cmd += ' -out_new_genes {}'.format(os.path.join(self.dir['merged'], 'new_genes.gtf'))
        cmd += ' -out_old_trans {}'.format(os.path.join(self.dir['merged'], 'old_trans.gtf'))
        cmd += ' -out_old_genes {}'.format(os.path.join(self.dir['merged'], 'old_genes.gtf'))
        command = self.add_command("run_assembly_stat_merged", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            return
        elif command.return_code == 1:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state("memory_limit", "memory is low!")
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in running {}".format(self.script['assembly_stat']))
                return
        else:
            self.set_error("fail to run {}".format(self.script['assembly_stat']))

    def run_assembly_stat_changed(self):
        if not os.path.isdir(self.dir['changed']):
            os.mkdir(self.dir['changed'])
        cmd = '{} {}'.format(self.program['python'], self.script['assembly_stat'])
        cmd += ' -tmapfile {}'.format(self.option('tmap').path)
        cmd += ' -transcript_file {}'.format(self.file['change_id_merged_gtf'])
        cmd += ' -out_new_trans {}'.format(os.path.join(self.dir['changed'], 'new_trans.gtf'))
        cmd += ' -out_new_genes {}'.format(os.path.join(self.dir['changed'], 'new_genes.gtf'))
        cmd += ' -out_old_trans {}'.format(os.path.join(self.dir['changed'], 'old_trans.gtf'))
        cmd += ' -out_old_genes {}'.format(os.path.join(self.dir['changed'], 'old_genes.gtf'))
        command = self.add_command("run_assembly_stat_changed", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            return
        elif command.return_code == 1:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state("memory_limit", "memory is low!")
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in running {}".format(self.script['assembly_stat']))
                return
        else:
            self.set_error("fail to run {}".format(self.script['assembly_stat']))

    def concat_ref_and_new(self):
        with open(self.file['ref_and_new_gtf'], 'w') as handle:
            handle.write(open(self.option('ref_gtf').path).read())
            handle.write(open(os.path.join(self.dir['changed'], 'new_trans.gtf')).read())

    def cal_gtf_to_trans2gene(self):
        gtf_to_trans2gene(self.file['ref_and_new_gtf'], self.file['trans2gene'])

    def run_gffread_new(self):
        cmd = '{} {}'.format(self.program['gffread'], os.path.join(self.dir['changed'], 'new_trans.gtf'))
        cmd += ' -g {}'.format(self.option('ref_fa').path)
        cmd += ' -w {}'.format(self.file['new_transcripts_fa'])
        command = self.add_command("run_gffread_new", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            return
        elif command.return_code == 1:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state("memory_limit", "memory is low!")
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in running gffread")
                return
        else:
            self.set_error("fail to run gffread")

    def run_gffread_all(self):
        cmd = '{} {}'.format(self.program['gffread'], self.file['ref_and_new_gtf'])
        cmd += ' -g {}'.format(self.option('ref_fa').path)
        cmd += ' -w {}'.format(self.file['all_transcripts_fa'])
        command = self.add_command("run_gffread_all", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            return
        elif command.return_code == 1:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state("memory_limit", "memory is low!")
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in running gffread")
                return
        else:
            self.set_error("fail to run gffread")

    def set_output(self):
        for path in glob.glob(os.path.join(self.output_dir, '*.gtf')) + ['']:
            if os.path.isfile(path):
                os.remove(path)
        os.link(os.path.join(self.dir['changed'], 'new_trans.gtf'),
                os.path.join(self.output_dir, 'new_transcripts.gtf'))
        os.link(os.path.join(self.dir['changed'], 'new_genes.gtf'),
                os.path.join(self.output_dir, 'new_genes.gtf'))
        os.link(os.path.join(self.dir['merged'], 'old_trans.gtf'),
                os.path.join(self.output_dir, 'old_transcripts.gtf'))
        os.link(os.path.join(self.dir['merged'], 'old_genes.gtf'),
                os.path.join(self.output_dir, 'old_genes.gtf'))
        os.link(os.path.join(self.file['add_code_merged_gtf']),
                os.path.join(self.output_dir, 'add_code_merged.gtf'))
        os.link(os.path.join(self.file['change_id_merged_gtf']),
                os.path.join(self.output_dir, 'change_id_merged.gtf'))
        os.link(os.path.join(self.file['ref_and_new_gtf']),
                os.path.join(self.output_dir, 'ref_and_new.gtf'))

        self.option('add_code_merged').set_path(os.path.join(self.output_dir, 'add_code_merged.gtf'))
        self.option('all_transcripts').set_path(self.file['all_transcripts_fa'])
        self.option('new_transcripts').set_path(os.path.join(self.output_dir, 'new_transcripts.gtf'))
        self.option('new_genes').set_path(os.path.join(self.output_dir, 'new_genes.gtf'))
        self.option('old_transcripts').set_path(os.path.join(self.output_dir, 'old_transcripts.gtf'))
        self.option('old_genes').set_path(os.path.join(self.output_dir, 'old_genes.gtf'))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "new_transcripts_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "ref_rna_v2.assembly.new_transcripts",
            "instant": False,
            "options": dict(
                tmap="/mnt/ilustre/users/sanger-dev/workspace/20200324/Single_gffcompare_2222_4978/Gffcompare/output"
                     "/gffcmp.input.gtf.tmap",
                merged_gtf="/mnt/ilustre/users/sanger-dev/workspace/20200324/Single_stringtie_merge_6944_9677"
                           "/StringtieMerge/output/out.gtf",
                ref_gtf="/mnt/ilustre/users/sanger-dev/workspace/20200313/Refrna_tsg_37039/FileCheck"
                        "/Arabidopsis_thaliana.TAIR10.43.gtf",
                ref_fa="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana"
                       "/TAIR10_Ensembl_43/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)

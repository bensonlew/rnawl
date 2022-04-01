# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import glob
import os
import unittest
import shutil

from biocluster.workflow import Workflow


class CircBrushWorkflow(Workflow):
    '''
    last_modify: 2019.10.14
    '''

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CircBrushWorkflow, self).__init__(wsheet_object)
        CIRC_METHOD = ('ciri2,find_circ', 'ciri2', 'find_circ', 'circ_finder', 'circexplorer')
        options = [
            {'name': 'circ_method', 'type': 'string', 'default': CIRC_METHOD[0]},
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'annotate', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            {'name': 'task_id', 'type': 'string', 'default': ''},
            {'name': 'project_sn', 'type': 'string', 'default': ''},

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()


    def run(self):
        self.run_circ_brush()
        super(CircBrushWorkflow, self).run()


    def run_circ_brush(self):
        self.step.add_steps('circ_brush')
        self.circ_brush = self.add_module('whole_transcriptome.circ_brush')
        self.circ_brush.set_options({
            'annotate': self.option('annotate'),
            'genome': self.option('genome'),
            'fastq_dir': self.option('fastq_dir'),
            'circ_method': self.option('circ_method')
        })
        self.circ_brush.on('start', self.set_step, {'start': self.step.circ_brush})
        self.circ_brush.on('end', self.set_step, {'end': self.step.circ_brush})
        self.circ_brush.on('end', self.set_output)
        self.circ_brush.run()

    def set_output(self):
        p = self.circ_brush.option('details').path
        link_names1 = os.path.join(self.output_dir, os.path.basename(p))
        shutil.copy(p, link_names1)
        # os.link(p, link_names1)
        #self.circ_brush.option('details').set_path(link_names1)
        #
        # q = self.circrpm.option('rpms').path
        # link_names2 = os.path.join(self.output_dir, os.path.basename(q))
        # shutil.copy(q, link_names2)
        # # os.link(q, link_names2)
        # self.circrpm.option('rpms').set_path(link_names2)
        #
        # d = self.circdetail.option('details').path
        # link_names3 = os.path.join(self.output_dir, os.path.basename(d))
        # shutil.copy(d, link_names3)
        # # os.link(d, link_names3)
        # self.circdetail.option('details').set_path(link_names3)



    def set_output(self, event):
        for basename in os.listdir(self.circ_brush.output_dir):
            source = os.path.join(self.circ_brush.output_dir, basename)
            link_name = os.path.join(self.output_dir, basename)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.set_db()

    def set_db(self):
        # self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        self.database = self.api.api('whole_transcriptome.circrna')
        self.database.add_circrna(
            task_id=self.option('task_id'),
            project_sn=self.option('project_sn'),
            circ_detail=os.path.join(self.output_dir, 'detail.txt'),
            # tabular=os.path.join(self.output_dir, 'lncRNA_ortholog.tabular'),
            # main_id=self.option('main_id')
        )
        # self.logger.info('finish set_db at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        super(CircBrushWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.whole_transcriptome.report.circrna import CircBrushWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'circrna_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.report.circrna',
            'options': {
                'circ_method': 'ciri2',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/workspace/20191010/Longrna_workflow_2617_9288/FastpRna/output/fastq',
                'annotate': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/gtf/Homo_sapiens.GRCh38.96.gtf',
                'genome': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa',
                'task_id' : 'whole_transcriptome',
                'project_sn' : 'whole_transcriptome'
            }
        }
        wsheet = Sheet(data=data)
        wf = CircBrushWorkflow(wsheet)
        wf.sheet.id = 'whole_transcriptome'
        wf.sheet.project_sn = 'whole_transcriptome'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

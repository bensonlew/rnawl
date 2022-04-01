# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.agent import Agent
from mbio.packages.ref_rna_v2.functions import toolfuncdeco, runcmd
from biocluster.tool import Tool
import os
import shutil
import unittest

class DiffDoEnrichAgent(Agent):
    '''
    last_modify: 2020.08.25
    '''
    def __init__(self, parent):
        super(DiffDoEnrichAgent, self).__init__(parent)
        options = [
            {'name': 'diff_list', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'do_list', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'alpha', 'type': 'float', 'default': 0.05},
            {'name': 'pval', 'type': 'float', 'default': 0.05},
            {'name': 'method', 'type': 'string', 'default': 'bh'},
            {'name': 'result', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {"name": "do_version", "type": "string", "default": ""},
        ]
        self.add_option(options)
        self.step.add_steps('do_enrich')
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.do_enrich.start()
        self.step.update()

    def step_finish(self):
        self.step.do_enrich.finish()
        self.step.update()

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        else:
            method = self.option('method').lower()

    def set_resource(self):
        self._cpu = 1
        self._memory = '8G'

    @toolfuncdeco
    def end(self):
        super(DiffDoEnrichAgent, self).end()

class DiffDoEnrichTool(Tool):
    def __init__(self, config):
        super(DiffDoEnrichTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }

        obo_fp=self.config.SOFTWARE_DIR + "/database/Annotation/all/DO/version_20200813/doid.obo"
        self.script = {
            'do_linkage': os.path.join(
                self.config.PACKAGE_DIR, 'medical_transcriptome/geneset/do_linkage.py'
            ),
            'do_enrichment': os.path.join(
                self.config.PACKAGE_DIR, 'medical_transcriptome/geneset/do_enrichment.py'
            ),
            'do_annot': os.path.join(
                self.config.PACKAGE_DIR, 'medical_transcriptome/geneset/do_add_annot.py'
            ),

        }
        self.file = {
            'obo': obo_fp,
            'outfile': os.path.join(self.work_dir, 'outfile.txt'),
            'study': os.path.join(self.work_dir, 'diff.list'),
            'population': os.path.join(self.work_dir, 'all.list'),
            'association': os.path.join(self.work_dir, 'background.txt'),
            'result': os.path.join(self.output_dir, 'do_enrich_geneset_list_gene.xls')
        }

    @toolfuncdeco
    def run(self):
        super(DiffDoEnrichTool, self).run()
        self.run_do_linkage()
        self.run_find_enrichment()
        self.run_do_enrich_stats()
        self.set_output()
        self.end()

    @toolfuncdeco
    def run_do_linkage(self):

        cmd = '{} {}'.format(self.program['python'], self.script['do_linkage'])
        cmd += ' {}'.format(self.option('do_list').prop['path'])
        cmd += ' {}'.format("do_linkage.list")

        cmd_name = 'run_do_linkage'
        runcmd(self, cmd_name, cmd)


    @toolfuncdeco
    def run_find_enrichment(self):
        correct = self.option("method")
        correct_method = 3
        if correct.lower() == "bh":
            correct_method = 3
        elif correct.lower() == 'bonferroni':
            correct_method = 1
        elif correct.lower() == 'holm':
            correct_method = 2
        elif correct.lower() == 'by':
            correct_method = 4
        else:
            print('correct method not exist, BH will be used')


        cmd = '{} {}'.format(self.program['python'], self.script['do_enrichment'])

        cmd += ' {} {}'.format("-deg", self.option("diff_list").prop['path'])
        cmd += ' {} {}'.format("-g2p", "do_linkage.list")
        cmd += ' {} {}'.format("-correct", correct_method)
        cmd_name = 'run_do_enrichment'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def run_do_enrich_stats(self):
        cmd = '{} {}'.format(self.program['python'], self.script['do_annot'])
        cmd += ' {}'.format(self.option("diff_list").prop['path'] + ".do_enrichment.xls")
        cmd += ' {}'.format(self.output_dir + "/do_enrichment.xls")
        cmd_name = 'run_do_enrich_addannot'
        
        runcmd(self, cmd_name, cmd)


    @toolfuncdeco
    def set_output(self):
        self.option('result').set_path(self.output_dir + "/do_enrichment.xls")
        pass


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'do_enrich_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'medical_transcriptome.geneset.do_enrich',
            'instant': False,
            'options': {
                'diff_list': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/ref_annot_class/do/x.list',
                'do_list': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/ref_annot_class/do/id2terms.G.tsv',
                'method': 'BH'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

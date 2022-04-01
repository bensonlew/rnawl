# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.agent import Agent
from mbio.packages.ref_rna_v2.functions import toolfuncdeco, runcmd
from biocluster.tool import Tool
import os
import shutil
import unittest

class DiffReactomeEnrichAgent(Agent):
    '''
    last_modify: 2020.08.25
    '''
    def __init__(self, parent):
        super(DiffReactomeEnrichAgent, self).__init__(parent)
        options = [
            {'name': 'diff_list', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'reactome_list', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'alpha', 'type': 'float', 'default': 0.05},
            {'name': 'pval', 'type': 'float', 'default': 0.05},
            {'name': 'method', 'type': 'string', 'default': 'bonferroni'},
            {'name': 'result', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {"name": "reactome_version", "type": "string", "default": ""},
        ]
        self.add_option(options)
        self.step.add_steps('reactome_enrich')
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.reactome_enrich.start()
        self.step.update()

    def step_finish(self):
        self.step.reactome_enrich.finish()
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
        super(DiffReactomeEnrichAgent, self).end()

class DiffReactomeEnrichTool(Tool):
    def __init__(self, config):
        super(DiffReactomeEnrichTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }

        self.script = {
            'reactome_path': os.path.join(
                self.config.PACKAGE_DIR, 'medical_transcriptome/geneset/reactome_path.py'
            ),
            'reactome_enrichment': os.path.join(
                self.config.PACKAGE_DIR, 'medical_transcriptome/geneset/reactome_enrichment.py'
            ),
            'reactome_annot': os.path.join(
                self.config.PACKAGE_DIR, 'medical_transcriptome/geneset/reactome_add_annot.py'
            ),

        }
        self.file = {
            'outfile': os.path.join(self.work_dir, 'outfile.txt'),
            'study': os.path.join(self.work_dir, 'diff.list'),
            'population': os.path.join(self.work_dir, 'all.list'),
            'association': os.path.join(self.work_dir, 'background.txt'),
            'result': os.path.join(self.output_dir, 'reactome_enrich_geneset_list_gene.xls')
        }

    @toolfuncdeco
    def run(self):
        super(DiffReactomeEnrichTool, self).run()
        self.run_reactome_linkage()
        self.run_find_enrichment()
        # self.run_reactome_enrich_stats()
        self.set_output()
        self.end()

    @toolfuncdeco
    def run_reactome_linkage(self):

        cmd = '{} {}'.format(self.program['python'], self.script['reactome_path'])
        cmd += ' {}'.format(self.option('reactome_list').prop['path'])
        cmd += ' {}'.format("reactome_path.list")

        cmd_name = 'run_reactome_linkage'
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


        cmd = '{} {}'.format(self.program['python'], self.script['reactome_enrichment'])

        cmd += ' {} {}'.format("-deg", self.option("diff_list").prop['path'])
        cmd += ' {} {}'.format("-g2p", "reactome_path.list")
        cmd += ' {} {}'.format("-correct", correct_method)
        cmd_name = 'run_reactome_enrichment'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def run_reactome_enrich_stats(self):
        cmd = '{} {}'.format(self.program['python'], self.script['reactome_annot'])
        cmd += ' {}'.format(self.option("diff_list").prop['path'] + ".reactome_enrichment.xls")
        cmd += ' {}'.format(self.output_dir + "/reactome_enrichment.xls")
        cmd_name = 'run_reactome_enrich_addannot'

        runcmd(self, cmd_name, cmd)


    @toolfuncdeco
    def set_output(self):
        name = os.path.basename(self.option("diff_list").prop['path'])
        outfiles = [name + ".reactome_enrichment.xls"]
        for item in outfiles:
            linkfile = self.output_dir + '/' + item
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/' + item, linkfile)

        #self.option('result').set_path(self.file['result'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to reactome test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'reactome_enrich_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'medical_transcriptome.geneset.reactome_enrich',
            'instant': False,
            'options': {
                'diff_list': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/ref_annot_class2/reactome/reactome.geneset',
                'reactome_list': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/ref_annot_class2/reactome/reactome.list',
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

# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.workflow import Workflow
from function_enrichment import FunctionEnrichmentWorkflow
import json
import os


class CorepanEnrichmentWorkflow(FunctionEnrichmentWorkflow):
    def __init__(self, wsheet):
        super(CorepanEnrichmentWorkflow, self).__init__(wsheet)
        
        self.option('functype', 'KEGG')
        self.option('levelenriched', 'Pathway')
        self.option('graph', True)
        self.cpsets = self.add_tool('bac_comp_genome.comp_function')
        self.option('testset', self.cpsets.work_dir + '/corepan_sets.xls')

    def interface(self):
        opts = [
            {'name': 'corepan', 'type': 'infile', 'format': 'sequence.profile_table'},
            {'name': 'pancat', 'type': 'string'},
        ] 
        opts += super(CorepanEnrichmentWorkflow, self).interface()
        return opts

    def run_first(self):
        self.get_cpsets()

    def run(self):
        self.cpsets.on('end', self.run_enrich)
        super(CorepanEnrichmentWorkflow, self).run()

    def get_cpsets(self):
        options = {
                'corepan': self.option('corepan').path,
                'pancat': self.option('pancat'),
                'corepansets': 'Y'
                }
        self.cpsets.set_options(options)

        self.cpsets.run()


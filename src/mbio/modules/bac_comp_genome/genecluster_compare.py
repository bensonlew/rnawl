# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'

from biocluster.module import Module
import os, sys

class GeneclusterCompareModule(Module):
    def __init__(self, work_id):
        super(GeneclusterCompareModule, self).__init__(work_id)
        options = [
            {'name': 'gene_dir', 'type': 'string', },
            {'name': 'seq_dir', 'type': 'string', },
            {'name': 'region', 'type': 'string'},
            {'name': 'cu_name', 'type': 'string'},
            {'name': 'samples', 'type': 'string'},
        ]
        self.add_option(options)

        self.regions = self.add_tool('bac_comp_genome.region_compare')
        self.homo_pairs = self.add_tool('bac_comp_genome.homo_pairs')
        self.distance = self.add_tool('bac_comp_genome.region_distance')

    def check_options(self):
        pass

    def run(self):
        super(GeneclusterCompareModule, self).run()
        self.on_rely([self.regions, self.homo_pairs], self.run_distance)
        self.distance.on('end', self.set_output)
        self.run_regions()
        self.run_homo_pairs()

    def run_regions(self):
        options = {
            'seq_dir': self.option('seq_dir'),
            'samples': self.option('samples'),
            'region': self.option('region'),
        }
        self.regions.set_options(options)
        self.regions.run()

    def run_homo_pairs(self):
        options = {
            'gene_dir': self.option('gene_dir'),
            'samples': self.option('samples'),
        }
        self.homo_pairs.set_options(options)
        self.homo_pairs.run()

    def run_distance(self):
        options = {
            'regions': self.regions.work_dir + '/regions.xls',
            'gene_dir': self.option('gene_dir'),
            'samples': self.option('samples'),
            'pairs': self.homo_pairs.work_dir + '/homo_pairs.xls',
        }
        self.distance.set_options(options)
        self.distance.run()

    def set_output(self):
        regions_file = self.regions.work_dir + '/regions.xls'
        self.link(regions_file, self.option('cu_name') + '_regions.xls')
        homo_pairs_file = self.homo_pairs.work_dir + '/homo_pairs.xls'
        self.link(homo_pairs_file, self.option('cu_name') + '_homo_pairs.xls')
        distance_file = self.distance.work_dir + '/genecluster_compare.xls'
        self.link(distance_file, self.option('cu_name') + '_genecluster_compare.xls')
        self.end()

    def end(self):
        super(GeneclusterCompareModule, self).end()

# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
import pandas as pd
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class RegionDistanceAgent(Agent):
    def __init__(self, parent):
        super(RegionDistanceAgent, self).__init__(parent)
        options = [
            {'name': 'gene_dir', 'type': 'string'},
            {'name': 'samples', 'type': 'string'},
            {'name': 'regions', 'type': 'string'},
            {'name': 'pairs', 'type': 'string'}
        ]
        self.add_option(options)
        self.step.add_steps('region_dist')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.region_dist.start()
        self.step.update()

    def stepfinish(self):
        self.step.region_dist.finish()
        self.step.update()

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = len(self.option('samples').split(',')) + 1
        self._memory = str(self._cpu) + 'G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', '',''],
        ])
        super(RegionDistanceAgent, self).end()


class RegionDistanceTool(Tool):
    def __init__(self, config):
        super(RegionDistanceTool, self).__init__(config)

    def run(self):
        super(RegionDistanceTool, self).run()
        self.dist()
        self.set_output()
        self.end()

    def dist(self):
        regions = pd.read_csv(self.option('regions'), sep='\t', index_col=0)
        pairs = pd.read_csv(self.option('pairs'), sep='\t', header=None)
        gff = self.get_gff()
        dist = []
        for sp in regions.index:
            region = regions.loc[sp, :]
            df = pd.read_csv(gff[sp], sep='\t', index_col=0)
            df['Sequence id'] = map(lambda x: x.split('_ORF')[0], df['Sequence id'].values)
            t = df['Start'] > df['End']
            df.loc[t, ['Start', 'End']] = df.loc[t, ['End', 'Start']].values
            df['Sequence id'] = region[0]
            df['sample'] = sp
            df['distance'] = df[['Start', 'End', 'Sequence id']].apply(self.d, args=(region,), axis=1)
            col = ['sample', 'Sequence id', 'Start', 'End', 'Strand', 'distance']
            dist.append(df[col])
        dist = pd.concat(dist)
        dist.to_csv('dist.txt', sep='\t')
        in_regions = list(df[df['distance'] < 0].index)
        pairs = pairs[pairs[0].isin(in_regions) | pairs[1].isin(in_regions)]
        pairs.to_csv('pairs.txt', sep='\t')

        col = [
            'ref_geneid', 'ref', 'ref_location', 'ref_start', 'ref_end',
            'ref_strand', 'ref_distance',
            'query_geneid', 'query', 'query_location', 'query_start', 'query_end',
            'query_stand', 'query_distance', 'identity'
        ]
        col = '\t'.join(col) + '\n'
        out_put = pairs.apply(self.get_output, args=(dist,), axis=1)

        with open('genecluster_compare.xls', 'w') as w:
            w.write(col)
            map(w.write, out_put)

    def set_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for f in files:
                os.remove(os.path.join(root, f))
        pass

    def end(self):

        super(RegionDistanceTool, self).end()

    def get_output(self, row, dist):
        o = [row[1], ] + list(dist.loc[row[1], :]) +\
                [row[0], ] + list(dist.loc[row[0], :]) + [row[2], ]
        return '\t'.join(map(str, o)) + '\n'

    def d(self, row, region):
        if row[2] != region[0]:
            return 'NaN'
        if row[0] > region[2]:
            return row[0] - region[2]
        elif row[1] < region[1]:
            return region[1] - row[1]
        else:
            return -1

    def get_gff(self):
        sp_list = self.option('samples').split(',')
        gff = {}
        for f in sp_list:
            gff_dir = os.path.join(self.option('gene_dir'), f)
            gff[f] = gff_dir + '/' + f + '_CDS.gff'
        return gff

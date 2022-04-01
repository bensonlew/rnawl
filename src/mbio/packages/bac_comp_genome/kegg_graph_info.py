# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'

import argparse
import pandas as pd
import re
from table_format import *


def _main(args):
    graph = pd.read_csv(args.graph_info, sep='\t')
    graph = graph['KO'].str.split(',', expand=True).stack().reset_index(level=0).set_index('level_0').rename(columns={0:'KO'}).join(graph.drop('KO', axis=1))
    levels = ['Level3', 'Pathway', 'KO']
    annot = annot_format(args.annotable, 'KEGG', levels)
    annot['Pathway'] = map(lambda x: re.sub(r'ko', 'map', x), annot['Pathway'])

    comb = pd.merge(annot, graph, on=['Pathway', 'KO'])

    #comb = reshape(comb, 'sample', 'genelist')

    comb.to_csv('kegg_graph_info.xls', sep='\t', index=None)


if __name__ == '__main__':
    parse = argparse.ArgumentParser()
    parse.add_argument('-g', '--graph_info', required=True,
                       help='')
    parse.add_argument('-a', '--annotable',)

    args = parse.parse_args()

    _main(args)

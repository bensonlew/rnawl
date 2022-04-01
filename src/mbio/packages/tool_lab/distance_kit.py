# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging

import numpy as np
import pandas as pd
from cogent.maths import distance_transform

METHOD = {'abund_jaccard': 'dist_abund_jaccard',
          'binary_chisq': 'binary_dist_chisq',
          'binary_chord': 'binary_dist_chord',
          'binary_euclidean': 'binary_dist_euclidean',
          'binary_hamming': 'binary_dist_hamming',
          'binary_jaccard': 'binary_dist_jaccard',
          'binary_lennon': 'binary_dist_lennon',
          'binary_ochiai': 'binary_dist_ochiai',
          'binary_otu_gain': 'binary_dist_otu_gain',
          'binary_pearson': 'binary_dist_pearson',
          'binary_sorensen_dice': 'binary_dist_sorensen_dice',
          'bray_curtis': 'dist_bray_curtis',
          'bray_curtis_faith': 'dist_bray_curtis_faith',
          'bray_curtis_magurran': 'dist_bray_curtis_magurran',
          'canberra': 'dist_canberra',
          'chisq': 'dist_chisq',
          'chord': 'dist_chord',
          'euclidean': 'dist_euclidean',
          'gower': 'dist_gower',
          'hellinger': 'dist_hellinger',
          'kulczynski': 'dist_kulczynski',
          'manhattan': 'dist_manhattan',
          'morisita_horn': 'dist_morisita_horn',
          'pearson': 'dist_pearson',
          'soergel': 'dist_soergel',
          'spearman_approx': 'dist_spearman_approx',
          'specprof': 'dist_specprof'}


def main(args):
    func = getattr(distance_transform, METHOD[args.metrics])
    sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[args.sep]
    idf = pd.read_csv(args.input, sep=sep, index_col=0)
    datamtx = np.array(idf.T)
    dists = func(datamtx)
    odf = pd.DataFrame(dists, index=idf.columns, columns=idf.columns)
    odf.to_csv(args.output, sep='\t')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='To be described')
    parser.add_argument('--input', action='store', required=True,
                        help='text table file', metavar='<FILE>', dest='input')
    parser.add_argument('--sep', action='store', choices=['tab', 'comma', 'semicolon', 'space'], required=True,
                        help='separator type', metavar='<STR>', dest='sep')
    parser.add_argument('--metrics', action='store', choices=METHOD.keys(), required=True,
                        help='metric to use', dest='metrics')
    parser.add_argument('--output', action='store', required=True,
                        help='dist matrix file', metavar='<FILE>', dest='output')

    logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

    args = parser.parse_args()

    main(args)

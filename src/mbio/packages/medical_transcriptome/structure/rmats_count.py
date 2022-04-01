# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.packages.medical_transcriptome.functions import pkgsfuncdeco
import pickle
import pandas as pd

@pkgsfuncdeco
def main(args):
    for i, o in {args.jc: args.ojc, args.jcec: args.ojcec}.items():
        io(i, o)

@pkgsfuncdeco
def io(ifile, ofile):
    dct = dict()
    for pk in [line.strip() for line in open(ifile)]:
        for s, e2c in pickle.load(open(pk)).items():
            print pk, s
            if s in dct:
                for e in dct[s]:
                    dct[s][e].update(e2c[e])
            else:
                dct[s] = e2c
    else:
        df = pd.DataFrame(
            pd.Series(dict(map(lambda t: (t[0], len(t[1])), e2c.items())), name=s) for s, e2c in dct.items(),
        )
        df.index.name = 'SAMPLE'
        df['TOTAL'] = df.sum(axis=1)
        df.to_csv(ofile, sep='\t')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Script for counting event number at sample level')
    parser.add_argument('--jc', action='store', required=True, dest='jc',
                        help='JC pickle list file')
    parser.add_argument('--jcec', action='store', required=True, dest='jcec',
                        help='JCEC pickle list file')
    parser.add_argument('--ojc', action='store', required=True, dest='ojc',
                        help='JC result table file')
    parser.add_argument('--ojcec', action='store', required=True, dest='ojcec',
                        help='JCEC result table file')
    args = parser.parse_args()

    main(args)

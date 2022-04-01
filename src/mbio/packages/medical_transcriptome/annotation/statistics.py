# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Conduct statistics of annotation results')
parser.add_option('--input', dest='input', help='input relationship file among location, database and expression level')
parser.add_option('--map', dest='map', type=str, help='input T2G2R2L2P file')
parser.add_option('--output', dest='output', type=str, help='output statistics file')
(opts, args) = parser.parse_args()

def main(l2d2t_tsv, map_tsv, stat_tsv):
    t2g_dict = dict(line.strip().split('\t')[:2] for line in open(map_tsv))
    t2g_dict_unigene = dict(line.strip().split('\t')[:2] for line in open(map_tsv) if line.strip().split('\t')[2] == "yes")
    num_dict = {'transcript': len(set(t2g_dict.keys())), 'gene': len(set(t2g_dict_unigene.values()))}
    print 'INFO: start reading {}'.format(l2d2t_tsv)
    df = pd.read_table(l2d2t_tsv, header=None)
    dct = dict((db, {'type': db}) for db in set(df[1]))
    dct.update({'annotation': {'type': 'annotation'}})
    for name, group in df.groupby(2):
        group[0] = group[0].apply(lambda f: set(line.strip() for line in open(f) if line.strip()))
        ids_set = set()
        ids_set.update(*group[0])
        group = group.append({0: ids_set, 1: 'annotation', 2: name}, ignore_index=True)
        for index, row in group.iterrows():
            dct[row[1]].update({
                name: len(row[0]), '{}_percent'.format(name): '{:.4f}'.format(float(len(row[0]))/(num_dict[name] if num_dict[name] else 1))
            })
    else:
        dct.update({'total': {
            'type': 'total',
            'transcript': num_dict['transcript'],
            'transcript_percent': '1.000',
            'gene': num_dict['gene'],
            'gene_percent': '1.000',
        }})
    odf = pd.DataFrame(dct.values()).reindex([
        'type', 'transcript', 'gene', 'transcript_percent', 'gene_percent'
    ], axis=1)
    odf = odf.set_index('type')
    odf = odf.reindex(['go', 'kegg', 'cog', 'nr', 'uniprot', 'pfam', 'annotation', 'total'])
    odf = odf.rename({'go': 'GO', 'kegg': 'KEGG', 'cog': 'COG', 'nr': 'NR', 'uniprot': 'Uniprot', 'pfam': 'Pfam',
                      'annotation': 'Total_anno', 'total': 'Total'})
    odf.to_csv(stat_tsv, sep='\t')
    if os.path.getsize(stat_tsv) > 0:
        print 'INFO: succeed in exporting {}'.format(stat_tsv)

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 3, ['input', 'map', 'output'])):
        main(opts.input, opts.map, opts.output)
    else:
        parser.print_help()

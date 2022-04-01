# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
from Bio.Blast import NCBIXML
import re
import os

parser = OptionParser(description='Convert BLAST XML to specific table')
parser.add_option('-i', '--input', dest='input', help='input XML file')
parser.add_option('-a', '--header', dest='header', action='store_true', help='add header to file')
parser.add_option('-o', '--output', dest='output', help='output tabular file')
(opts, args) = parser.parse_args()

default_header = [
    'Score', 'E-Value', 'HSP-Len', 'Identity-%', 'Similarity-%',
    'Query-Name', 'Q-Len', 'Q-Begin', 'Q-End', 'Q-Frame',
    'Hit-Name', 'Hit-Len', 'Hsp-Begin', 'Hsp-End', 'Hsp-Frame', 'Hit-Description'
]

all_header = [
    'Score', 'E-Value', 'HSP-Len', 'Identity-%', 'Similarity-%', 'Identity', 'Positives',
    'Query-Name', 'Q-Len', 'Q-Begin', 'Q-End', 'Q-Frame',
    'Hit-Name', 'Hit-Len', 'Hsp-Begin', 'Hsp-End', 'Hsp-Frame', 'Hit-Description',
    'Q-Strand', 'Hsp-Strand', 'Mismatch', 'Gapopen_num'
]

def main(xml_in, tabular_out, columns=None, header=True):
    if columns:
        print 'INFO: define columns from specific input'
        for i in columns:
            if i not in all_header:
                raise Exception('ERROR: find unvalid item in header -> {}'.format(i))
    else:
        print 'INFO: define columns from default header'
        columns = default_header
    print 'INFO: {}'.format(columns)
    with open(tabular_out, 'w') as f:
        if header:
            print 'INFO: add header to {} according to input options'.format(tabular_out)
            f.write('{}\n'.format('\t'.join(columns)))
        print 'INFO: start parsing {} as generator object by {}'.format(xml_in, NCBIXML)
        records = NCBIXML.parse(open(xml_in))
        values = {i: 'N/A' for i in all_header}
        print 'INFO: start processing {}'.format(records)
        for record in records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    one_hsp = values.copy()
                    one_hsp['Query-Name'] = str(re.split(' ', record.query, maxsplit=1)[0])
                    one_hsp['Hit-Name'] = str(alignment.hit_id)
                    one_hsp['Hit-Description'] = str(alignment.hit_def)
                    one_hsp['Score'] = str(hsp.score)
                    one_hsp['E-Value'] = str(hsp.expect)
                    one_hsp['HSP-Len'] = str(hsp.align_length)
                    one_hsp['Identity'] = str(hsp.identities)
                    one_hsp['Positives'] = str(hsp.positives)
                    one_hsp['Q-Len'] = str(record.query_length)
                    one_hsp['Q-Begin'] = str(hsp.query_start)
                    one_hsp['Q-End'] = str(hsp.query_end)
                    one_hsp['Q-Frame'] = str(hsp.frame[0])
                    one_hsp['Hit-Len'] = str(alignment.length)
                    one_hsp['Hsp-Begin'] = str(hsp.sbjct_start)
                    one_hsp['Hsp-End'] = str(hsp.sbjct_end)
                    one_hsp['Hsp-Frame'] = str(hsp.frame[1])
                    one_hsp['Q-Strand'] = str(hsp.strand[0])
                    one_hsp['Hsp-Strand'] = str(hsp.strand[1])
                    one_hsp['Identity-%'] = str(round(float(hsp.identities) / hsp.align_length, 3) * 100)
                    one_hsp['Similarity-%'] = str(round(float(hsp.positives) / hsp.align_length, 3) * 100)
                    one_hsp['Mismatch'] = str(int(hsp.align_length) - int(len(hsp.match)))
                    one_hsp['Gapopen_num'] = str(hsp.gaps)
                    line = list()
                    for i in columns:
                        line.append(one_hsp[i])
                    f.write('{}\n'.format('\t'.join(line)))
    if os.path.getsize(tabular_out) > 0:
        print 'INFO: succeed in exporting {}'.format(tabular_out)

if __name__ == '__main__':
    if opts.input and opts.output:
        main(
            xml_in=opts.input,
            tabular_out=opts.output,
            columns=None,
            header=opts.header,
        )
    else:
        parser.print_help()
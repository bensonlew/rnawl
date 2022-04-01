# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import xml.etree.ElementTree as ET
import os

parser = OptionParser(description='Generate gene BLAST XML file')
parser.add_option('--t2g', dest='t2g', help='input T2G file')
parser.add_option('--xml', dest='xml', type=str, help='input transcript BLAST XML file')
parser.add_option('--output', dest='output', type=str, help='output gene BLAST XML file')
(opts, args) = parser.parse_args()

head = '''<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
'''

def main(t2g_tsv, txpt_xml, gene_xml):
    print 'INFO: start reading {}'.format(t2g_tsv)
    t2g_dict = dict(line.strip().split('\t')[:2] for line in open(t2g_tsv))
    print 'INFO: start reading {}'.format(txpt_xml)
    tree = ET.parse(txpt_xml)
    root = tree.getroot()
    BlastOutput_iterations = root.find('BlastOutput_iterations')
    for Iteration in BlastOutput_iterations.findall('Iteration'):
        query_def_obj = Iteration.find('Iteration_query-def')
        query_def = query_def_obj.text.split()[0]
        if query_def not in t2g_dict:
            BlastOutput_iterations.remove(Iteration)
        else:
            query_def_obj.text = t2g_dict[query_def]
    else:
        tree.write(gene_xml)
        print 'INFO: start adding file header to {}'.format(gene_xml)
        with open(gene_xml, 'r+') as f:
            body = f.read()
            f.seek(0, 0)
            f.write(head + body)
        if os.path.getsize(gene_xml) > 0:
            print 'INFO: succeed in exporting {}'.format(gene_xml)

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 3, ['t2g', 'xml', 'output'])):
        main(opts.t2g, opts.xml, opts.output)
    else:
        parser.print_help()

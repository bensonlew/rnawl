# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import os
import re

parser = OptionParser(description='Merge ref and new GTF file and export trans2gene file')
parser.add_option('-r', '--ref', dest='ref', help='Input ref GTF file')
parser.add_option('-n', '--new', dest='new', help='Input new GTF file')
parser.add_option('-a', '--all', dest='all', help='Output ref_and_new GTF file')
parser.add_option('-o', '--t2g', dest='t2g', help='Output trans2gene file')
(opts, args) = parser.parse_args()

def main(ref, new, ran, t2g):
    ret = merge(ref, new, ran)
    if ret:
        ret = extract(ran, t2g)
        if ref:
            print 'INFO: finish merging GTF file and exporting trans2gene file'
        else:
            raise Exception('ERROR: fail to create {}'.format(t2g))
    else:
        raise Exception('ERROR: fail to create {}'.format(ran))

def merge(ref, new, ran):
    print 'INFO: start merging {} and {}'.format(ref, new)
    with open(ref) as ref_handle, open(new) as new_handle:
        ref_lines = ref_handle.readlines()
        new_lines = new_handle.readlines()
    with open(ran, 'w') as ran_handle:
        ran_handle.writelines(ref_lines + new_lines)
    if os.path.isfile(ran) and os.path.getsize(ran) > 0:
        print 'INFO: succeed in creating {}'.format(ran)
        return True
    else:
        print 'WARNING: fail to create {}'.format(ran)
        return False

def extract(ran, t2g):
    print 'INFO: start extracting {} from {}'.format(t2g, ran)
    t2g_list = list()
    for n, line in enumerate(open(ran)):
        if line[0] != '#' and 'gene_id' in line and 'transcript_id' in line:
            transcript_id = re.search(r'transcript_id\s*"(\S+)";', line).group(1)
            gene_id = re.search(r'gene_id\s*"(\S+)";', line).group(1)
            t2g_list.append('{}\t{}\n'.format(transcript_id, gene_id))
        else:
            print 'WARNING: drop {} at line {}'.format(line, n + 1)
    with open(t2g, 'w') as t2g_handle:
        t2g_handle.writelines(sorted(list(set(t2g_list))))
    if os.path.isfile(t2g) and os.path.getsize(t2g) > 0:
        print 'INFO: succeed in creating {}'.format(t2g)
        return True
    else:
        print 'WARNING: fail to create {}'.format(t2g)
        return False

if __name__ == '__main__':
    if opts.ref and opts.new and opts.all and opts.t2g:
        main(opts.ref, opts.new, opts.all, opts.t2g)
    else:
        parser.print_help()
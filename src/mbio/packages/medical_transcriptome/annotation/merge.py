# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import numpy as np
import os

parser = OptionParser(description='Merge results from annotation classification module')
parser.add_option('--ref', dest='ref', type=str, help='input REF file')
parser.add_option('--new', dest='new', type=str, help='input NEW file')
parser.add_option('--header', dest='header', choices=['yes', 'no'], help='whether there are headers in files')
parser.add_option('--output', dest='output', type=str, help='output merged file')
(opts, args) = parser.parse_args()

def main(file_ref, file_new, header, file_out):
    print 'INFO: start reading {}'.format(file_ref)
    ref_lines = open(file_ref).readlines()
    if os.path.exists(file_new) and os.path.isfile(file_new):
        print 'INFO: start reading {}'.format(file_new)
        new_lines = open(file_new).readlines()
    else:
        new_lines = []

    if os.path.exists(os.path.dirname(file_out)):
        pass
    else:
        os.makedirs(os.path.dirname(file_out))
    open(file_out, 'w').writelines(merge(ref_lines, new_lines, True if header == 'yes' else False))
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

def merge(ref_lines, new_lines, header=True):
    if header:
        # arr = np.intersect1d(ref_lines[0], new_lines[0])
        if 1:
            if len(new_lines) >= 1:
                head = new_lines[0]
            else:
                head = ref_lines[0]
            lines = [head] + ref_lines[1:]
            if len(new_lines) >= 2:
                lines += new_lines[1:]
        else:
            raise Exception('ERROR: find inconsistent betwenn {} and {}'.format(ref_lines[0], new_lines[0]))
    else:
        lines = ref_lines + new_lines
    return lines

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 4, ['ref', 'new', 'header', 'output'])):
        main(opts.ref, opts.new, opts.header, opts.output)
    else:
        parser.print_help()

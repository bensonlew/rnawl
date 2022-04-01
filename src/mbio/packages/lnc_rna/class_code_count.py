# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, qinjincheng'

from optparse import OptionParser
import re
import os

parser = OptionParser(description='Export statistical results of transcripts count')
parser.add_option('-i', '--input', dest='input', help='Input GTF file of merged transcripts containing class code')
parser.add_option('-o', '--output', dest='output', help='Output class code statistical file')
(opts, args) = parser.parse_args()

def main(file_in, file_out):
    print 'INFO: start processing {}'.format(file_in)
    fr = open(file_in)
    txpt_cls_content = list()
    cls_content = set()
    for line in fr:
        m = re.match('#.*', line)
        if not m:
            nine_line = line.strip().split('\t')[-1]
            n = re.search(r'\s*transcript_id\s+\"(\S+)\";.*\s*class_code\s+\"(\S+)\";*', nine_line)
            if n:
                cls_content.add(n.group(2))
                new_line = 'transcript_id "{}";\t class_code "{}"\n'.format(n.group(1), n.group(2))
                txpt_cls_content.append(new_line)

    cls_txpt_set_dic = dict()
    for cls in cls_content:
        cls = cls.strip()
        if cls:
            cls_txpt_set_dic[cls] = {'txpt_set': set(), 'count': 0}

    for record in txpt_cls_content:
        m = re.search(r'\s*transcript_id\s+\"(\S+)\";\s*class_code\s+\"(\S+)\"', record.strip())
        if m:
            cls = m.group(2)
            txpt = m.group(1)
            cls_txpt_set_dic[cls]['txpt_set'].add(txpt)

    fw = open(file_out, 'wb')
    for cls in cls_txpt_set_dic.keys():
        cls_txpt_set_dic[cls]['count'] = len(cls_txpt_set_dic[cls]['txpt_set'])
        newline = '{}\t{}\t{}\n'.format(cls, ','.join(cls_txpt_set_dic[cls]['txpt_set']), str(cls_txpt_set_dic[cls]['count']))
        fw.write(newline)
    fw.close()
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.output:
        main(opts.input, opts.output)
    else:
        parser.print_help()
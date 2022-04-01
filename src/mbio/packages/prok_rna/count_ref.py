# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'

import sys


if len(sys.argv) != 3:
    exit('USAGE: %s ref.fa ref.fa.hdrs' % sys.argv[0])

ref = sys.argv[1]
hdrs = sys.argv[2]
with open(ref, 'r') as ref_r, open(hdrs, 'w') as hdrs_w:
    ref_info = ref_r.read().split('\n>')
    for block in ref_info:
        block = block.strip().lstrip('>').split('\n')
        id = block[0].split(' ')[0].strip()
        seq = ''.join(block[1:]).strip()
        hdrs_w.write('>' + id + ' ' + '/len=' + str(len(seq)) + ' ' + '/nonNlen=' + str(len(seq) - seq.lower().count('n')) + ' ' + '/org=ref' + '\n')

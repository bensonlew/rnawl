# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan, qinjincheng'

import sys
import os

def main(file_in, dir_out):
    print 'INFO: start reading {}'.format(file_in)
    lines = open(file_in).readlines()
    body = lines[1:]
    all_query = set()

    print 'INFO: start iterating and exporting files to {}'.format(dir_out)
    level2 = open(os.path.join(dir_out, 'go_level2.xls'), 'w')
    level3 = open(os.path.join(dir_out, 'go_level3.xls'), 'w')
    level4 = open(os.path.join(dir_out, 'go_level4.xls'), 'w')
    for i in range(1, 4):
        d = dict()
        for line in body:
            items = line.strip().split('\t')
            if len(items) < 2:
                continue
            term_type = items[0]
            go_term = items[i * 2 - 1]
            go_id = items[i * 2]
            seq_num = int(items[7])
            seq_list = items[8].split(';')
            if d.has_key(go_term):
                for seq in seq_list:
                    posi1 = seq.find('(')
                    query1 = seq[:posi1]
                    all_query.add(query1)
                    gos = seq[posi1 + 1: -1].split(',')
                    if d[go_term]['SeqList'].has_key(query1):
                        for go_item in gos:
                            d[go_term]['SeqList'][query1].add(go_item)
                    else:
                        go_set = set()
                        for go_item in gos:
                            go_set.add(go_item)
                        d[go_term]['SeqList'][query1] = go_set
            else:
                new_query_dict = dict()
                for seq in seq_list:
                    posi2 = seq.find('(')
                    query2 = seq[:posi2]
                    gos2 = seq[posi2 + 1: -1].split(',')
                    new_go_set = set()
                    all_query.add(query2)
                    for go_id2 in gos2:
                        new_go_set.add(go_id2)
                    new_query_dict[query2] = new_go_set
                d[go_term] = {'term_type': term_type, 'term': go_id, 'SeqList': new_query_dict}
        if i == 1:
            level2.write('term_type\tterm\tGO\tnumber\tpercent\tsequence\n')
        if i == 2:
            level3.write('term_type\tterm\tGO\tnumber\tpercent\tsequence\n')
        if i == 3:
            level4.write('term_type\tterm\tGO\tnumber\tpercent\tsequence\n')
        for k in d:
            s = d[k]['SeqList']
            pct = float(len(s)) / len(all_query)
            query_list = list()
            for t in s:
                query_list.append('{}({})'.format(t, ','.join(list(s[t]))))
            if i == 1:
                level2.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    d[k]['term_type'], k, d[k]['term'], len(s), str('%.8f' % pct), ';'.join(query_list)
                ))
            if i == 2:
                level3.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    d[k]['term_type'], k, d[k]['term'], len(s), str('%.8f' % pct), ';'.join(query_list)
                ))
            if i == 3:
                level4.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    d[k]['term_type'], k, d[k]['term'], len(s), str('%.8f' % pct), ';'.join(query_list)
                ))
    else:
        level2.close()
        level3.close()
        level4.close()
        print 'INFO: succeed in exporting 3 files to {}'.format(dir_out)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])

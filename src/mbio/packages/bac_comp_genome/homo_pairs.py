# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
import argparse
from multiprocessing import Pool


def diamond(arg):
    (diamond, ref, query) = arg
    db = os.path.splitext(os.path.basename(ref))[0]
    q = os.path.splitext(os.path.basename(query))[0]
    db_cmd = diamond + ' makedb --in {} -d diamond_db/{} -p 1'.format(ref, db)
    ret_code = os.system(db_cmd)
    if ret_code != 0:
        return db_cmd + ' 运行失败'

    align_cmd = diamond + ' blastp  -q {} -d diamond_db/{} -e 1e-10 --id 30 -o {}-{}.blast -p 2 -k 5'\
        .format(query, db, q, db)
    ret_code = os.system(align_cmd)
    if ret_code != 0:
        return align_cmd + ' 运行失败'


def get_pairs():
    blasts = filter(lambda x: x.endswith('blast'), os.listdir('.'))
    w = open('homo_pairs.xls', 'w')
    w.write('query_geneid\tref_geneid\tidentity\n')
    for b in blasts:
        with open(b, 'r') as r:
            for line in r.readlines():
                l = line.split('\t')
                w.write('{}\t{}\t{}\n'.format(l[0], l[1], l[2]))
    w.close()


def _main():
    if not os.path.exists('diamond_db'):
        os.mkdir('diamond_db')
    sp_list = args.samples.split(',')
    args_run = [(args.diamond, sp_list[i], sp_list[i+1]) for i in range(len(sp_list) - 1)]

    pools = Pool(len(sp_list))
    pret = pools.map_async(diamond, args_run)
    pools.close()
    pools.join()
    pret = pret.get()
    if any(pret):
        exit(pret)
    get_pairs()
    # diamond(args.diamond, args.ref, args.query)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--diamond',)
    parser.add_argument('-s', '--samples')

    args = parser.parse_args()

    _main()

# -*- coding: utf-8 -*-
# __author__: xieshichang
import subprocess
from biocluster.config import Config
import time


client = Config().get_mongo_client(mtype="metagenomic", ref=True)
db = client[Config().get_mongo_dbname("metagenomic", ref=True)]
nr_collection = db.NR_sequence
taxon_collection = db.species_taxon


def get_taxid(gis, taxid):
    t1 = time.time()
    col_result = nr_collection.find({'_id': {'$in': gis}}, {'taxid': 1})
    for n in col_result:
        if ('taxid' in n) and n['taxid']:
            taxid[str(n['_id'])] = int(n['taxid'])
    print("get taxid time used one run: {}, len: {}".format(time.time() - t1, len(taxid)))


def get_taxinfo(taxid_list, tax_infos):
    t1 = time.time()
    taxons = taxon_collection.find({'_id': {'$in': taxid_list}})
    for tax_info in taxons:
        tax_id = tax_info['_id']
        detail = tax_info['detail']
        for i in range(len(detail)):
            detail[i] = '{2}{{{0}}}'.format(*detail[i].partition(':'))
        tax_infos[tax_id] = [str(tax_id),';'.join(detail)]
    print("get taxinfo time used one run: {}, len: {}".format(time.time() - t1, len(tax_infos)))


def gi_taxon(gilist):
    if not isinstance(gilist, list) and not isinstance(gilist, set):
        raise TypeError('传入参数必须是gi的列表')
    gilist = map(int, gilist)
    size = 5000
    gilist_list = [gilist[i:i+size] for i in range(0, len(gilist), size)]

    gi2taxid = {}
    tax_infos = {}

    t1 = time.time()
    [get_taxid(n_list, gi2taxid) for n_list in gilist_list]
    print("get gi2taxid time used: {}, len: {}".format(time.time() - t1, len(gi2taxid)))

    taxids = list(set(gi2taxid.values()))
    ids_list = [taxids[i:i+size] for i in range(0, len(taxids), size)]

    t1 = time.time()
    [get_taxinfo(ids, tax_infos) for ids in ids_list]
    print("get tax_infos time used: {}, len: {}".format(time.time() - t1, len(tax_infos)))

    taxons = {name: tax_infos[tax_id] for name, tax_id in gi2taxid.items() if tax_id in tax_infos}
    client.close()
    return taxons


if __name__ == '__main__':  # for test
    import sys
    gilist = []
    with open(sys.argv[1], 'r') as r, open('o_new2.txt', 'w') as w:
        [gilist.append(l.strip('\n')) for l in r]

        rets = gi_taxon(gilist)
        for k in rets:
            w.write(('{}\t{}\n'.format(k, rets[k])))

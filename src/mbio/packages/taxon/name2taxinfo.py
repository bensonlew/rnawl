# -*- coding: utf-8 -*-
# __author__: xieshichang
import subprocess
from biocluster.config import Config
import time


client = Config().get_mongo_client(mtype="metagenomic", ref=True)
db = client[Config().get_mongo_dbname("metagenomic", ref=True)]
nr_collection = db.NR_sequence_20200604
taxon_collection = db.species_taxon_20200604


def get_taxid(namelist, taxid):
    t1 = time.time()
    col_result = nr_collection.find({'origin_sequence_id': {'$in': namelist}}, {'origin_sequence_id':1, 'taxid':1})
    for n in col_result:
        if 'taxid' in n:
            taxid[n['origin_sequence_id']] = int(n['taxid'])
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


def name_taxon(namelist):
    if not isinstance(namelist, list) and not isinstance(namelist, set):
        raise TypeError('传入参数必须是gi的列表')
    namelist = list(namelist)
    size = 5000
    namelist_list = [namelist[i:i+size] for i in range(0, len(namelist), size)]

    name2taxid = {}
    tax_infos = {}

    t1 = time.time()
    [get_taxid(n_list, name2taxid) for n_list in namelist_list]
    print("get name2taxid time used: {}, len: {}".format(time.time() - t1, len(name2taxid)))

    taxids = list(set(name2taxid.values()))
    ids_list = [taxids[i:i+size] for i in range(0, len(taxids), size)]

    t1 = time.time()
    [get_taxinfo(ids, tax_infos) for ids in ids_list]
    print("get tax_infos time used: {}, len: {}".format(time.time() - t1, len(tax_infos)))

    taxons = {name: tax_infos[tax_id] for name, tax_id in name2taxid.items() if tax_id in tax_infos}
    client.close()
    return taxons


if __name__ == '__main__':  # for test
    import sys
    gilist = []
    with open(sys.argv[1], 'r') as r, open('o_new2.txt', 'w') as w:
        [gilist.append(l.strip('\n')) for l in r]

        rets = name_taxon(gilist)
        for k in rets:
            w.write(('{}\t{}\n'.format(k, rets[k])))

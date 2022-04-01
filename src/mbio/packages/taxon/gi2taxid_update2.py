# -*- coding: utf-8 -*-
import subprocess
from biocluster.config  import Config
## __author__: xieshichang

db = Config().get_mongo_client(mtype="metagenomic", ref=True)[Config().get_mongo_dbname("metagenomic", ref=True)]
nr_collection = db.NR_sequence
taxon_collection = db.species_taxon


def __get_taxinfo(gis, taxons):
    gi_infos = nr_collection.find({'_id': {'$in': gis}})
    if gi_infos:
        for gi_info in gi_infos:
            if "taxid" in gi_info:
                tax_id = int(gi_info['taxid'])
            else:
                continue
            tax_info = taxon_collection.find_one({'_id': tax_id})
            if not tax_info:
                print 'WARNNING:taxid:{}在表中不存在'.format(tax_id)
            else:
                detail = tax_info['detail']
                for i in range(len(detail)):
                    detail[i] = '{2}{{{0}}}'.format(*detail[i].partition(':'))
                taxons[str(gi_info['_id'])] = ';'.join(detail)


def gi_taxon(gilist):
    if not isinstance(gilist, list) and not isinstance(gilist, set):
        raise TypeError('传入参数必须是gi的列表')
    gilist = map(str, gilist)
    size = 5000
    gilist_list = [gilist[i:i+size] for i in range(0, len(gilist), size)]
    taxons = {}
    [__get_taxinfo(gis, taxons) for gis in gilist_list]
    return taxons

if __name__ == '__main__':  # for test
    import sys
    gilist = []
    with open(sys.argv[1], 'r') as r, open('o_new.txt', 'w') as w:
        [gilist.append(l.strip('\n')) for l in r]

        rets = gi_taxon(gilist)
        for k in rets:
            w.write(('{}\t{}\n'.format(k, rets[k])))


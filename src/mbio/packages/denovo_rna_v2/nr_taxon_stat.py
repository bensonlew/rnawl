# -*- coding: utf-8 -*-

from mbio.packages.taxon.gi2taxon import gi_taxon
from mbio.packages.rna.accession2taxon import taxon
import sys

def filter_query(fp):
    query_acc = dict()
    openfp = open(fp)
    openfp.readline()
    for i in openfp:
        line_sp = i.split('\t')
        if line_sp[5] in query_acc:
            continue
        else:
            hitname = line_sp[10].split('|')
            if line_sp[10].startswith("gi|"):
                query_acc[line_sp[5]] = hitname[3]
            else:
                query_acc[line_sp[5]] = hitname[0]
    return query_acc

if __name__ == '__main__':
    
    fp = sys.argv[1]
    out = sys.argv[2]
    db = sys.argv[3]
    acc2tax_db = sys.argv[4]
    taxon = taxon(db=db, acc2tax=acc2tax_db)
    query_acc_all = filter_query(fp)
    print(query_acc_all.items()[:5])
    query_acc = list(set(query_acc_all.values()))
    acc_taxons = taxon.get_accession_ranks(query_acc, ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
    acc2name = dict(zip(query_acc, acc_taxons))

    with open(out, 'w') as w:
        for item in query_acc_all.iteritems():
            acc = query_acc_all[item[0]]
            species_names = acc2name[acc]
            w.write(item[0] + '\t' + acc  + '\t' + species_names + '\n')

    '''
        # 生成物种分类树 tree
        cmd = "{} {} {}".format(self.python, self.taxon_tree, self.output_dir + '/query_taxons_detail.xls')
        command = self.add_command("taxon_tree", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("生成物种分类树")
        else:
            self.set_error("生成分类树出错")

        return True
    '''



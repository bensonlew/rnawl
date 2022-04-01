# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
# last_modify: 20170606

import sys
# from pymongo import MongoClient
from biocluster.config import Config
from Bio.Blast import NCBIXML


class meta_kegg_anno(object):
    """
    meta宏基因的kegg注释需要的信息和已有的kegg注释结果不符，所以重写了这个部分
    """
    def __init__(self):
        # self.client = MongoClient('mongodb://10.100.200.129:27017')
        # self.kegg_gene = self.client.sanger_biodb.kegg_gene
        # self.kegg_ko = self.client.sanger_biodb.kegg_ko
        self.client = Config().get_mongo_client(mtype="metagenomic")
        self.db = client[Config().get_mongo_dbname("metagenomic")]
        self.kegg_gene = self.db.kegg_gene
        self.kegg_ko = self.db.kegg_ko

    def kegg_by_mongo(self, kegg_table, out_dir):
        kegg_anno = out_dir + '/kegg_anno.xls'
        kegg_enzyme = out_dir + '/kegg_enzyme_list.xls'
        kegg_module = out_dir + '/kegg_module_list.xls'
        kegg_pathway = out_dir + '/kegg_pathway_list.xls'
        with open(kegg_table, 'r') as table, open(kegg_anno, 'wb') as anno, open(kegg_enzyme, 'wb') as e, \
                open(kegg_module, 'wb') as m, open(kegg_pathway, 'wb') as p:
            anno.write('#Query\tGene\tKO\tDefinition\tPathway\tEnzyme\tModules\tHyperlink\n')
            query_list = []
            for line in table:
                line = line.strip("\n").split("\t")
                if line[0] != "Score":
                    query_name = ('_').join(line[5].split("_")[:-1])  # Query
                    if query_name in query_list:
                        pass
                    else:
                        query_list.append(query_name)
                        koids = self.kegg_gene.find_one({"gene_id": line[10]})  # Gene = line[10]
                        if koids:
                            ko_id = koids["koid"]                               # KO
                            result = self.kegg_ko.find_one({"ko_id": ko_id})
                            if result:
                                des = result["ko_desc"]  # Gene = line[10]
                                path = result["pathway_id"]  # Pathway 的列表形式
                                path_ = (',').join(path)
                                if len(path) == 0:
                                    path_ = "-"
                                else:
                                    pa_list = []
                                    path_des = result["pathway_category"]
                                    for i in path_des:
                                        pa_list.append(i[-1])
                                    pa_des = (';').join(pa_list)
                                    p.write('{}\t{}\t{}\n'.format(query_name, path_, pa_des))
                                _enzyme = result["enzyme_id"]  # Enzyme
                                enzyme = (',').join(_enzyme)
                                if len(result["enzyme_id"]) == 0:
                                    enzyme = "-"
                                else:
                                    en_list = []
                                    enzyme_des = result["enzyme_category"]
                                    for i in enzyme_des:
                                        en_list.append(i[-1])
                                    en_des = (';').join(en_list)
                                    e.write('{}\t{}\t{}\n'.format(query_name, enzyme, en_des))
                                module = result["module_id"]
                                module_ = (',').join(module)  # Module
                                if len(module) == 0:
                                    module_ = "-"
                                else:
                                    mo_list = []
                                    module_des = result["module_category"]
                                    for i in module_des:
                                        mo_list.append(i[-1])
                                    mo_des = (';').join(mo_list)
                                    m.write('{}\t{}\t{}\n'.format(query_name, module_, mo_des))
                            else:
                                des = '-'
                                path_ = '-'
                                enzyme = '-'
                                module_ = '-'
                                # enzyme_des = '-'
                            h = 'http://www.genome.jp/dbget-bin/www_bget?ko:' + ko_id
                            anno.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.
                                       format(query_name, line[10], ko_id, des, path_, enzyme, module_, h))
                        else:
                            continue
                else:
                    continue

if __name__ == '__main__':
    meta_kegg_anno = meta_kegg_anno()
    meta_kegg_anno.kegg_by_mongo(kegg_table=sys.argv[1], out_dir=sys.argv[2])
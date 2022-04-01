# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modify: 20180226

import sys
from biocluster.config import Config
from pymongo import MongoClient
from Bio.Blast import NCBIXML
from bson.objectid import ObjectId
import re
import logging

class kegg_xlm_mongo(object):
    """
    利用KEGG注释的xml从mongo库获取层级信息
    """
    def __init__(self):
        #self.client = MongoClient('mongodb://10.100.200.129:27017')
        #self.mongo = self.client.sanger_biodb
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.kegg_gene = self.mongodb['kegg_v94.2_gene']
        self.kegg_ko = self.mongodb['kegg_v94.2_ko']
        #self.kegg_png = self.client.sanger_biodb.kegg_pathway_png_v1

    def kegg_by_mongo(self, kegg_table, out_dir):
        kegg_anno = out_dir + '/gene_kegg_anno.xls'
        kegg_enzyme = out_dir + '/kegg_enzyme_list.xls'
        kegg_module = out_dir + '/kegg_module_list.xls'
        kegg_pathway = out_dir + '/kegg_pathway_list.xls'
        kegg_level = out_dir + '/kegg_level.xls'
        #fs = gridfs.GridFS(self.mongodb)
        #f = open(pidpath)
        with open(kegg_table, 'r') as table, open(kegg_anno, 'wb') as anno, open(kegg_enzyme, 'wb') as e, \
            open(kegg_module, 'wb') as m, open(kegg_pathway, 'wb') as p, open(kegg_level, 'wb') as outfile:
            anno.write('#Query\tGene\tKO\tDefinition\tPathway\tEnzyme\tModule\tHyperlink\tIdentity(%)\tAlign_len\tLevel1\tLevel2\tLevel3\n')
            outfile.write('#Query\tlevel1\tlevel2\tlevel3\tpathway_id\tKO\n')
            query_list = []
            for line in table:
                line = line.strip("\n").split("\t")
                align_len = line[2]
                iden = line[3]
                if line[0] != "Score":
                    #query_name = ('_').join(line[5].split("_")[:-1])  # Query
                    if line[5].endswith("_1"):
                        query_name = line[5].rsplit("_1",1)[0]
                    else:
                        query_name = line[5]
                    if query_name in query_list:
                        pass
                    else:
                        query_list.append(query_name)
                        koids = self.kegg_gene.find_one({"gene_id": line[10]})  # Gene = line[10]
                        if koids:
                            ko_id = koids["koid"]                               # KO
                            result = self.kegg_ko.find_one({"ko_id": ko_id})
                            if result:
                                ko_name = result["ko_name"]
                                des = result["ko_desc"]  # Gene = line[10]
                                path = result["pathway_id"]  # Pathway 的列表形式
                                path_ = ';'.join(path)
                                if len(path) == 0:
                                    path_ = "-"
                                else:
                                    pa_list = []
                                    path_des = result["pathway_category"]
                                    for i in path_des:
                                        pa_list.append(i[-1])
                                    pa_des = ';'.join(pa_list)
                                    if path_ != "-":
                                        p.write('{}\t{}\t{}\n'.format(query_name, path_, pa_des))
                                _enzyme = result["enzyme_id"]  # Enzyme
                                enzyme = ';'.join(_enzyme)
                                if len(result["enzyme_id"]) == 0:
                                    enzyme = "-"
                                else:
                                    en_list = []
                                    enzyme_des = result["enzyme_category"]
                                    for i in enzyme_des:
                                        en_list.append(i[-1])
                                    logging.info(ko_id)
                                    en_des = ';'.join(en_list)
                                    e.write('{}\t{}\t{}\n'.format(query_name, enzyme, en_des))
                                module = result["module_id"]
                                module_ = ';'.join(module)  # Module
                                if len(module) == 0:
                                    module_ = "-"
                                else:
                                    mo_list = []
                                    module_des = result["module_category"]
                                    for i in module_des:
                                        mo_list.append(i[-1])
                                    mo_des = ';'.join(mo_list)
                                    #mo_des = ';'.join(module_des)
                                    m.write('{}\t{}\t{}\n'.format(query_name, module_, mo_des))
                                levels_paths = result["pathway_category"]
                                level1 = []
                                level2 = []
                                level3 = []
                                i = 0
                                if len(levels_paths) != 0:
                                    for each in levels_paths:
                                        levels = "\t".join(each)
                                        level1.append(each[0])
                                        level2.append(each[1])
                                        level3.append(each[2])
                                        #pathway_id = path_.split(";")[i]
                                        pathway_id = result["pathway_id"][i]
                                        i = i + 1
                                        outfile.write(query_name + "\t" + levels + "\t" + pathway_id + "\t" + ko_id + "\n")
                                    l1 = ";".join(level1)
                                    l2 = ";".join(level2)
                                    l3 = ";".join(level3)
                                else:
                                    l1 = "-"
                                    l2 = "-"
                                    l3 = "-"
                                h = 'http://www.genome.jp/dbget-bin/www_bget?ko:' + ko_id
                                anno.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query_name, ko_name, ko_id, des, path_, enzyme, module_, h, iden, align_len, l1, l2, l3))
                            else:
                                continue
                                #des = '-'
                                #path_ = '-'
                                #enzyme = '-'
                                #module_ = '-'
                                # enzyme_des = '-'
                            #h = 'http://www.genome.jp/dbget-bin/www_bget?ko:' + ko_id
                            #anno.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.
                            #           format(query_name, ko_name, ko_id, des, path_, enzyme, module_, h, iden, align_len, l1, l2, l3))
                        else:
                            continue
                else:
                    continue

if __name__ == '__main__':
    meta_kegg_anno = kegg_xlm_mongo()
    meta_kegg_anno.kegg_by_mongo(kegg_table=sys.argv[1], out_dir=sys.argv[2])

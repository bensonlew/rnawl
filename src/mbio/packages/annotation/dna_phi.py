# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'


import sys
from pymongo import MongoClient
import argparse
import os
from biocluster.config import Config
from pymongo.errors import ServerSelectionTimeoutError
from pymongo.errors import NetworkTimeout
import time

class phi_anno(object):
    """
    dna phi数据详细注释信息
    ardb.parse.anno.xls
    """

    def __init__(self):
        # self.client = MongoClient('mongodb://10.100.200.129:27017')
        # self.mongodb = self.client.sanger_biodb
        # self.mongodb = Config().biodb_mongo_client.sanger_biodb
        self.process_rerun = 0
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.phi = self.mongodb.PHI_v49

    def phi_by_mongo(self, align_table, anno_table):
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile, open('phi_stat', 'wb') as outfile1:
            outfile.write(
                'Gene ID\tPHI ID\tProtein ID\tHit Gene Name\tNCBI Tax ID\tPathogen Species\tPhenotype\tHost Descripton\tHost NCBI ID\tHost Species\tGene Function\tIdentity(%)\tEvalue\tDisease\tCoverage(%)\tScore\n')
            # query_list = []
            # done_list = []
            phi = {}
            des = {}
            type = {}
            for line in infile:
                line = line.strip().split("\t")
                if not "Score" in line[0]:
                    query = line[5]
                    hit = line[10]
                    coverge = round(abs(float(line[8])-float(line[7]))/float(line[6]),3)
                    coverge = coverge*100
                    id_list = hit.split("#")[1].split("__")
                    evalue = line[1]
                    score = line[0]
                    align_len = line[2]
                    act_iden = line[3]
                    pro_list = []
                    hitgene_list = []
                    ncbi_tax_list = []
                    path_spe_list = []
                    phenotype_list = []
                    host_des_list = []
                    host_id_list = []
                    exp_list = []
                    function_list = []
                    disease_list = []
                    for id in id_list:
                        try:
                            detail = self.phi.find_one({"PHI": id})
                        except (ServerSelectionTimeoutError, NetworkTimeout):  # 捕获因为mongo服务器问题导致的异常后重运行此方法
                            if self.process_rerun < 5:
                                self.process_rerun += 1
                                # self.logger.info("检测到TimeoutError, 第{}次重运行方法".format(self.process_rerun))
                                time.sleep(5)
                                self.phi_by_mongo(align_table, anno_table)
                            else:
                                raise Exception("重运行5次仍未成功连接mongo")
                        if detail:
                            pro = detail["ProteinID"]
                            hitgene = detail["Gene_name"]
                            ncbi_tax = detail["Pathogen_NCBI_species_Taxonomy_ID"]
                            path_spe = detail["Pathogen_species"]
                            phenotype = detail["Phenotype_of_mutant"]
                            host_des = detail["Host_description"]
                            host_id = detail["Host_NCBI_Taxonomy_ID"]
                            exp = detail["Experimental_host_species"]
                            function = detail["Function"]
                            disease = detail["disease"]
                            if id in phi:
                                phi[id] += 1
                            else:
                                phi[id] = 1
                                des[id] = phenotype
                            pro_list = list(set(pro_list).union(set(pro.split('; '))))
                            hitgene_list = list(set(hitgene_list).union(set(hitgene.split('; '))))
                            ncbi_tax_list = list(set(ncbi_tax_list).union(set(ncbi_tax.split('; '))))
                            path_spe_list = list(set(path_spe_list).union(set(path_spe.split('; '))))
                            phenotype_list = list(set(phenotype_list).union(set(phenotype.split('; '))))
                            host_des_list = list(set(host_des_list).union(set(host_des.split('; '))))
                            host_id_list = list(set(host_id_list).union(set(host_id.split('; '))))
                            exp_list = list(set(exp_list).union(set(exp.split('; '))))
                            function_list = list(set(function_list).union(set(function.split('; '))))
                            disease_list = list(set(disease_list).union(set(disease.split("; "))))
                        else:
                            print(id + "wrong ID")
                    for one_type in phenotype_list:
                        if one_type in type:
                            type[one_type] += 1
                        else:
                            type[one_type] = 1
                    outfile.write(
                        '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, ";".join(id_list),
                                                                                      ";".join(pro_list),
                                                                                      ";".join(hitgene_list),
                                                                                      ";".join(ncbi_tax_list),
                                                                                      ";".join(path_spe_list),
                                                                                      ";".join(phenotype_list),
                                                                                      ";".join(host_des_list),
                                                                                      ";".join(host_id_list),
                                                                                      ";".join(exp_list),
                                                                                      ";".join(function_list),
                                                                                       act_iden, evalue,
                                                                                      ";".join(disease_list),
                                                                                       coverge, score))
        phi_file = os.path.dirname(anno_table) + '/phi_stat.xls'
        with open(phi_file, 'wb') as outfile1:
            outfile1.write('PHI ID\tPHI Description\tPHI No.\n')
            for key in phi:
                outfile1.write('{}\t{}\t{}\n'.format(key, des[key], phi[key]))
        phenotype_file = os.path.dirname(anno_table) + '/phi_phenotype.xls'
        with open(phenotype_file, 'wb') as outfile2:
            outfile2.write('phenotype\tGene No.\n')
            for key in type:
                outfile2.write('{}\t{}\n'.format(key, type[key]))
                #
                # outfile1.write()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[xml_table]', required=True, help='Input xml table')
    parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    dna_phi_anno = phi_anno()
    dna_phi_anno.phi_by_mongo(align_table, anno_table)

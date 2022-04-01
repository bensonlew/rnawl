# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'


import sys
from pymongo import MongoClient
import argparse
import os
from biocluster.config import Config


class phi_anno(object):
    """
    meta phi数据详细注释信息
    ardb.parse.anno.xls
    """

    def __init__(self):
        # self.client = MongoClient('mongodb://10.100.200.129:27017')
        # self.mongodb = self.client.sanger_biodb
        # self.mongodb = Config().biodb_mongo_client.sanger_biodb
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.phi = self.mongodb.PHI_v49 #by zhaozhigang 20200929

    def phi_by_mongo(self, align_table, anno_table):
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:  # , open('phi_stat', 'wb') as outfile1:
            outfile.write(
                '#Query\tPHI ID\tprotein\tHit Gene Name\tNCBI Tax ID\tPathogen Species\tPhenotype\tHost Description\tHost NCBI ID\tExperimental Host Species\tFunction\tDisease\tIdentity(%)\tAlign_len\n')
            # 表头有区别，Gene ID改为 #Query, Protein ID 改为Protein方便查表计算丰度
            phi = {}
            des = {}
            type = {}
            for line in infile:
                line = line.strip().split("\t")
                if not "Score" in line[0]:
                    query = line[5].rsplit("_1", 1)[0]  # 宏基因组的序列名称需要去掉最后面的_1
                    hit = line[10]
                    id_list = hit.split("#")[1].split("__")
                    pro = hit.split("#")[0]  # 加入protein id
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
                        detail = self.phi.find_one({"PHI": id, "ProteinID": pro})
                        if detail:
                            # pro = detail["ProteinID"]  # 不要这个proteinID
                            hitgene = detail["Gene_name"]
                            ncbi_tax = detail["Pathogen_NCBI_species_Taxonomy_ID"]
                            path_spe = detail["Pathogen_species"]
                            phenotype = detail["Phenotype_of_mutant"]
                            host_des = detail["Host_description"]
                            host_id = detail["Host_NCBI_Taxonomy_ID"]
                            exp = detail["Experimental_host_species"]
                            function = detail["Function"]
                            disease = detail['disease']
                            if id in phi:
                                phi[id] += 1
                            else:
                                phi[id] = 1
                                des[id] = phenotype
                            # pro_list = list(set(pro_list).union(set(pro.split('; '))))
                            hitgene_list = list(set(hitgene_list).union(set(hitgene.split('; '))))
                            # ncbi_tax_list = list(set(ncbi_tax_list).union(set(ncbi_tax.split('; '))))
                            # path_spe_list = list(set(path_spe_list).union(set(path_spe.split('; '))))
                            pathogen_dict = dict(zip(ncbi_tax.split('; '), path_spe.split("; ")))
                            ncbi_tax_list = pathogen_dict.keys()
                            path_spe_list = pathogen_dict.values()
                            phenotype_list = list(set(phenotype_list).union(set(phenotype.split('; '))))
                            # host_des_list = list(set(host_des_list).union(set(host_des.split('; '))))
                            # host_id_list = list(set(host_id_list).union(set(host_id.split('; '))))
                            # exp_list = list(set(exp_list).union(set(exp.split('; '))))
                            id_exp_dict = dict(zip(host_id.split('; '), exp.split('; ')))
                            id_des_dict = dict(zip(host_id.split('; '), host_des.split('; ')))
                            host_des_list = id_des_dict.values()
                            exp_list = id_exp_dict.values()
                            host_id_list = id_exp_dict.keys()
                            function_list = list(set(function_list).union(set(function.split('; '))))
                            disease_list = list(set(disease_list).union(set(disease.split('; '))))
                        else:
                            print("wrong ID:" + id + "\t" + pro)
                    for one_type in phenotype_list:
                        if one_type in type:
                            type[one_type] += 1
                        else:
                            type[one_type] = 1
                    outfile.write(
                        '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, ";".join(id_list),
                                                                                      # ";".join(pro_list),
                                                                                      pro,
                                                                                      ";".join(hitgene_list),
                                                                                      ";".join(ncbi_tax_list),
                                                                                      ";".join(path_spe_list),
                                                                                      ";".join(phenotype_list),
                                                                                      ";".join(host_des_list),
                                                                                      ";".join(host_id_list),
                                                                                      ";".join(exp_list),
                                                                                      ";".join(function_list),
                                                                                      ";".join(disease_list),
                                                                                      act_iden, align_len))
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

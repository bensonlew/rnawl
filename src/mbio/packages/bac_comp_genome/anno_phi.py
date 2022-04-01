# -*- coding: utf-8 -*-
# __author__ = 'gaohao'


import sys
from pymongo import MongoClient
import argparse
import os
from biocluster.config import Config


class phi_anno(object):
    """
    dna phi数据详细注释信息
    ardb.parse.anno.xls
    """
    def phi_by_anno(self, ref, align_table, anno_table):
        db = self.save_dict(ref)
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write(
                'Gene ID\tPHI ID\tProtein ID\tHit Gene Name\tNCBI Tax ID\tPathogen Species\tPhenotype\tHost Descripton\tHost NCBI ID\tHost Species\tGene Function\tIdentity(%)\tEvalue\tDisease\tCoverage(%)\tScore\n')
            infile = infile.readlines()
            for line in infile[1:]:
                line = line.strip().split("\t")
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
                    if id in db.keys():
                        anno = db[id].split("\t")
                        pro = anno[0]
                        hitgene = anno[8]
                        ncbi_tax = anno[6]
                        path_spe = anno[1]
                        phenotype = anno[3]
                        host_des = anno[2]
                        host_id = anno[7]
                        exp = anno[5]
                        function = anno[4]
                        disease = anno[9]
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

    def save_dict(self, file):
        dict = {}
        with open (file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                des = "\t".join(lin[1:])
                dict[lin[0]] = des
        return dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[xml_table]', required=True, help='Input xml table')
    parser.add_argument('-f', metavar='[function_database]', help='phi function', default="function_database")
    parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    db = args.f
    dna_phi_anno = phi_anno()
    dna_phi_anno.phi_by_anno(db, align_table, anno_table)
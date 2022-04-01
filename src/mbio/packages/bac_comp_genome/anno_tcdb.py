# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import sys
import argparse

class meta_tcdb_anno(object):
    """
    dna tcdb数据详细注释信息
    ardb.parse.anno.xls
    """
    def tcdb_by_anno(self, ref, align_table, anno_table):
        db =self.save_dict(ref)
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write(
                'Gene ID\tTCDB ID\tTCDB Description\tTCDB Family\tTCDB Subclass\tTCDB Class\tIdentity(%)\tEvalue\tScore\tAlign_len\n')
            # query_list = []
            # done_list = []
            infile = infile.readlines()
            for line in infile[1:]:
                line = line.strip().split("\t")
                query = line[5]
                hit = line[10]
                id = hit.split("|")[-1]
                evalue = line[1]
                score = line[0]
                align_len = line[11]
                act_iden = line[3]
                if id in db.keys():
                    anno =db[id].split("\t")
                    subclass = anno[3]
                    class_name = anno[1]
                    des = anno[0]
                    family = anno[2]
                    outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, id, des,
                                                                            family, subclass, class_name, act_iden,
                                                                            evalue, score, align_len))
                else:
                    print "wrong ID"

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
    parser.add_argument('-f', metavar='[function_database]', help='tcdb function', default="function_database")
    parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    db = args.f
    meta_tcdb_anno = meta_tcdb_anno()
    meta_tcdb_anno.tcdb_by_anno(db, align_table, anno_table)
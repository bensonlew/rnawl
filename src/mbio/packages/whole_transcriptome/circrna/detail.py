# -*- coding: utf-8 -*-
import pandas as pd

def circ_detail(listid,detail):
    circ_list = list()
    with open(listid) as lists:
        for line in lists:


            id_list,circ_path = line.rstrip().split('\t')
            f = pd.read_table(circ_path)
            f['Junction_site_type'] = f['Junction_site_type'].map(lambda x:x.upper())
            f_detail = f[['circRNA_name','host_gene_id','chr','strand','circRNA_start','circRNA_end','Junction_site_type','circRNA_type']]
            f_detail.rename(columns={'circRNA_name': 'circrna_id', 'Junction_site_type': 'signal','circRNA_start': 'circrna_start','circRNA_end': 'circrna_end','circRNA_type':'circrna_type'}, inplace=True)
            a = f_detail.set_index("circrna_id")
            circ_list.append(a)


        detail_matrix = pd.concat(circ_list,join='outer')
        # detail_matrix_none = detail_matrix.fillna(0)
        detail_matrix_none = detail_matrix.fillna("")


        # detail_matrix_drop = detail_matrix_none.drop_duplicates()
        # 按circ id 去除冗余
        detail_matrix_drop = detail_matrix_none[~detail_matrix_none.index.duplicated(keep='first')]
        detail_matrix_drop.index.name = 'circrna_id'

        detail_matrix_drop.to_csv(detail,sep = '\t')

def main(args):
    circ_detail(args.list_id,args.detail)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='circ detail')
    parser.add_argument('-i', action='store', required=True,
                        help='list_id', dest='list_id')
    parser.add_argument('-o', action='store', required=True,
                        help='detail ', dest='detail')

    args = parser.parse_args()
    main(args)

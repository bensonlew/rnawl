import pandas as pd
import numpy as np

def circ_count(listid,counts):
    circ_list = list()
    with open(listid) as lists:
        for line in lists:
            circ_count = dict()
            id_list,circ_path = line.rstrip().split('\t')
            f = pd.read_table(circ_path)
            seq_id = f['circRNA_name']
            count = f['junction_reads'].astype(np.int64)
            circ_count['seq_id'] = seq_id
            circ_count[id_list] = count
            circcount = pd.DataFrame(circ_count)
            circcount.drop_duplicates(subset=['seq_id'], keep='first', inplace=True)
            a = circcount.set_index("seq_id")
            circ_list.append(a)
        print circ_list

        count_matrix = pd.concat(circ_list,axis=1)
        count_matrix_none = count_matrix.fillna(0)
        count_matrix_int = count_matrix_none.astype(np.int64)
        count_matrix_int.index.name = 'seq_id'

        count_matrix_int.to_csv(counts,sep = '\t')

def main(args):
    circ_count(args.list_id,args.counts)



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='merge_circ_count')
    parser.add_argument('-i', action='store', required=True,
                        help='list_id', dest='list_id')
    parser.add_argument('-o', action='store', required=True,
                        help='counts ', dest='counts')

    args = parser.parse_args()
    main(args)
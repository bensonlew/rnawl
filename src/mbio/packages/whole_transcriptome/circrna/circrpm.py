import pandas as pd

def circ_rpm(listid,rpms):
    circ_list = list()
    with open(listid) as lists:
        for line in lists:
            circ_rpm = dict()
            id_list,circ_path = line.rstrip().split('\t')
            f = pd.read_table(circ_path)
            seq_id = f['circRNA_name']
            rpm = f['RPM']
            circ_rpm['seq_id'] = seq_id
            circ_rpm[id_list] = rpm
            circrpm = pd.DataFrame(circ_rpm)
            circrpm.drop_duplicates(subset=['seq_id'], keep='first', inplace=True)
            a = circrpm.set_index("seq_id")
            circ_list.append(a)


        rpm_matrix = pd.concat(circ_list,axis=1)
        rpm_matrix_none = rpm_matrix.fillna(0)
        rpm_matrix_none.index.name = 'seq_id'

        rpm_matrix_none.to_csv(rpms,sep = '\t')

def main(args):
    circ_rpm(args.list_id,args.rpms)



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='merge_circ_rpm')
    parser.add_argument('-i', action='store', required=True,
                        help='list_id', dest='list_id')
    parser.add_argument('-o', action='store', required=True,
                        help='rpm ', dest='rpms')

    args = parser.parse_args()
    main(args)
from Bio import SeqIO
import pandas as pd
import re


def basetomine(database, minecirc, detail, detail_new):
    circbase_seq_list = list()
    circmine_seq_list = list()
    circmine_id_seq = dict()
    circbase_id_seq = dict()
    for seq_record_circbase in SeqIO.parse(database, "fasta"):
        circbase_seq = seq_record_circbase.seq
        circbase_id = seq_record_circbase.id
        circbase_seq_upper = circbase_seq.upper()
        if circbase_seq_upper not in circbase_id_seq:
            circbase_id_seq[circbase_seq_upper] = [circbase_id]
        else:
            circbase_id_seq[circbase_seq_upper].append(circbase_id)
        circbase_seq_list.append(circbase_seq_upper)

    for seq_record_circmine in SeqIO.parse(minecirc, "fasta"):
        circmine_seq = seq_record_circmine.seq
        circmine_id = seq_record_circmine.id
        (chr, start, end) = re.split('[:-]', circmine_id)
        circmine_name = chr+'_'+start+'_'+end
        circmine_seq_upper = circmine_seq.upper()
        if circmine_seq_upper not in circmine_id_seq:
            circmine_id_seq[circmine_seq_upper] = [circmine_name]
        else:
            circmine_id_seq[circmine_seq_upper].append(circmine_name)
        circmine_seq_list.append(circmine_seq_upper)

    circrna_id = list()
    for i in circmine_seq_list:
        if i in circbase_id_seq.keys():
            circrna_id.extend(circmine_id_seq[i])
    circbase = list()
    df = pd.read_table(detail)
    for l in range(len(df)):
        if df['circrna_id'][l] in circrna_id:
            circbase.append('yes')
        else:
            circbase.append('no')

    df['circbase']=circbase
    df.to_csv(detail_new, index=False, sep='\t', header=True)

def main(args):
    basetomine(args.database, args.minecirc, args.detail, args.detail_new)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='circBase')
    parser.add_argument('-i', action='store', required=True,
                        help='detail', dest='detail')
    parser.add_argument('-d', action='store', required=True,
                        help='database', dest='database')
    parser.add_argument('-m', action='store', required=True,
                        help='minecirc', dest='minecirc')
    parser.add_argument('-o', action='store', required=True,
                        help='output ', dest='detail_new')
    args = parser.parse_args()

    main(args)

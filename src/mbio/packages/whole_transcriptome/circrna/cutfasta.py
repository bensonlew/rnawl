import pandas as pd
from Bio import SeqIO


def cutfasta(details, circfasta, ref, new):
    id_list_ref = list()
    id_list_new = list()
    detail = pd.read_table(details, sep='\t')
    for i in range(len(detail)):
        if detail['circbase'][i] == 'yes':
            a = detail['circrna_id'][i]
            # chr, start, end = a.split('_')
            # b = str(chr) + ':' + str(start) + '-' + str(end)
            b = a
            if b not in id_list_ref:
                id_list_ref.append(b)
        else:
            c = detail['circrna_id'][i]
            # chr, start, end = c.split('_')
            # d = str(chr) + ':' + str(start) + '-' + str(end)
            d = c
            if d not in id_list_new:
                id_list_new.append(d)

    id_set = set()
    sequences = list()
    for record in SeqIO.parse(circfasta, 'fasta'):
        if record.id not in id_set:
            id_set.add(record.id)
            sequences.append(record)
    SeqIO.write(sequences, circfasta, 'fasta')

    fasta = SeqIO.index(circfasta, "fasta")
    handle_ref = open(ref, "w")
    for circ_id1 in id_list_ref:
        handle_ref.write(fasta.get_raw(circ_id1))
    handle_ref.close()

    handle_new = open(new, "w")
    for circ_id2 in id_list_new:
        handle_new.write(fasta.get_raw(circ_id2))
    handle_new.close()


def main(args):
    cutfasta(args.details, args.circfasta, args.ref, args.new)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='cut_fasta')
    parser.add_argument('-d', action='store', required=True,
                        help='detail', dest='details')
    parser.add_argument('-f', action='store', required=True,
                        help='circfasta', dest='circfasta')
    parser.add_argument('-r', action='store', required=True,
                        help='ref', dest='ref')
    parser.add_argument('-n', action='store', required=True,
                        help='new', dest='new')

    args = parser.parse_args()

    main(args)

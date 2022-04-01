from Bio import SeqIO

class palindrome(object):
    def __init__(self, fasta, length_min, length_max, output):
        self.fasta = fasta
        self.min = int(length_min)
        self.max = int(length_max)
        self.output = output
    def read_fasta(self):
        fasta = {}
        fa = SeqIO.parse(self.fasta, 'fasta')
        for record in fa:
            fasta[record.id] = str(record.seq).upper()
        return fasta.items()

    def complement(self, s):
        basecomplement = {
            "A":"T",
            "T":"A",
            "G":"C",
            "C":"G",
            'U':'A',
              }
        letters = list(s)
        letters = [basecomplement[base] for base in letters]
        return ''.join(letters)

    def run(self):
        with open(self.output, 'w') as out:
            out.write('id' + '\t' + 'start' + '\t' + 'end' + '\t' + 'sequence' + '\t' + 'length' + '\n')

            for seqname, seq in self.read_fasta():
                for i in range(len(seq)):
                    for j in range(self.min, self.max):
                        selectseq = seq[i:i + j]
                        if 'U' in selectseq:
                            selectseq_new = selectseq.replace('U', 'T')
                        else:
                            selectseq_new = selectseq
                        reverse_complement_seq = self.complement(selectseq_new)[::-1]
                        if selectseq_new == reverse_complement_seq:
                            maxlen = len(selectseq)
                            start = i + 1
                            end = i + j + 1
                            CPSeq = selectseq
                            out.write(seqname + '\t' + str(start) + '\t' + str(end) + '\t' + CPSeq +'\t' + str(maxlen) + '\n')
                            start = None
                            end = None
                            CPSeq = None

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='This script is used to identify pirna from piRbase database.')
    parser.add_argument('-f', type=str, help='fasta file')
    parser.add_argument('-mi', type=int, help='min length')
    parser.add_argument('-ma', type=int, help='max length')
    parser.add_argument('-o', type=str, help='output path')
    args = parser.parse_args()
    palindrome = palindrome(args.f, args.mi, args.ma, args.o)
    palindrome.run()
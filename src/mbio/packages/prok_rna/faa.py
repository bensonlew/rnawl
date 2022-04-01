from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import argparse


def run_faa(cds, faa):
    with open(cds, 'r') as fa_r, open(faa, 'w') as faa_w:
        for block in fa_r.read().split('\n>'):
            block = block.lstrip('>').split('\n')
            coding_dna = Seq(''.join(block[1:]), IUPAC.ambiguous_dna)
            protein = coding_dna.translate(gap='X')
            faa_w.write('>' + block[0].strip() + '\n' + str(protein) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', type=str, metavar='cds_file', required=True)
    parser.add_argument('-o', type=str,metavar='output file',required=True)
    args = parser.parse_args()
    run_faa(args.c, args.o)
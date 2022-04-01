# coding=utf-8
from Bio import SeqIO
import unittest
author = 'gdq'


def split_fasta(fasta, chunk_size=100, out_prefix='my_fasta'):
    with open(fasta) as f:
        seq_num = 0
        chunk_num = 0
        out_name = out_prefix + '_' + str(chunk_num)
        file_objects = {chunk_num: open(out_name, 'w')}
        for line in f:
            if line.startswith('>'):
                seq_num += 1
                chunk_num = seq_num//(chunk_size+1)
                if chunk_num not in file_objects:
                    file_objects[chunk_num-1].close()
                    out_name = out_prefix + '_' + str(chunk_num)
                    file_objects[chunk_num] = open(out_name, 'w')
            file_objects[chunk_num].write(line)
        else:
            file_objects[chunk_num].close()


def split_fasta2(fasta, chunk_size=100):
    """
    """
    line = 1
    i = 1
    w = open('fasta_1', 'wb')
    for seq_record in SeqIO.parse(fasta, "fasta"):
        if line <= chunk_size:
            w.write('>{}\n{}\n'.format(seq_record.id, seq_record.seq))
            line += 1
        else:
            i += 1
            w.close()
            line = 1
            w = open('fasta_%s' % i, 'wb')
    w.close()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_my(self):
        split_fasta('ref_index.transcripts.fa')

    def test_other(self):
        split_fasta2('ref_index.transcripts.fa')


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_my'))
    suite.addTest(TestFunction('test_other'))
    unittest.TextTestRunner(verbosity=2).run(suite)

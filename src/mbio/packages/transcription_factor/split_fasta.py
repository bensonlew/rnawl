# coding=utf-8
import argparse
import unittest
author = 'gdq'


def split_fasta(fasta, chunk_size=100, out_prefix='split', strip_char=".*", replace=('.', 'X')):
    with open(fasta) as f:
        seq_num = 0
        chunk_num = 0
        out_name = out_prefix + '_' + str(chunk_num) + '.fa'
        file_objects = {chunk_num: open(out_name, 'w')}
        for line in f:
            if line.startswith("#") or (not line.strip()):
                continue
            if line.startswith('>'):
                chunk_num = seq_num // chunk_size
                seq_num += 1
                if chunk_num not in file_objects:
                    file_objects[chunk_num-1].close()
                    out_name = out_prefix + '_' + str(chunk_num) + '.fa'
                    file_objects[chunk_num] = open(out_name, 'w')
            else:
                line = line.strip().strip(strip_char) + '\n' # 去除两端的*或.
                line = line.replace(replace[0], replace[1])  # 替换'.'为X
            file_objects[chunk_num].write(line)
        else:
            file_objects[chunk_num].close()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_my(self):
        split_fasta('ref_index.transcripts.fa')


if __name__ == '__main__':
    # suite = unittest.TestSuite()
    # suite.addTest(TestFunction('test_my'))
    # unittest.TextTestRunner(verbosity=2).run(suite)
    parser = argparse.ArgumentParser(description="split fasta file")
    parser.add_argument('-f', help="fasta file to be split")
    parser.add_argument('-size', metavar="split_size", type=int, help="specify sequence number of each split file")
    parser.add_argument('-prefix', help="prefix of each split file name")
    args = parser.parse_args()
    split_fasta(args.f, args.size, args.prefix)



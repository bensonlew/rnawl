# -*- coding: utf-8 -*-
# __author__ : 'shicaiping'
# __date__: 20200520


from Bio import SeqIO
from Bio.Seq import reverse_complement
import re
from string import maketrans


class ExciseCandidate(object):
    """
    This script excised potential microRNA precursor sequences
    from a genome using the positions of aligned reads as guidelines.
    """

    def __init__(self, file_fasta, file_blast_parsed, precursor_length, filtered_precursors):
        """
        The fasta file given as input should be the genome in question,
        and the file in blastparsed format should contain the alignments.
        the precursor_length is the length of excised precursors.
        """
        self.file_fasta = file_fasta
        self.file_blast_parsed = file_blast_parsed
        self.precursor_length = int(precursor_length)
        self.filtered_precursors = filtered_precursors
        self.hash_fasta = dict()
        self.hash_align = dict()

    def run(self):
        self.parse_file_fasta(self.file_fasta, self.hash_fasta)
        self.parse_file_blast_parsed(self.file_blast_parsed)
        self.excise()

    def parse_file_fasta(self, file_fasta, hash_fasta):
        for seq_record in SeqIO.parse(file_fasta, "fasta"):
            seq_id = seq_record.id
            seq_squence = seq_record.seq
            hash_fasta[seq_id] = seq_squence

    def parse_file_blast_parsed(self, file_blast_parsed):
        pattern = r'^(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)$'
        with open(file_blast_parsed, "r") as f:
            num = 0
            for line in f:
                num += 1
                print "Process line{}: {}".format(num, line)
                re_result = re.match(pattern, line.strip())
                query = re_result.group(1)
                query_lng = int(re_result.group(2))
                query_beg = int(re_result.group(3))
                query_end = int(re_result.group(4))
                subject = re_result.group(5)
                subject_lng = int(re_result.group(6))
                subject_beg = int(re_result.group(7))
                subject_end = int(re_result.group(8))
                e_value = re_result.group(9)
                pid = re_result.group(10)
                bitscore = re_result.group(11)
                other = re_result.group(12)
                strand = self.find_strand(other)
                if subject not in self.hash_fasta:
                    print ("{} not in genome file.".format(subject))
                    continue
                else:
                    self.insertfeature(query, query_beg, query_end, query_lng, subject, subject_beg, subject_end,
                                       subject_lng, strand)

    def excise(self):
        short_length = self.precursor_length - 23
        subjects = sorted(self.hash_align.keys())
        for subject in subjects:
            count = 0
            strands = sorted(self.hash_align[subject].keys())
            for strand in strands:
                if strand == "+":
                    seq = self.hash_fasta[subject]
                else:
                    seq = reverse_complement(self.hash_fasta[subject])
                seq_lng = len(seq)
                subject_begs = sorted(self.hash_align[subject][strand].keys())
                for subject_beg in subject_begs:
                    queries = sorted(self.hash_align[subject][strand][subject_beg].keys())
                    for query in queries:
                        subject_end = self.hash_align[subject][strand][subject_beg][query]["subject_end"]
                        if subject_end - subject_beg > 30:
                            self.excise_position(subject, seq, seq_lng, strand, subject_beg - 22, subject_end + 22,
                                                 count)
                            count += 1
                        else:
                            self.excise_position(subject, seq, seq_lng, strand, subject_beg - 22,
                                                 subject_beg + short_length, count)
                            count += 1
                            self.excise_position(subject, seq, seq_lng, strand, subject_end - short_length,
                                                 subject_end + 22, count)
                            count += 1

    def excise_position(self, subject, seq, seq_lng, strand, excise_beg_old, excise_end_old, count):
        length_account = self.precursor_length + 30
        excise_beg = 1 if 1 > excise_beg_old else excise_beg_old
        excise_end = seq_lng if seq_lng < excise_end_old else excise_end_old
        excise_lng = excise_end - excise_beg + 1
        if length_account < excise_lng:
            return
        seq_sub = seq[(excise_beg - 1):(excise_beg - 1 + excise_lng)]
        with open(self.filtered_precursors, "a") as w:
            w.write(">" + subject + "_" + str(count) + " strand:" + strand + " excise_beg:" + str(excise_beg) +
                    " excise_end:" + str(excise_end) + "\n")
            w.write(str(seq_sub) + "\n")
        return

    def find_strand(self, strand):
        strand_result = "+"
        if re.search(r'-', strand):
            strand_result = "-"
        if re.search(r'minus', strand, re.I):
            strand_result = "-"
        return strand_result

    def contained(self, begin1, end1, begin2, end2):
        if begin2 <= begin1 and end1 <= end2:
            return 1
        else:
            return 0

    def overlapping(self, begin1, end1, begin2, end2):
        if ((begin1 <= begin2 <= end1 + 1) or (begin1 <= end2 + 1 and end2 <= end1) or self.contained(
                begin1, end1, begin2, end2) or self.contained(begin2, end2, begin1, end1)):
            return 1
        else:
            return 0

    def find_overlap(self, begin1, end1, begin2, end2):
        self.test_begin_end(begin1, end1, begin2, end2)
        begin = begin1 if begin1 < begin2 else begin2
        end = end1 if end1 > end2 else end2
        return begin, end

    def test_begin_end(self, begin1, end1, begin2, end2):
        if not(begin1 < end1 and begin2 < end2):
            print "begin1: {}, end1: {}, begin2: {}, end2: {}".format(begin1, end1, begin2, end2)
            raise Exception("Begin positions must be numerically smaller than endpositions.")

    def delete_feature(self, subject, strand, subject_beg, query):
        begins = len(self.hash_align[subject][strand][subject_beg].keys())
        if begins == 1:
            del (self.hash_align[subject][strand][subject_beg])
        else:
            del (self.hash_align[subject][strand][subject_beg][query])

    def insertfeature(self, query, query_beg, query_end, query_lng, subject, subject_beg, subject_end, subject_lng,
                      strand):
        if self.hash_align and subject in self.hash_align and strand in self.hash_align[subject]:
            prev_begs = sorted(self.hash_align[subject][strand].keys())
            for prev_beg in prev_begs:
                distance = subject_beg - prev_beg
                if distance > 1000:
                    continue
                if distance < -1000:
                    break
                prev_queries = sorted(self.hash_align[subject][strand][prev_beg].keys())
                for prev_query in prev_queries:
                    prev_end = self.hash_align[subject][strand][prev_beg][prev_query]["subject_end"]
                    flank_beg = 1 if 1 > (subject_beg - 30) else (subject_beg - 30)
                    flank_end = subject_lng if subject_lng < (subject_end + 30) else subject_end + 30
                    if self.overlapping(flank_beg, flank_end, prev_beg, prev_end):
                        new_beg, new_end = self.find_overlap(subject_beg, subject_end, prev_beg, prev_end)
                        subject_beg = new_beg
                        subject_end = new_end
                        self.delete_feature(subject, strand, prev_beg, prev_query)
        if subject not in self.hash_align:
            self.hash_align[subject] = dict()
        if strand not in self.hash_align[subject]:
            self.hash_align[subject][strand] = dict()
        if subject_beg not in self.hash_align[subject][strand]:
            self.hash_align[subject][strand][subject_beg] = dict()
        if query not in self.hash_align[subject][strand][subject_beg]:
            self.hash_align[subject][strand][subject_beg][query] = dict()
        self.hash_align[subject][strand][subject_beg][query]["subject_end"] = subject_end
        self.hash_align[subject][strand][subject_beg][query]["subject_lng"] = subject_lng
        self.hash_align[subject][strand][subject_beg][query]["query_beg"] = query_beg
        self.hash_align[subject][strand][subject_beg][query]["query_end"] = query_end
        self.hash_align[subject][strand][subject_beg][query]["query_lng"] = query_lng


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='This script is used to convert rna-seq expression units.')
    parser.add_argument('file_fasta', type=str,
                        help='The fasta file given as input should be the genome in question.')
    parser.add_argument('file_blast_parsed', type=str,
                        help='The file in blastparsed format should contain the alignments.')
    parser.add_argument('precursor_length', type=int, default=300, help='The precursor_length is the length of '
                                                                        'excised precursors.')
    parser.add_argument('filtered_precursors', type=str, help='output filtered_precursors.fa')

    args = parser.parse_args()

    ExciseCandidate = ExciseCandidate(args.file_fasta, args.file_blast_parsed, args.precursor_length,
                                      args.filtered_precursors)
    ExciseCandidate.run()

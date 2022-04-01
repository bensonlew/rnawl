from Bio.Blast import NCBIXML
import re
import os

class id2nr(object):
    def __init__(self, xml_in, id, output):
        self.xml_in = xml_in
        self.id = id
        self.output = output

    def blast_nr_xml(self):

        records = NCBIXML.parse(open(self.xml_in))
        blast_dict = dict()
        for record in records:
            for alignment in record.alignments:
                query_id = str(re.split(' ', record.query, maxsplit=1)[0])
                if query_id not in blast_dict:
                    blast_dict[query_id] = [str(alignment.hit_id), str(alignment.hit_def)]
                else:
                    break

        with open(self.id, 'r') as ens_id, open(self.output, 'w') as out:
            for line in ens_id.readlines():
                if line.strip() in blast_dict:
                    out.write(line.strip() + '\t' + blast_dict[line.strip()][0] + '\t' + blast_dict[line.strip()][1] + '\n')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='This script is used to id to NR description')
    parser.add_argument('-xml', type=str, help='blast_nr.xml')
    parser.add_argument('-id', type=str, help='id file')
    parser.add_argument('-o', type=str, help='output path')
    args = parser.parse_args()
    id2nr = id2nr(args.xml, args.id, args.o)
    id2nr.blast_nr_xml()
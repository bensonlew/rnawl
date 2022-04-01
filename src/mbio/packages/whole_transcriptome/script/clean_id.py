import sys

from Bio import SeqIO

list_file = sys.argv[1]

fastqs = [line.strip().split('\t')[0] for line in open(list_file) if line.strip()]

for fastq in fastqs:
    sequences = list()
    for record in SeqIO.parse(fastq, 'fastq'):
        old_id = record.id
        new_id = old_id[:old_id.rfind('.')]
        record.id = record.id.replace(old_id, new_id)
        record.name = record.name.replace(old_id, new_id)
        record.description = record.description.replace(old_id, new_id)
        sequences.append(record)
    SeqIO.write(sequences, fastq, 'fastq')
    print fastq

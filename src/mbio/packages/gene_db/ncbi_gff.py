# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from BCBio import GFF
import urllib
import sys
import gzip
from mbio.api.database.gene_db import genome
import argparse


class NcbiGff(object):
    def __init__(self, genome_acc):
        self.genome_api = genome.Genome(None)
        self.genome_acc = genome_acc
        self.gene_schema = ['sg_gene_id', 'gene_id', 'gene_type', 'GeneID', 'HGNC', 'ID', 'Name',
                            'description', 'gbkey', 'gene', 'gene_biotype']
        self.tran_schema = self.gene_schema + ['sg_tran_id', 'tran_id', 'transcript_id', 'protein_id', 'Genbank']

    def pipe(self, gff, out_pre):
        genome_dict = self.get_genome_acc(self.genome_acc)
        # print genome_dict
        self.gff2table(gff, out_pre, genome_dict)

    def get_genome_acc(self, genome_acc):
        # print "genome acc is {}".format(genome_acc)
        chr_records = self.genome_api.get_genome_chr(acc_id=genome_acc,
                                                     assembly_unit={"$in": ["Primary Assembly", "non-nuclear", 'C57BL/6J']})
        chr_list = list()
        for record in chr_records:
            # print record
            # 大鼠的染色体前加了chr
            if record['sequence_name'].startswith("chr"):
                    record['sequence_name'] = record['sequence_name'].split("chr")[1]
            if record['sequence_role'] == 'assembled-molecule':
                # print record
                chr_list.append((record['refseq_accn'], record['sequence_name']))
            elif record['sequence_role'] in ['unlocalized-scaffold', 'unplaced-scaffold']:
                chr_list.append((record['refseq_accn'], record['genbank_accn']))
        return dict(chr_list)


    def gff2table(self, in_file, out_pre, genome_dict):
        print genome_dict

        in_handle = gzip.open(in_file)
        rec_list = []
        for rec in GFF.parse(in_handle):
            try:
                rec_list.append(rec)
            except ValueError:
                print "坐标错误"

        in_handle.close()
        '''
        gene2id_handle = open(out_file + ".gene2enterz", 'w')
        gene2name_handle = open(out_file + ".tran2name", 'w')
        gene2des_handle = open(out_file + ".tran2des", 'w')
        tran_id_list = []
        '''
        gene_list = list()
        tran_list = list()

        strand_dict= {"1": "1", "-1": "0"}

        for seq_record in rec_list:
            chromosome_acc = seq_record.id
            if chromosome_acc in genome_dict:
                chromosome_name = genome_dict[chromosome_acc]
            else:
                continue

            for gene_feature in seq_record.features:
                if gene_feature.type.lower() == "gene" or gene_feature.type.lower() == "pseudogene":
                    gene_dict =dict()
                    gene_dict['sg_gene_id'] = "_".join([chromosome_name,
                                                    str(gene_feature.location.start + 1),
                                                    str(gene_feature.location.end),
                                                    strand_dict[str(gene_feature.location.strand)]])
                    gene_dict['gene_id'] = gene_feature.id
                    gene_dict['gene_type'] = gene_feature.type

                    for k,v in gene_feature.qualifiers.items():
                        if k == "Dbxref":
                            for i in v:
                                gene_dict[i.split(":")[0]] = ":".join(i.split(":")[1:])
                        else:
                            gene_dict[k] = "|".join(v)

                    gene_list.append(gene_dict)

                    for rna_feature in gene_feature.sub_features:
                        rna_dict = gene_dict
                        rna_dict['sg_tran_id'] = "_".join([chromosome_name,
                                                    str(gene_feature.location.start + 1) + "." + str(gene_feature.location.end),
                                                    strand_dict[str(gene_feature.location.strand)]])
                        rna_dict['tran_id'] = rna_feature.id
                        rna_dict['exons'] = list()


                        for k,v in rna_feature.qualifiers.items():
                            if k == "Dbxref":
                                for i in v:
                                    if i.split(":")[0] not in ['GeneID', 'HGNC']:
                                        rna_dict[i.split(":")[0]] = ":".join(i.split(":")[1:])
                            elif k in ["transcript_id", "protein_id"]:
                                rna_dict[k] = "|".join(v)
                            elif k == "gbkey":
                                rna_dict["tran_gbkey"] = "|".join(v)
                            else:
                                pass

                        if rna_feature.type.lower() in ["cds", 'exon']:
                            pass
                        else:
                            for exon_feature in rna_feature.sub_features:
                                if exon_feature.type == "exon":
                                    rna_dict['exons'].append(str(exon_feature.location.start + 1) + "." + str(exon_feature.location.end))
                                if exon_feature.type.lower() == "cds":
                                    for k,v in exon_feature.qualifiers.items():
                                        if k == "protein_id":
                                            rna_dict[k] = "|".join(v)
                            if len(rna_dict['exons']) == 0:
                                pass
                            else:
                                rna_dict['sg_tran_id'] =  "_".join([chromosome_name] +
                                                                    rna_dict['exons'] +
                                                                    [strand_dict[str(gene_feature.location.strand)]])
                        tran_list.append(rna_dict)
                else:
                    pass
        with open(out_pre + ".gene.tsv", 'w') as fo:
            fo.write("\t".join(self.gene_schema) + "\n")
            for gene_dict in gene_list:
                fo.write("\t".join([gene_dict.get(field, "") for field in self.gene_schema]) + "\n")
        with open(out_pre + ".tran.tsv", 'w') as fo:
            fo.write("\t".join(self.tran_schema) + "\n")
            for tran_dict in tran_list:
                fo.write("\t".join([tran_dict.get(field, "") for field in self.tran_schema]) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-acc', type=str, required=True,
                        help='genome acc')
    parser.add_argument('-gff', type=str, required=True,
                        help='refseq gff')
    parser.add_argument('-out_pre', type=str, required=False, default="table",
                        help="output pre")

    # ----------------------------------------------------------------------------------------------
    args = parser.parse_args()
    mg = NcbiGff(args.acc)
    mg.pipe(args.gff, args.out_pre)

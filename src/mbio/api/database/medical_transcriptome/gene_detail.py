# -*- coding: utf-8 -*-
# __author__ 'gudeqing,qinjincheng'

from mbio.api.database.medical_transcriptome.api_base import ApiBase
from biocluster.api.database.base import report_check
import re
from collections import defaultdict
import datetime
import unittest

class GeneDetail(ApiBase):
    def __init__(self, bind_object):
        super(GeneDetail, self).__init__(bind_object)

    @report_check
    def add_gene_detail(self, refrna_seqdb, t2g_file, txpt_bed, txpt_fa, gene_bed, gene_fa,
                        biomart_file, biomart_type, species_urls, new_cds, new_pep):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        name = 'Gene_detail_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        main_info = dict([
            ('task_id', task_id),
            ('project_sn', project_sn),
            ('desc', 'gene_detail_information'),
            ('created_ts', created_ts),
            ('status', 'start'),
            ('name', name),
            ('refrna_seqdb', refrna_seqdb)
        ])
        main_id = self.create_db_table('sg_genes', [main_info])

        t2g_dict = dict([x.strip().split()[:2] for x in open(t2g_file) if x.strip()])
        txpt_loc = self.bed_to_loc(txpt_bed)
        txpt_seq = self.fasta2dict(txpt_fa)
        gene_loc = self.bed_to_loc(gene_bed)
        gene_seq = self.fasta2dict(gene_fa)
        biomart_info, pep2trans = self.parse_biomart(biomart_file, biomart_type)
        new_cds_pep = self.get_new_cds_pep(new_cds, new_pep) if new_cds and new_pep else dict()
        gene_detail_list = list()
        g2t = dict()
        for t, g in t2g_dict.items():
            g2t.setdefault(g, set())
            g2t[g].add(t)
        transcripts_total_length = 0
        for g, t_set in g2t.items():
            trans_info = list()
            for t in t_set:
                if t not in txpt_loc:
                    self.bind_object.set_error('%s not found in %s', variables=(t, txpt_bed), code="55600107")
                t_length = len(txpt_seq[t])
                transcripts_total_length += t_length
                trans_info.append(dict(
                    transcript_id=t,
                    start=txpt_loc[t]['start'],
                    end=txpt_loc[t]['end'],
                    length=t_length,
                ))
            if g not in gene_loc:
                self.bind_object.set_error('%s not found in %s', variables=(g, gene_bed), code="55600108")
            start = gene_loc[g]['start']
            end = gene_loc[g]['end']
            location = '{}-{}'.format(start, end)
            gene_length = len(gene_seq[g])
            if g in biomart_info:
                gene_type = biomart_info[g]['gene_type'][0]
            else:
                gene_type = 'predicted_non-coding'
                for t in t_set:
                    if t in new_cds_pep:
                        gene_type = 'protein_coding'
                        break
            if 'ensembl' in species_urls:
                gene_ensembl = "{}/Gene/Summary?g={}".format(species_urls[:-11], g)
            else:
                gene_ensembl = None
            gene_detail_list.append(dict(
                gene_id=g,
                transcripts=trans_info,
                start=start,
                end=end,
                strand=gene_loc[g]['strand'],
                chrom=gene_loc[g]['chr'],
                location=location,
                gene_length=gene_length,
                gene_type=gene_type,
                gene_ensembl=gene_ensembl,
                main_id=main_id,
                transcript_num=len(trans_info),
            ))
        self.create_db_table('sg_genes_detail', gene_detail_list)
        self.update_db_record('sg_genes', main_id, status='end', main_id=main_id,
                              transcripts_total_length=transcripts_total_length)

    @staticmethod
    def parse_biomart(biomart_path, biomart_type):
        if biomart_type == "type1":
            gene_id_ind = 0
            trans_id_ind = 1
            gene_name_ind = 2
            chromosome_ind = 8
            gene_type_ind = 16
            desc_ind = 7
            strand_ind = 11
            start_ind = 9
            end_ind = 10
            pep_id_ind = 6
        elif biomart_type == 'type2':
            gene_id_ind = 0
            trans_id_ind = 1
            gene_name_ind = 2
            chromosome_ind = 6
            gene_type_ind = 14
            desc_ind = 5
            strand_ind = 9
            start_ind = 7
            end_ind = 8
            pep_id_ind = 4
        elif biomart_type == 'type3':
            gene_id_ind = 0
            trans_id_ind = 1
            gene_name_ind = 0
            chromosome_ind = 4
            gene_type_ind = 12
            desc_ind = 3
            strand_ind = 7
            start_ind = 5
            end_ind = 6
            pep_id_ind = 2
        else:
            self.bind_object.set_error('biomart_type should be one of type1, type2, type3', code="55600109")

        biomart_info = dict()
        pep2transcript = dict()
        with open(biomart_path) as f:
            for line in f:
                if not line.strip():
                    continue
                line = line.replace('\t\t', '\t-\t')
                tmp_list = line.strip('\n').split("\t")
                gene_id = tmp_list[gene_id_ind]
                trans_id = tmp_list[trans_id_ind]
                gene_name = tmp_list[gene_name_ind]
                if biomart_type == 'type3':
                    gene_name = '-'
                chromosome = tmp_list[chromosome_ind]
                gene_type = tmp_list[gene_type_ind]
                desc = tmp_list[desc_ind]
                strand_tmp = tmp_list[strand_ind]
                if strand_tmp == "1":
                    strand = "+"
                elif strand_tmp == "-1":
                    strand = "-"
                elif strand_tmp == "0":
                    strand = "."
                else:
                    strand = strand_tmp
                start = tmp_list[start_ind]
                end = tmp_list[end_ind]
                pep_id = tmp_list[pep_id_ind]

                biomart_info.setdefault(gene_id, defaultdict(list))
                biomart_info[gene_id]['trans_id'].append(trans_id)
                biomart_info[gene_id]['gene_name'].append(gene_name)
                biomart_info[gene_id]['chromosome'].append(chromosome)
                biomart_info[gene_id]['gene_type'].append(gene_type)
                biomart_info[gene_id]['description'].append(desc)
                biomart_info[gene_id]['pep_id'].append(pep_id)
                biomart_info[gene_id]['strand'].append(strand)
                biomart_info[gene_id]['start'].append(start)
                biomart_info[gene_id]['end'].append(end)
                pep2transcript[pep_id] = trans_id

        if not biomart_info:
            self.bind_object.set_error("biomart information is None", code="55600110")
        else:
            # raw check
            if (not start.isdigit()) or (not end.isdigit()):
                self.bind_object.set_error('we find "start" or "end" is not digit. Maybe biomart_type is wrong', code="55600111")
        print('Information of {} genes was parsed from biomart file'.format(len(biomart_info)))
        return biomart_info, pep2transcript

    @staticmethod
    def bed_to_loc(gene_bed_path):
        gene_info = dict()
        with open(gene_bed_path, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                line = line.strip('\n').split("\t")
                if gene_info.get(line[3]) is None:
                    gene_info[line[3]] = {"chr": line[0], "strand": line[5], "start": str(int(line[1]) + 1),
                                          "end": line[2]}
                else:
                    print("lalallala{}".format(str(line[3])))
                    # self.bind_object.set_error("%s中基因%s在多行出现", variables=(gene_bed_path, line[3]), code="55600112")
        return gene_info

    @staticmethod
    def fasta2dict(fasta_file):
        seq = dict()
        match_name = re.compile(r'>([^\s]+)').match
        with open(fasta_file, 'r') as fasta:
            j = 0
            seq_id, sequence = '', ''
            for line in fasta:
                if line.startswith('#') or not line.strip():
                    continue
                if line.startswith('>'):
                    j += 1
                    if j > 1:
                        seq[seq_id] = sequence
                        sequence = ''
                    # get seq name
                    seq_id = match_name(line).group(1)
                    # gene id和转录本id可能包含括号
                    if '(+)' in seq_id:
                        seq_id = seq_id.split('(+)')[0]
                    elif '(-)' in seq_id:
                        seq_id = seq_id.split('(-)')[0]
                    elif '(.)' in seq_id:
                        seq_id = seq_id.split('(.)')[0]
                    elif '|' in seq_id:
                        seq_id = seq_id.split('|')[0]
                    else:
                        pass
                else:
                    sequence += line.strip()
            else:
                # save the last sequence
                seq[seq_id] = sequence
        if not seq:
            print('提取序列信息为空')
        print("从{}共统计出{}条序列信息".format(fasta_file, len(seq)))
        return seq

    @staticmethod
    def parse_seq_file(seq_file):
        with open(seq_file, 'r') as f:
            j = 0
            seq_tag = tuple()
            sequence = ""
            for line in f:
                if not line.strip():
                    continue
                if line.startswith('>'):
                    j += 1
                    if j > 1:
                        yield seq_tag, sequence
                    seq_tag = line.strip()
                    sequence = ''
                else:
                    sequence += line.strip()
            else:
                yield seq_tag, sequence

    def get_new_cds_pep(self, cds_file, pep_file):
        pep_parser = self.parse_seq_file(pep_file)
        pep_dict = dict()
        for pep_desc, pep_seq in pep_parser:
            pep_dict[pep_desc] = pep_seq
        t2cds_pep = dict()
        match = re.compile(r'>(.*?)::(.*?)::.*type:(.*?)\s+len:\d+.*:(\d+-\d+\(.\)).*').match
        cds_parser = self.parse_seq_file(cds_file)
        for cds_desc, cds_seq in cds_parser:
            g_id, t_id, _type, cds_pos = match(cds_desc).groups()
            pep_seq = pep_dict[cds_desc]
            seq_id = re.match(r'>(.*?)\s.*', cds_desc).groups()[0]
            t2cds_pep.setdefault(t_id, list())
            t2cds_pep[t_id].append((seq_id, _type, cds_seq, pep_seq))
        return t2cds_pep

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """



    def test_gene_detail(self):
        from mbio.workflows.medical_transcriptome.medical_transcriptome_test_api import MedicalTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
            "id": "medical_transcriptome",
            "project_sn": "medical_transcriptome",
            "type": "workflow",
            "name": "medical_transcriptome.medical_transcriptome_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = MedicalTranscriptomeTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("medical_transcriptome.gene_detail")
        refrna_seqdb = "/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/Detail/output/refrna_seqs.db"
        t2g_file = "/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/RefrnaAssemble/output/NewTranscripts/trans2gene"
        txpt_fa = "/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/RefrnaAssemble/output/NewTranscripts/all_transcripts.fa"
        new_cds ="/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/AnnotOrfpfam/output/new_transcripts.fa.transdecoder.cds"
        new_pep = "/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/AnnotOrfpfam/output/new_transcripts.fa.transdecoder.pep"
        txpt_bed = "/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/GeneFa/output/transcript.bed"
        gene_bed = "/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/GeneFa/output/gene.bed"
        gene_fa = "/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/GeneFa/output/gene.fasta"
        biomart_file ="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/biomart/biomart.txt"
        biomart_type ="type1"
        species_urls = "http://asia.ensembl.org/Homo_sapiens/Info/Index"
        wf.test_api.add_gene_detail(refrna_seqdb, t2g_file, txpt_bed, txpt_fa, gene_bed, gene_fa,
                            biomart_file, biomart_type, species_urls, new_cds, new_pep)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_gene_detail'))
    unittest.TextTestRunner(verbosity=2).run(suite)

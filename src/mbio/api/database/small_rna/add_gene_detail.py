# -*- coding: utf-8 -*-
# __author__ = 'shicaiping, qinjincheng'

from mbio.api.database.small_rna.api_base import ApiBase
import sqlite3
import re
from collections import defaultdict
import datetime
import unittest
import os
import pandas as pd

class AddGeneDetail(ApiBase):
    def __init__(self, bind_object):
        super(AddGeneDetail, self).__init__(bind_object)

    def add_gene_detail(self, db_path, gene_fasta, transcript_fasta, t2g, biomart_file, biomart_type,
                        cds_fasta, pep_fasta, gene_bed, transcript_bed, species_urls, gene_stat):
        '''
        Build a seqs.db (sqlite3 file) for querying and dump related information into mongo.
        :param db_path: absolute path of sqlite3 file
        :param gene_fasta: absolute path of FASTA file of gene
        :param transcript_fasta: absolute path of FASTA file of transcript
        :param t2g: absolute path of t2g file
        :param biomart_file: absolute path of biomart file
        :param biomart_type: the type of biomart file
        :param cds_fasta: absolute path of FASTA file of cds
        :param pep_fasta: absolute path of FASTA file of pep
        :param gene_bed: absolute path of BED file of gene
        :param transcript_bed: absolute path of BED file of transcript
        :param species_urls: the url of specified species in ensembl
        :return:
        '''

        # prepare variables for build_seq_db
        gene_seq = self.fasta2dict(gene_fasta)
        trans_seq = self.fasta2dict(transcript_fasta)
        t2g_dict = dict([x.strip().split()[:2] for x in open(t2g) if x.strip()])
        # g2t_dict = dict()
        # with open(t2g, 'r') as t:
        #     for i in t.readlines():
        #         trans_id, gene_id = i.strip().split()
        #         if gene_id not in g2t_dict:
        #             g2t_dict[gene_id] = [trans_id]
        #         else:
        #             g2t_dict[gene_id].append(trans_id)

        biomart_gene_detail, pep2trans = self.biomart(biomart_file, biomart_type)
        known_cds = self.get_cds_seq(cds_fasta)
        known_pep = self.get_pep_seq(pep_fasta, pep2trans)


        # gene_id_list = list()
        # trans_num = list()
        # cds_num = list()
        # pep_num = list()
        # for i in gene_seq.keys():
        #     count_cds = 0
        #     count_pep = 0
        #     gene_id_list.append(i)
        #     trans_num.append(len(g2t_dict[i]))
        #     for j in g2t_dict[i]:
        #         if j in known_cds:
        #             count_cds += 1
        #         if j in known_pep:
        #             count_pep += 1
        #     cds_num.append(count_cds)
        #     pep_num.append(count_pep)
        # gene_stat_df = pd.DataFrame({'gene_id': gene_id_list, 'tran_num': trans_num, 'cds_num': cds_num, 'pep_num': pep_num})

        # build_seq_db and check return value
        boolean = self.build_seq_db(db_path, gene_seq, trans_seq, t2g_dict, known_cds, known_pep)
        if boolean:
            self.bind_object.logger.info('succeed in build_seq_db')
        else:
            self.bind_object.set_error('fail to build_seq_db')

        # prepare variables for dump_mongo
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        gene_loc = self.gene_location(gene_bed)
        trans_loc = self.gene_location(transcript_bed)

        # dump_mongo and check return value
        boolean = self.dump_mongo(task_id, project_sn, db_path, t2g_dict, trans_loc, transcript_bed, trans_seq,
                                  gene_loc, gene_bed, gene_seq, biomart_gene_detail, species_urls)
        if boolean:
            self.bind_object.logger.info('succeed in dump_mongo')
        else:
            self.bind_object.set_error('fail to dump_mongo')


        # seq_stat

        self.add_seq_stat(gene_stat)

    def fasta2dict(self, fasta_file):
        '''
        Get sequence from provided FASTA file.
        :param fasta_file: absolute path of FASTA file
        :return: dict, {seq_id: sequence}
        '''
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
                    # gene id and transcript id may contain bracket
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
            self.bind_object.set_error('sequence information is None')
        self.bind_object.logger.info('information of {} sequences were parsed from {}'.format(len(seq), fasta_file))
        return seq

    def biomart(self, biomart_file, biomart_type):
        '''
        Get description and gene_type of known gene from provided biomart file.
        :param biomart_file: absolute path of biomart file, this file must be tab separated
        :param biomart_type: the type of biomart file
        :return: dict1, {gene_id: {"trans_id": [trans_id], "gene_name": [gene_name], "chromosome": [chromosome],
                                   "gene_type": [gene_type], "description": [desc], "strand": [strand],
                                   "pep_id": [pep_id], "start": [start], "end": [end]}}
                 dict2, {pep_id: trans_id}
        '''
        if biomart_type == 'type1':
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
            self.bind_object.set_error('biomart_type should be one of type1, type2, type3')

        biomart_info = dict()
        pep2transcript = dict()
        with open(biomart_file) as f:
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
                if strand_tmp == '1':
                    strand = '+'
                elif strand_tmp == '-1':
                    strand = '-'
                elif strand_tmp == '0':
                    strand = '.'
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
            self.bind_object.set_error('biomart information is None')
        else:
            # raw check
            if (not start.isdigit()) or (not end.isdigit()):
                self.bind_object.set_error('value of start and end is not digit')
        self.bind_object.logger.info('information of {} genes were parsed from biomart file'.format(len(biomart_info)))
        return biomart_info, pep2transcript

    def get_cds_seq(self, cds_path):
        '''
        Get cds information of transcript from provided FASTA file containing nucleotide sequence about cds.
        :param cds_path: absolute path of FASTA file about cds
        :return: dict, {trans_id：{"name": cds_id, "sequence": cds_sequence, "sequence_length": len(cds_sequence)}}
        '''
        cds_dict = dict()
        cds_pattern_match = re.compile(r'>([^\s]+)').match
        with open(cds_path) as f:
            j = 0
            trans_id = ''
            cds_id = ''
            cds_sequence = ''
            for line in f:
                if not line.strip():
                    continue
                if line.startswith('>'):
                    j += 1
                    if j > 1:
                        seq_len = len(cds_sequence)
                        cds_dict[trans_id] = dict(name=cds_id, sequence=cds_sequence, sequence_length=seq_len)
                        cds_dict[cds_id] = dict(name=cds_id, sequence=cds_sequence, sequence_length=seq_len)
                    cds_id = cds_pattern_match(line).group(1)
                    if '.' in cds_id:
                        trans_id = cds_id[:cds_id.rfind('.')]
                    else:
                        trans_id = cds_id
                    # cds_id and trans_id will be both saved as gtf may use either one of them
                    cds_sequence = ''
                else:
                    cds_sequence += line.strip()
            else:
                seq_len = len(cds_sequence)
                cds_dict[trans_id] = dict(name=cds_id, sequence=cds_sequence, sequence_length=seq_len)
                cds_dict[cds_id] = dict(name=cds_id, sequence=cds_sequence, sequence_length=seq_len)
        if not cds_dict:
            self.bind_object.set_error('cds sequence information is None')
        self.bind_object.logger.info('information of {} cds from transcript were parsed from FASTA file about cds'.format(len(cds_dict)))
        return cds_dict

    def get_pep_seq(self, pep_path, p2t):
        '''
        Get pep information of transcript from provided FASTA file containing protein sequence about pep.
        :param pep_path: absolute path of FASTA file about pep
        :param p2t: dict, {pep_id: transcript_id}
        :return: dict, {trans_id: {"name": pep_id, "sequence": pep_sequence, "sequence_length": len(pep_sequence)}}
        '''
        pep_dict = dict()
        j, trans_id, trans_id_else, pep_sequence, pep_id = 0, '', '', '', ''
        pep_pattern = re.compile(r'>([^\s]+)')
        trans_pattern = re.compile(r'transcript:([^\s]+)')
        with open(pep_path) as f:
            for line in f:
                if not line.strip():
                    continue
                if line.startswith('>'):
                    j += 1
                    if j > 1:
                        seq_len = len(pep_sequence)
                        if trans_id:
                            pep_dict[trans_id] = dict(name=pep_id, sequence=pep_sequence, sequence_length=seq_len)
                            pep_dict[trans_id_else] = dict(name=pep_id, sequence=pep_sequence, sequence_length=seq_len)
                    pep_id = pep_pattern.match(line).group(1)
                    try:
                        trans_id = trans_pattern.search(line).group(1)
                    except Exception:
                        if pep_id not in p2t:
                            if '.' in pep_id:
                                pep_id_else = pep_id[:pep_id.rfind('.')]
                                if pep_id_else not in p2t:
                                    self.bind_object.logger.debug('transcript id -> protein {} failed in biomart'.format(pep_id_else))
                                    trans_id = None
                                    continue
                                else:
                                    trans_id = p2t[pep_id_else]
                            else:
                                self.bind_object.logger.debug('transcript id -> protein {} failed in biomart'.format(pep_id))
                                trans_id = None
                                continue
                        else:
                            trans_id = p2t[pep_id]
                    if '.' in trans_id:
                        trans_id_else = trans_id[:trans_id.rfind('.')]
                    else:
                        trans_id_else = trans_id
                    pep_sequence = ''
                else:
                    pep_sequence += line.strip()
            else:
                seq_len = len(pep_sequence)
                pep_dict[trans_id] = dict(name=pep_id, sequence=pep_sequence, sequence_length=seq_len)
                pep_dict[trans_id_else] = dict(name=pep_id, sequence=pep_sequence, sequence_length=seq_len)
        if not pep_dict:
            self.bind_object.set_error('pep sequence information is None')
        self.bind_object.logger.info('information of {} pep from transcript were parsed from FASTA file about pep'.format(len(pep_dict)))
        return pep_dict

    def build_seq_db(self, db_path, gene_seq, trans_seq, t2g_dict, known_cds, known_pep):
        conn = sqlite3.connect(db_path)
        c = conn.cursor()
        # create gene_seq and trans_seq table
        seq_dicts = {'gene_seq': gene_seq, 'trans_seq': trans_seq}
        for table_name, seq_dict in seq_dicts.iteritems():
            c.execute('DROP TABLE IF EXISTS {}'.format(table_name))
            c.execute('CREATE TABLE {} (seq_id text, sequence text)'.format(table_name))
            for seq_id, seq in seq_dict.items():
                c.execute("INSERT INTO {} VALUES ('{}', '{}')".format(table_name, seq_id, seq))
        # create trans_annot table
        table_name = 'trans_annot'
        c.execute('DROP TABLE IF EXISTS {}'.format(table_name))
        c.execute('CREATE TABLE {} (transcript_id text, cds_id text, pep_id text, cds_seq text, pep_seq text, orf_type text)'.format(table_name))
        for transcript in t2g_dict:
            if transcript in known_cds:
                cds_id = known_cds[transcript]['name']
                cds_seq = known_cds[transcript]['sequence']
                orf_type = 'complete'
                if transcript in known_pep:
                    pep_id = known_pep[transcript]['name']
                    pep_seq = known_pep[transcript]['sequence']
                else:
                    pep_id = 'None'
                    pep_seq = 'None'
                c.execute("INSERT INTO {} VALUES ('{}', '{}', '{}', '{}', '{}', '{}')".format(table_name, transcript, cds_id, pep_id, cds_seq, pep_seq, orf_type))
            else:
                self.bind_object.logger.info('{} has no CDS and PEP annotation.'.format(transcript))
        # commit changes and return
        conn.commit()
        conn.close()
        return True

    def gene_location(self, bed_path):
        '''
        Get gene location from provided BED file, this file should have at least 6 columns.
        :param bed_path: absolute path of file with BED format
        :return: dict, {gene: {'chr': chr, 'strand': strand, 'start': start, 'end': end}}
        '''
        gene_info = dict()
        with open(bed_path) as f:
            for line in f:
                if not line.strip():
                    continue
                line = line.strip('\n').split("\t")
                if gene_info.get(line[3]) is None:
                    gene_info[line[3]] = {"chr": line[0], "strand": line[5], "start": str(int(line[1]) + 1), "end": line[2]}
                else:
                    self.bind_object.set_error('{} was found in multi-line at {}'.format(line[3], bed_path))
        return gene_info

    def dump_mongo(self, task_id, project_sn, db_path, t2g_dict, trans_loc, transcript_bed, trans_seq,
                   gene_loc, gene_bed, gene_seq, biomart_gene_detail, species_urls):
        # create main table
        time_now = datetime.datetime.now()
        main_info = dict([
            # essential keys
            ('task_id', task_id),
            ('project_sn', project_sn),
            ('created_ts', time_now.strftime('%Y-%m-%d %H:%M:%S')),
            ('name', 'GeneDetail_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))),
            ('desc', 'gene_detail_main_table'),
            ('params', None),
            # alternative keys
            ('seqdb', db_path),
            # status of process
            ('status', 'start'),
        ])
        main_id = self.create_db_table('sg_genes', [main_info])
        # add detail table
        gene_detail_list = list()
        g2t = dict()
        for t, g in t2g_dict.items():
            g2t.setdefault(g, set())
            g2t[g].add(t)
        total_trans_length = 0
        for g, t_set in g2t.items():
            trans_info = list()
            for t in t_set:
                if t not in trans_loc:
                    self.bind_object.set_error('{} was not found in {}'.format(t, transcript_bed))
                t_length = len(trans_seq[t])
                total_trans_length += t_length
                trans_info.append({
                    'transcript_id': t,
                    'start': trans_loc[t]['start'],
                    'end': trans_loc[t]['end'],
                    'length': t_length,
                })
            if g not in gene_loc:
                self.bind_object.set_error('{} was not found in {}'.format(g, gene_bed))
            start = gene_loc[g]['start']
            end = gene_loc[g]['end']
            location = '{}-{}'.format(start, end)
            length = len(gene_seq[g])
            if g in biomart_gene_detail:
                gene_type = biomart_gene_detail[g]['gene_type'][0]
            else:
                gene_type = 'predicted_non-coding'
            if 'ensembl' in species_urls:
                ensembl = "{}/Gene/Summary?g={}".format(species_urls[:-11], g)
            else:
                ensembl = None
            gene_detail_list.append({
                'gene_id': g,
                'transcripts': trans_info,
                'start': start,
                'end': end,
                'strand': gene_loc[g]['strand'],
                'chrom': gene_loc[g]['chr'],
                'location': location,
                'gene_length': length,
                'gene_type': gene_type,
                'gene_ensembl': ensembl,
                'main_id': main_id,
                'transcript_num': len(trans_info),
            })
        self.create_db_table('sg_genes_detail', gene_detail_list)
        self.update_db_record('sg_genes', main_id, status='end', main_id=main_id, transcripts_total_length=total_trans_length)
        return True

    def add_seq_stat(self, gene_stat):
        # task_id = self.bind_object.sheet.id
        # project_sn = self.bind_object.sheet.project_sn
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')

        # 导入基因序列统计
        name = 'Seq_gene_stat_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        main_info = dict([
            ('task_id', task_id),
            ('project_sn', project_sn),
            ('desc', 'seq_gene_stat'),
            ('created_ts', created_ts),
            ('status', 'start'),
            ('level', 'G'),
            ('name', name),
        ])
        gene_stat_main_id = self.create_db_table('sg_seq_stat', [main_info])
        gene_seq_stat_df = pd.read_table(gene_stat)
        gene_seq_stat_df["seq_stat_id"] = gene_stat_main_id
        gene_seq_stats = gene_seq_stat_df.to_dict("records")
        self.create_db_table('sg_seq_stat_detail', gene_seq_stats)
        self.update_db_record('sg_seq_stat', gene_stat_main_id, status='end', main_id=gene_stat_main_id)
        self.bind_object.logger.info("导入seq_stat成功!")
        # 导入转录本序列统计


class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.small_rna.small_rna_test_api import SmallRnaTestApiWorkflow
        from biocluster.wsheet import Sheet

        data = {
            'id': 'add_gene_detail_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'small_rna.small_rna_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = SmallRnaTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('small_rna.add_gene_detail')

        db_path = os.path.join(wf.output_dir, 'seqs.db')
        gene_fasta = '/mnt/ilustre/users/sanger-dev/workspace/20181211/Refrna_tsg_33013/GeneFa/output/gene.fa'
        transcript_fasta = '/mnt/ilustre/users/sanger-dev/workspace/20181211/Refrna_tsg_33013/RefrnaAssemble/output/NewTranscripts/all_transcripts.fa'
        t2g = '/mnt/ilustre/users/sanger-dev/workspace/20181211/Refrna_tsg_33013/RefrnaAssemble/output/NewTranscripts/trans2gene'
        biomart_file = '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/biomart/Mus_musculus.GRCm38.biomart_gene.txt'
        biomart_type = 'type1'
        cds_fasta = '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/cds/Mus_musculus.GRCm38.cds.all.fa'
        pep_fasta = '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/cds/Mus_musculus.GRCm38.pep.all.fa'
        gene_bed = '/mnt/ilustre/users/sanger-dev/workspace/20181211/Refrna_tsg_33013/GeneFa/ref_new_bed'
        transcript_bed = '/mnt/ilustre/users/sanger-dev/workspace/20181211/Refrna_tsg_33013/GeneFa/ref_new_trans_bed'
        species_urls = 'http://www.ensembl.org/Mus_musculus/Info/Index'

        wf.test_api.add_gene_detail(db_path, gene_fasta, transcript_fasta, t2g, biomart_file, biomart_type,
                        cds_fasta, pep_fasta, gene_bed, transcript_bed, species_urls)

if __name__ == '__main__':
    unittest.main()

#!/usr/bin/python
# -*- coding: utf-8 -*-
import re
import os
import time
import datetime
import types
import unittest
import sqlite3
from collections import OrderedDict, defaultdict
from bson.objectid import ObjectId
from bson.son import SON
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config


class RefrnaGeneDetail(Base):
    def __init__(self, bind_object):
        super(RefrnaGeneDetail, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + '_ref_rna'

    @staticmethod
    def biomart(biomart_path, biomart_type='type1'):
        """
        为了获得已知基因的description, gene_type信息
        :param biomart_path: this file must be tab separated.
        :param biomart_type: the type of biomart file
        :return: dict1. gene_id:  {"trans_id": [trans_id], "gene_name": [gene_name],
        "chromosome": [chromosome], "gene_type": [gene_type], "description": [desc],
         "strand": [strand], "pep_id": [pep_id], "start": [start], "end": [end]}
        dict2. pep->transcript
        """
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
            raise ValueError('biomart_type should be one of type1, type2, type3')

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
            raise Exception("biomart information is None")
        # raw check
        if (not start.isdigit()) or (not end.isdigit()):
            raise NameError('we find "start" or "end" is not digit. Maybe biomart_type is wrong')
        print('Information of {} genes was parsed from biomart file'.format(len(biomart_info)))
        return biomart_info, pep2transcript

    @staticmethod
    def get_cds_seq(cds_path):
        """
        从已经下载好的cds序列文件中提取转录本对应的cds序列。
        :param cds_path: cds序列文件的绝对路径
        :return: dict，转录本id：{"name": cds_id, "sequence": cds_sequence,
                                "sequence_length": len(cds_sequence)}
        """
        cds_dict = dict()
        cds_pattern_match = re.compile(r'>([^\s]+)').match
        with open(cds_path, 'r') as f:
            j = 0
            trans_id, cds_id, cds_sequence = '', '', ''
            for line in f:
                if not line.strip():
                    continue
                if line.startswith('>'):
                    j += 1
                    if j > 1:
                        seq_len = len(cds_sequence)
                        cds_dict[trans_id] = dict(name=cds_id, sequence=cds_sequence,
                                                  sequence_length=seq_len)
                        cds_dict[cds_id] = dict(name=cds_id, sequence=cds_sequence,
                                                sequence_length=seq_len)
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
                cds_dict[trans_id] = dict(name=cds_id, sequence=cds_sequence,
                                          sequence_length=seq_len)
                cds_dict[cds_id] = dict(name=cds_id, sequence=cds_sequence,
                                        sequence_length=seq_len)
        if not cds_dict:
            print('提取cds序列信息为空')
        print("共统计出{}条转录本的cds信息".format(len(cds_dict)))
        return cds_dict

    @staticmethod
    def get_pep_seq(pep_path, p2t):
        """
        get transcript's pep info, including protein sequence
        :param pep_path:
        :param p2t: dict of pep_id:transcript_id
        :return: dict, trans_id={"name": pep_id, "sequence": pep_sequence,
                                 "sequence_length": len(pep_sequence)}
        """
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
                            pep_dict[trans_id] = dict(name=pep_id, sequence=pep_sequence,
                                                      sequence_length=seq_len)
                            pep_dict[trans_id_else] = dict(name=pep_id, sequence=pep_sequence,
                                                           sequence_length=seq_len)
                    pep_id = pep_pattern.match(line).group(1)
                    try:
                        trans_id = trans_pattern.search(line).group(1)
                    except Exception:
                        if pep_id not in p2t:
                            print('transcript id -> protein {} failed in biomart'.format(pep_id))
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
                pep_dict[trans_id] = dict(name=pep_id, sequence=pep_sequence,
                                          sequence_length=seq_len)
                pep_dict[trans_id_else] = dict(name=pep_id, sequence=pep_sequence,
                                               sequence_length=seq_len)
        if not pep_dict:
            print('提取蛋白序列信息为空')
        print("共统计出{}条转录本的蛋白序列信息".format(len(pep_dict)))
        return pep_dict

    @staticmethod
    def parse_biomart_enterz(biromart_entrez):
        """
        parse *biomart_enterz.txt into dict
        :param biromart_entrez: example, no header, tab as separator
        -------------------------------------------------------
        ENSMUSG00000064341      ENSMUST00000082392      17716
        ENSMUSG00000064342      ENSMUST00000082393
        ENSMUSG00000064345      ENSMUST00000082396      17717
        -------------------------------------------------------
        :return: {gene_id: set(entrez_id1, en2)}
        """
        ensembl2entrez = dict()
        with open(biromart_entrez) as f:
            for line in f:
                if not line.strip():
                    continue
                tmp_list = line.strip().split('\t')
                if len(tmp_list) >= 3:
                    ensembl2entrez.setdefault(tmp_list[0], set())
                    ensembl2entrez[tmp_list[0]].add(tmp_list[2])
        return ensembl2entrez

    @staticmethod
    def parse_gene2ensembl(gene2ensembl):
        """
        parse file gene2ensembl into dict, contain information across species
        :param gene2ensembl: column info
        1	#tax_id
        2	GeneID
        3	Ensembl_gene_identifier
        4	RNA_nucleotide_accession.version
        5	Ensembl_rna_identifier
        6	protein_accession.version
        7	Ensembl_protein_identifier
        :return: dict, {prot: set{entrez_id1,id2,..}}
        """
        prot2entrez = dict()
        with open(gene2ensembl) as f:
            _ = f.readline()
            for line in f:
                if not line.strip():
                    continue
                tmp_list = line.strip().split('\t')
                if not tmp_list[5] == '-':
                    prot2entrez.setdefault(tmp_list[5], set())
                    prot2entrez[tmp_list[5]].add(tmp_list[1])
        return prot2entrez

    def query_pep_from_blast_result(self, query_id=None, blast_xls=None, blast_id=None):
        """根据新基因nr比对的结果 找出蛋白的accession号,方便后续根据蛋白号找到entrez id
        :param query_id: 待查询的新转录本id，if blast_xls used, this param will not be used
        :param blast_id: sg_annotaion_blast_detail中的blast_id，是ObjectId对象
        :param blast_xls: the blast result table. if used, other params will be not used.
        :return: if blast_xls, a dict with (pep_accession, description) as value will be returned.
                 if blast_id,  (pep_accession, description) will be returned.
        """
        if blast_xls:
            blast2pep = dict()
            with open(blast_xls) as f:
                header = f.readline().strip('\n').split('\t')
                query_ind = header.index("Query-Name")
                hit_name_ind = header.index("Hit-Name")
                desc_ind = header.index("Hit-Description")
                for line in f:
                    if not line.strip():
                        continue
                    tmp_list = line.strip('\n').split('\t')
                    query = tmp_list[query_ind]
                    if not blast2pep.get(query):
                        # Only the first hit_name will be record.
                        hit_name = tmp_list[hit_name_ind].strip('|').split('|')[-1]
                        description = tmp_list[desc_ind]
                        blast2pep[query] = (hit_name, description)
            return blast2pep
        elif blast_id:
            if not isinstance(blast_id, ObjectId):
                if isinstance(blast_id, types.StringTypes):
                    blast_id = ObjectId(blast_id)
                else:
                    raise Exception('blast_id 必须为ObjectId对象或其对应的字符串！')
            collection = self.db["sg_annotation_blast_detail"]
            data = collection.find_one({"blast_id": blast_id, "database": "nr",
                                        "anno_type": "gene", 'query_id': query_id})
            if data:
                hit_name = data["hit_name"]
                pep_accession = hit_name.strip('|').split('|')[-1]
                description = data['description']
                return pep_accession, description
            else:
                return '-', '-'
        else:
            raise Exception('blast_id or blast_xls must be specified')

    @staticmethod
    def gene_location(gene_bed_path):
        """
        Get gene location from provided bed file, this bed should have at least 6 columns.
        :param gene_bed_path: absolute path of novel gene file.
        :return: dict, gene:{chr, strand, start, end}
        """
        gene_info = dict()
        with open(gene_bed_path, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                line = line.strip('\n').split("\t")
                if gene_info.get(line[3]) is None:
                    gene_info[line[3]] = {"chr": line[0], "strand": line[5],
                                          "start": str(int(line[1])+1), "end": line[2]}
                else:
                    raise Exception("{}中基因{}在多行出现".format(gene_bed_path, line[3]))
        return gene_info

    @staticmethod
    def parse_class_code_info(class_code):
        """
        根据输入文件提取基因和转录本对应的关系，gene-id和gene name关系，转录本对应classcode关系。
        :param class_code: 文件名，tab分割.
        First line will be skipped. Example：
        -------------------------------------------------------------------------------------
        #assemble_txpt_id       assemble_gene_id        class_code      ref_gene_name
        ENSMUST00000166088      ENSMUSG00000009070      =       Rsph14
        ENSMUST00000192671      ENSMUSG00000038599      =       Capn8
        ENSMUST00000110582      ENSMUSG00000079083      =       Jrkl
        --------------------------------------------------------------------------------------
        :return: 3 dict, gene2trans_dict, gene2name_dict, trans2class_code
        """
        gene2trans_dict = dict()
        gene2name_dict = dict()
        trans2class_code = dict()
        with open(class_code, 'r') as f1:
            f1.readline()
            for line in f1:
                if not line.strip():
                    continue
                line = line.strip('\n').split("\t")
                gene = line[1]
                gene2trans_dict.setdefault(gene, list())
                gene2trans_dict[gene].append(line[0])
                trans2class_code[line[0]] = line[2]
                if len(line) >= 4:
                    gene2name_dict[gene] = line[3]
        return gene2trans_dict, gene2name_dict, trans2class_code

    @staticmethod
    def fasta2dict(fasta_file):
        """
        提取gene和transcript序列信息，都会用这个函数
        :param fasta_file:
        :return: dict, ｛seq_id: sequence｝
        """
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
                    if '(' in seq_id:
                        seq_id = seq_id.split('(')[0]
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
            print '提取序列信息为空'
        print "从{}共统计出{}条序列信息".format(fasta_file, len(seq))
        return seq

    @staticmethod
    def build_seq_database(seq_dicts, db_path):
        """
        :param seq_dicts: {table_name, seq_dict}
        :param db_path: abs path of seq db to be build.
        :return:
        """
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        # Create table

        for table_name in seq_dicts:
            cursor.execute('DROP TABLE IF EXISTS {}'.format(table_name))
            cursor.execute('CREATE TABLE {} (seq_id text, sequence text)'.format(table_name))
            seq_dict = seq_dicts[table_name]
            for seq_id in seq_dict:
                seq = seq_dict[seq_id]  # seq_id is transcript_id or gene_id
                if type(seq) == dict:
                    seq = seq_dict[seq_id]['sequence']
                cursor.execute("INSERT INTO {} VALUES ('{}', '{}')".format(table_name, seq_id, seq))
        conn.commit()
        conn.close()

    def add_gene_detail_class_code_detail(self, class_code,
                                          biomart_path=None,
                                          biomart_type="type1",
                                          biomart_entrez_path=None,
                                          gene2ensembl_path=None,
                                          blast_id='',
                                          blast_xls='',
                                          gene_location_path=None,
                                          trans_location_path=None,
                                          cds_path=None,
                                          pep_path=None,
                                          transcript_path=None,
                                          gene_path=None,
                                          species=None,
                                          assembly_method='StringTie',
                                          test_this=False):
        """
        :param class_code: class_code文件. "transcript--gene--class_code--gene_name"
        :param biomart_path: 获取已知基因的description和gene_type信息
        :param biomart_type: biomart_type, which decide the way to parse biomart_path
        :param biomart_entrez_path: info from ensembl. For query of a known gene's entrez id.
        :param gene2ensembl_path: gene2ensembl, info from ncbi. For query of a new gene's entrez id.
        :param blast_id: nr统计表 sg_annotation_blast_detail的id. if blast_xls used, it could be None.
        :param blast_xls: Result table file of NR blast for novel genes. Preferred than blast_id.
        :param gene_location_path: 基因的bed文件, 计算产生的文件
        :param trans_location_path: 转录本的bed文件，计算产生的文件
        :param cds_path: cds文件路径， from biomart
        :param pep_path: pep文件路径， from biomart
        :param transcript_path: 转录本的fa文件, 计算产生的文件
        :param gene_path: 基因的fa文件， 计算产生的文件
        :param species: ensembl的物种主页，其包含物种名称(拉丁文)
        :param assembly_method: assemble method, StringTie or Cufflink
        :param test_this: used for test purpose， and task_id will be "demo_test"
        **: liubinxu provide files: biomart_path, biomart_entrez_path, gene2ensembl_path,
                                     cds_path, pep_path.
            zengjing provide files: blast_xls
        """
        print("Begin time: {}".format(''.join(datetime.datetime.now().ctime())))
        start_time = time.time()
        local_args = locals()
        for arg in local_args:
            if local_args[arg] is None:
                raise Exception('{} must be specified'.format(arg))

        # -----create main table-----
        if not test_this:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
        else:
            project_sn = task_id = 'demo_test'
        all_transcript_length = 0
        with open(transcript_path) as f1:
            for line in f1:
                if line.startswith('>') or line.startswith('#'):
                    continue
                else:
                    all_transcript_length += len(line.strip())
        create_time = datetime.datetime.now()
        data = [('task_id', task_id),
                ('project_sn', project_sn),
                ('assembly_method', assembly_method),
                ('desc', 'gene_detail_information'),
                ('created_ts', create_time.strftime('%Y-%m-%d %H:%M:%S')), ('status', 'end'),
                ('name', 'Gene_detail_' + create_time.strftime("%Y%m%d_%H%M%S")),
                ('transcripts_total_length', all_transcript_length)]
        if not test_this:
            data.append(('refrna_seqdb', self.bind_object.work_dir + '/refrna_seqs.db'))
        else:
            data.append(('refrna_seqdb', os.path.join(os.getcwd(), 'refrna_seqs.db')))
        collection = self.db["sg_express_class_code"]
        class_code_id = collection.insert_one(SON(data)).inserted_id
        print("导入Gene_detail主表信息完成！")
        # -----end of creating main table-----
        # -------------------add class code information------------------------
        class_code_info = list()
        with open(class_code) as f:
            for line in f:
                line = line.strip('\n').split("\t")
                tmp_data = [('assembly_trans_id', line[0]),
                            ('assembly_gene_id', line[1]),
                            ('class_code', line[2]),
                            ('gene_name', line[3]),
                            ('class_code_id', class_code_id),
                            ('type', "express_diff")]
                tmp_data = SON(tmp_data)
                class_code_info.append(tmp_data)
        try:
            collection = self.db["sg_express_class_code_detail"]
            collection.insert_many(class_code_info)
            print("导入class_code_information完成")
        except Exception as e:
            print("导入class_code_information出错:{}".format(e))
        # -------------------------------------------
        # ----------------add gene detail----------------------------------------
        gene2trans, gene2name, trans2class_code = self.parse_class_code_info(class_code)
        biomart_data, pep2transcript = self.biomart(biomart_path, biomart_type=biomart_type)
        gene_location_info = self.gene_location(gene_location_path)
        trans_location_info = self.gene_location(trans_location_path)
        trans_cds_info = self.get_cds_seq(cds_path)
        trans_pep_info = self.get_pep_seq(pep_path, pep2transcript)
        trans_sequence = self.fasta2dict(transcript_path)
        gene_sequence_dict = self.fasta2dict(gene_path)
        prot2entrez = self.parse_gene2ensembl(gene2ensembl_path)
        ensembl2entrez = self.parse_biomart_enterz(biomart_entrez_path)
        if blast_xls:
            blast2pep = self.query_pep_from_blast_result(blast_xls=blast_xls)

            def query_pep(query_id=None, **useless):
                result = blast2pep.get(query_id)
                if result:
                    return result
                else:
                    return '-', '-'
        elif blast_id:
            query_pep = self.query_pep_from_blast_result
        else:
            raise Exception('blast_id or blast_xls must be specified')

        # ----------build seq database----------
        seq_dicts = dict(gene=gene_sequence_dict, transcript=trans_sequence,
                         pep=trans_pep_info, cds=trans_cds_info)
        if not test_this:
            db_path = self.bind_object.work_dir + '/refrna_seqs.db'
        else:
            db_path = os.path.join(os.getcwd(), 'refrna_seqs.db')

        self.build_seq_database(seq_dicts, db_path)
        # -------------------------------------

        data_list = list()
        new_num = 0
        for gene_id in gene2trans:
            description, gene_type, gene_name, is_new = '-', '-', '-', False

            # get gene's entrez id
            trans_list = gene2trans[gene_id]
            u_appear = [1 for x in trans_list if trans2class_code[x] == "u"]
            entrez_ids = '-'
            if sum(u_appear) < 1:
                if gene_id.startswith('MSTRG') or gene_id.startswith('TCONS') or gene_id.startswith('XLOC'):
                    print("{} looks like a novel gene, but judged as a known gene".format(gene_id))
                if ensembl2entrez.get(gene_id):
                    entrez_ids = ensembl2entrez[gene_id]
            else:
                is_new = True
                new_num += 1
                if (not gene_id.startswith('MSTRG')) and (not gene_id.startswith('TCONS')) and (not gene_id.startswith('XLOC')):
                    print("{} looks like a known gene, but judged as a novel one".format(gene_id))
                pep_id, description = query_pep(query_id=gene_id, blast_id=blast_id)
                if prot2entrez.get(pep_id):
                    # print(pep_id)
                    entrez_ids = prot2entrez[pep_id]

            # get transcript info
            trans_list = gene2trans[gene_id]
            transcript_num = len(trans_list)
            if not trans_list:
                error_info = '{}基因没有对应的转录本信息'.format(gene_id)
                raise Exception(error_info)
            transcripts = ",".join(trans_list)
            transcript_info = OrderedDict()
            new_transcript_info = OrderedDict()
            for ind, trans_ll in enumerate(trans_list):
                # as "." is not allowed to be used as key in Mongodb. str(ind) will be used instead.
                new_ind = str(ind)
                if ('MSTRG' not in trans_ll) and ('TCONS' not in trans_ll) and ('XLOC' not in trans_ll):
                    transcript_info[new_ind] = dict()
                    # 已知转录本输入的是cds和pep信息, 其键值索引是转录本的ensembl编号.
                    # if trans_ll in trans_cds_info.keys():
                    #     transcript_info[new_ind]["cds"] = trans_cds_info[trans_ll]
                    # else:
                    #     transcript_info[new_ind]["cds"] = '-'
                    #
                    # if trans_ll in trans_pep_info.keys():
                    #     transcript_info[new_ind]["pep"] = trans_pep_info[trans_ll]
                    # else:
                    #     transcript_info[new_ind]["pep"] = '-'

                    if trans_ll in trans_sequence:
                        # transcript_info[new_ind]["sequence"] = trans_sequence[trans_ll]
                        transcript_info[new_ind]["length"] = len(trans_sequence[trans_ll])
                    else:
                        transcript_info[new_ind]["length"] = '-'
                        # transcript_info[new_ind]["sequence"] = '-'
                else:
                    if trans_ll in trans_location_info:
                        tmp_value = trans_location_info[trans_ll]
                        new_transcript_info[new_ind] = dict(start=tmp_value['start'],
                                                            end=tmp_value['end'])
                    else:
                        print('{} was not found in {}'.format(trans_ll, trans_location_path))
                        new_transcript_info[new_ind] = dict(start='-', end='-')

                    if trans_ll in trans_sequence:
                        # new_transcript_info[new_ind]['sequence'] = trans_sequence[trans_ll]
                        new_transcript_info[new_ind]['length'] = len(trans_sequence[trans_ll])
                    else:
                        print('{} was not found in {}'.format(trans_ll, transcript_path))
                        # new_transcript_info[new_ind]['sequence'] = '-'
                        new_transcript_info[new_ind]['length'] = '-'

            # get known gene description and gene name
            if not is_new:
                if gene_id in biomart_data:
                    # 获取已知基因的description，gene_type, gene_name信息
                    description = biomart_data[gene_id]["description"][0]
                    gene_type = biomart_data[gene_id]['gene_type'][0]
                    if len(gene2name[gene_id]) <= 1:
                        gene_name = biomart_data[gene_id]['gene_name'][0]
                    else:
                        gene_name = gene2name[gene_id]
                else:
                    print('Warning: known gene {} was not in {}'.format(gene_id, biomart_path))

            # get gene location info
            if gene_id in gene_location_info:
                start = gene_location_info[gene_id]['start']
                end = gene_location_info[gene_id]['end']
                strand = gene_location_info[gene_id]['strand']
                chrom = gene_location_info[gene_id]['chr']
            else:
                print('{} was not found in {}'.format(gene_id, gene_location_path))
                start, end, strand, chrom = '-', '-', '-', '-'

            # get sequence info
            if gene_id in gene_sequence_dict:
                gene_sequence = gene_sequence_dict[gene_id]
                gene_length = len(gene_sequence)
            else:
                print('{} was not found in {}'.format(gene_id, gene_path))
                gene_sequence, gene_length = '-', '-'

            # format http sites
            if entrez_ids and (not entrez_ids == "-"):
                ncbi = 'https://www.ncbi.nlm.nih.gov/gene/?term={}'.format(list(entrez_ids)[0])
            else:
                ncbi = None
            if not is_new:
                if 'ensembl' in species:
                    ensembl = "{}/Gene/Summary?g={}".format(species[:-11], gene_id)
                else:
                    ensembl = None
            else:
                ensembl = None
            data = [
                    ("is_new", is_new),
                    ("type", "gene_detail"),
                    ("class_code_id", class_code_id),
                    ("gene_id", gene_id),
                    ("entrez_id", ','.join(entrez_ids)),
                    ("description", description),
                    ("strand", strand),
                    ("start", start),
                    ("location", "{}-{}".format(str(start), str(end))),
                    ("end", end),
                    ("chrom", chrom),
                    ("gene_name", gene_name),
                    ("transcript", transcripts),
                    ("transcript_number", transcript_num),
                    ("gene_type", gene_type),
                    # ("gene_sequence", gene_sequence),
                    ("gene_length", gene_length),
                    ("gene_ncbi", ncbi),  # only display one site
                    ("gene_ensembl", ensembl),
                    ]
            if transcript_info:
                data.append(("trans_info", transcript_info))
            if new_transcript_info:
                data.append(("new_trans_info", new_transcript_info))
            data = SON(data)
            data_list.append(data)

        try:
            collection = self.db["sg_express_class_code_detail"]
            collection.insert_many(data_list)
            print("Gene_detail_log: {} genes were dumped into Mongodb. {} are novel genes".format(
                len(gene2name), new_num))
        except Exception as e:
            print("Dumping gene details to Mongodb failed. log:{}".format(e))
        else:
            end_time = time.time()
            duration = end_time - start_time
            m, s = divmod(duration, 60)
            h, m = divmod(m, 60)
            print('导入基因详情表耗时 {}h:{}m:{}s'.format(h, m, s))


class TestFunction(unittest.TestCase):
    class TmpRefrnaGeneDetail(RefrnaGeneDetail):
        def __init__(self):
            super(TmpRefrnaGeneDetail, self).__init__(None)
            #self._db_name = Config().MONGODB + '_ref_rna'
            #self._db = None
            #self._config = Config()

    def test(self):
        base_path = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/"
        biomart_path = base_path + "biomart/Mus_musculus.GRCm38.biomart_gene.txt"
        biomart_type = "type1"
        biomart_entrez_path = base_path+"NCBI/Mus_musculus.GRCm38.biomart_enterz.txt"
        pep_path = base_path + "cds/Mus_musculus.GRCm38.pep.all.fa"
        cds_path = base_path + "cds/Mus_musculus.GRCm38.cds.all.fa"
        gene2ensembl_path = "/mnt/ilustre/users/sanger-dev/app/database/refGenome/ncbi_gene2ensembl/gene2ensembl"
        gene_bed = '/mnt/ilustre/users/sanger-dev/workspace/20170724/Single_gene_fa_5/GeneFa/ref_new_bed'
        trans_bed = '/mnt/ilustre/users/sanger-dev/workspace/20170724/Single_gene_fa_5/GeneFa/ref_new_trans_bed'
        gene_path = "/mnt/ilustre/users/sanger-dev/workspace/20170724/Single_gene_fa_2/GeneFa/output/gene.fa"
        transcript_path = "/mnt/ilustre/users/sanger-dev/workspace/20170706/Single_rsem_stringtie_mouse_total_2/Express1/TranscriptAbstract/output/exons.fa"
        class_code_info = "/mnt/ilustre/users/sanger-dev/workspace/20170702/Single_rsem_stringtie_mouse_total_1/Express/MergeRsem/class_code"
        blast_xls = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/ref_rna/ref_anno/taxonomy/mouse/new/anno_stat/blast/nr.xls"

        a = self.TmpRefrnaGeneDetail()
        # blast_id = a.db["sg_annotation_blast"].find_one({"task_id": "demo_0711"})
        # blast_id = blast_id['_id']
        a.add_gene_detail_class_code_detail(class_code_info,
                                            assembly_method="StringTie",
                                            biomart_path=biomart_path,
                                            biomart_type=biomart_type,
                                            biomart_entrez_path=biomart_entrez_path,
                                            gene2ensembl_path=gene2ensembl_path,
                                            gene_location_path=gene_bed,
                                            trans_location_path=trans_bed,
                                            cds_path=cds_path,
                                            pep_path=pep_path,
                                            species='http://www.ensembl.org/Mus_musculus/Info/Index',
                                            transcript_path=transcript_path,
                                            gene_path=gene_path,
                                            # blast_id=blast_id,
                                            blast_xls=blast_xls,
                                            test_this=True)

if __name__ == '__main__':
    unittest.main()

# !/usr/bin/python
# -*- coding: utf-8 -*-
import re
from collections import defaultdict
import sqlite3
import datetime
from api_base import ApiBase
import unittest
from biocluster.config import Config

class AddGeneDetail(ApiBase):
    def __init__(self, bind_object):
        super(AddGeneDetail, self).__init__(bind_object)

    def add_gene_detail(self, db_path, gene_bed, transcript_bed, species_urls, biomart_file, biomart_type,
                        known_cds, known_pep,  new_cds, new_pep, transcript_fasta, gene_fasta, t2g, gene_type, trans_type, lnc_ids):

        print ", ".join([db_path, gene_bed, transcript_bed, species_urls, biomart_file, biomart_type,
                        known_cds, known_pep,  new_cds, new_pep, transcript_fasta, gene_fasta, t2g, gene_type, trans_type, lnc_ids])
        """
        :param db_path: 待创建的sqlite3数据库路径
        :param gene_bed: 基因bed文件
        :param transcript_bed: 转录本bed文件
        :param species_urls: 物种主页链接
        :param biomart_file: biomart注释文件
        :param biomart_type: type of biomart file
        :param known_cds: 已知cds文件路径，如下：
        /mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/cds/Mus_musculus.GRCm38.cds.all.fa
        :param known_pep: 已知pep文件路径，如下：
        /mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/cds/Mus_musculus.GRCm38.pep.all.fa
        :param new_cds: 新cds文件路径，如下：
        /mnt/ilustre/users/sanger-dev/workspace/20180531/Single_annot_new_mus_orf18-20-37/AnnotOrfpfam/output/Mus.new.exon.fa.transdecoder.cds
        :param new_pep: 新pep文件路径
        /mnt/ilustre/users/sanger-dev/workspace/20180531/Single_annot_new_mus_orf18-20-37/AnnotOrfpfam/output/Mus.new.exon.fa.transdecoder.pep
        :param transcript_fasta: 转录本序列文件路径
        :param gene_fasta: 基因序列文件路径
        :param t2g: 转录本和基因对应关系文件路径
        :param gene_type: 基因类型
        :param trans_type: 转录本类型
        :param lnc_ids: lncRNA列表，对应id关系
        :return:
        """
        gene_loc = self.gene_location(gene_bed)
        trans_loc = self.gene_location(transcript_bed)
        gene_type = self.get_gene_type(gene_type)
        t2g_dict = dict([(x.strip().split()[0], x.strip().split()[1]) for x in open(trans_type) if x.strip()])
        trans_type = self.get_trans_type(trans_type)
        lnc_ids, linc_id2link =self.get_lnc_ids(lnc_ids, species_urls)
        biomart_gene_detail, pep2trans = self.biomart(biomart_file, biomart_type)
        
        known_cds = self.get_cds_seq(known_cds)
        known_pep = self.get_pep_seq(known_pep, pep2trans)
        if new_cds and new_pep:
            new_cds_pep = self.get_predict_cds_pep(new_cds, new_pep)
        else:
            new_cds_pep = dict()
        gene_seq = self.fasta2dict(gene_fasta)
        trans_seq = self.fasta2dict(transcript_fasta)
        # ----------build seq database ---------------------------------------------------------
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        # Create gene_seq and trans_seq table
        seq_dicts = dict(gene_seq=gene_seq, trans_seq=trans_seq)
        for table_name, seq_dict in seq_dicts.items():
            cursor.execute('DROP TABLE IF EXISTS {}'.format(table_name))
            cursor.execute('CREATE TABLE {} (seq_id text, sequence text)'.format(table_name))
            for seq_id, seq in seq_dict.items():
                cursor.execute("INSERT INTO {} VALUES ('{}', '{}')".format(table_name, seq_id, seq))

        # create trans_annot table
        table_name = "trans_annot"
        cursor.execute('DROP TABLE IF EXISTS {}'.format(table_name))
        cursor.execute('CREATE TABLE {} (transcript_id text, cds_id text, pep_id text, cds_seq text, pep_seq text, orf_type text)'.format(table_name))
        for transcript in t2g_dict:
            if transcript in known_cds:
                cds_id = known_cds[transcript]['name']
                cds_seq = known_cds[transcript]['sequence']
                orf_type = 'complete'
                if transcript in known_pep:
                    pep_id = known_pep[transcript]['name']
                    pep_seq = known_pep[transcript]['sequence']
                    orf_type = 'complete'
                else:
                    print transcript
                    pep_id = 'None'
                    pep_seq = 'None'
                cursor.execute("INSERT INTO {} VALUES ('{}', '{}', '{}', '{}', '{}', '{}')".format(
                    table_name, transcript, cds_id, pep_id, cds_seq, pep_seq, orf_type))

            elif transcript in new_cds_pep:
                for seq_id, _type, cds_seq, pep_seq in new_cds_pep[transcript]:
                    cursor.execute("INSERT INTO {} VALUES ('{}', '{}', '{}', '{}', '{}', '{}')".format(
                        table_name, transcript, seq_id, seq_id, cds_seq, pep_seq, _type))
            else:
                print('{} has no CDS and PEP annotation.'.format(transcript))
        conn.commit()
        conn.close()
        # -------------------------------------------------------------------------------
        # ---------------dump gene location to Mongo ----
        # create main table
        try:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
        except Exception:
            project_sn = task_id = 'do_test'
        create_time = datetime.datetime.now()
        main_info = dict([
            ('task_id', task_id),
            ('project_sn', project_sn),
            ('desc', 'gene_detail_information'),
            ('created_ts', create_time.strftime('%Y-%m-%d %H:%M:%S')),
            ('status', 'start'),
            ('name', 'Gene_detail_' + create_time.strftime("%Y%m%d_%H%M%S")),
            ('refrna_seqdb', db_path)
        ])
        main_id = self.create_db_table('sg_genes', [main_info])

        # add detail
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
                    raise Exception("{} not found in {}".format(t, transcript_bed))
                t_length = len(trans_seq[t])
                total_trans_length += t_length
                trans_info.append(dict(
                    transcript_id=t,
                    start=trans_loc[t]['start'],
                    end=trans_loc[t]['end'],
                    length=t_length,
                ))
            if g not in gene_loc:
                raise Exception("{} not found in {}".format(g, gene_bed))
            start = gene_loc[g]['start']
            end = gene_loc[g]['end']
            location = start + '-' + end
            length = len(gene_seq[g])
            '''
            if g in biomart_gene_detail:
                gene_type = biomart_gene_detail[g]['gene_type'][0]
            else:
                gene_type = 'predicted_non-coding'
                for t in t_set:
                    if t in new_cds_pep:
                        gene_type = 'protein_coding'
                        break
            '''
            if g in linc_id2link:
                ensembl = linc_id2link[g]
            elif 'ensembl' in species_urls:
                ensembl = "{}/Gene/Summary?g={}".format(species_urls[:-11], g)
            else:
                ensembl = None

            gene_name= biomart_gene_detail[g]["gene_name"][0] if g in biomart_gene_detail else "",
            description=biomart_gene_detail[g]["description"][0] if g in biomart_gene_detail else "",
            gene_detail_list.append(dict(
                gene_id=g,
                transcripts=trans_info,
                start=start,
                end=end,
                strand=gene_loc[g]['strand'],
                chrom=gene_loc[g]['chr'],
                location=location,
                gene_length=length,
                gene_type=gene_type[g],
                gene_name=gene_name,
                description=description,
                gene_ensembl=ensembl,
                genes_id=main_id,
                level="G",
                rna_type = gene_type[g],
                transcript_num=len(trans_info),
            ))

        for t, g in t2g_dict.items():
            if t in linc_id2link:
                ensembl = linc_id2link[t]
            elif 'ensembl'  in species_urls:
                ensembl = "{}/Gene/Summary?g={}".format(species_urls[:-11], g)
            else:
                ensembl = None
            start=trans_loc[t]['start']
            end=trans_loc[t]['end']
            t_length = len(trans_seq[t])
            location = start + '-' + end
            length=t_length
            if t in lnc_ids:
                alias = lnc_ids[t]
            else:
                alias = dict()

            gene_name= biomart_gene_detail[g]["gene_name"][0] if g in biomart_gene_detail else "",
            description=biomart_gene_detail[g]["description"][0] if g in biomart_gene_detail else "",
            gene_detail_list.append(dict(
                gene_id=g,
                trans_id=t,
                start=start,
                end=end,
                strand=trans_loc[t]['strand'],
                chrom=trans_loc[t]['chr'],
                location=location,
                gene_length=length,
                gene_type=trans_type[t],
                gene_ensembl=ensembl,
                gene_name=gene_name,
                description=description,
                genes_id=main_id,
                level="T",
                rna_type = trans_type[t],
                alias = alias,
                exon_num = trans_loc[t]['exon_num'] if 'exon_num' in trans_loc[t] else 1,
            ))

        self.create_db_table('sg_genes_detail', gene_detail_list)
        self.update_db_record('sg_genes', main_id, status="end", main_id=main_id, transcripts_total_length=total_trans_length)



    def get_gene_type(self, gene_type):
        # 获取基因类型
        with open(gene_type, 'r') as f:
            gene2type = {line.strip().split("\t")[0]:line.strip().split("\t")[2] for line in f.readlines()}
        return gene2type

    def get_trans_type(self, trans_type):
        # 获取转录本类型
        with open(trans_type, 'r') as f:
            trans2type = {line.strip().split("\t")[0]:line.strip().split("\t")[2] for line in f.readlines()}
        return trans2type

    def get_lnc_ids(self, lnc_ids, species_urls):
        # 获取已知lncRNA数据库ID
        
        lnc_url = {
            "ENSEMBL" : '{}/Gene/Summary?g='.format(species_urls[:-11]) + '{}',
            "NCBI" : 'https://www.ncbi.nlm.nih.gov/search/all/?term={}',
            "NONCODE" : 'http://www.noncode.org/show_rna.php?id={}',
        }
        green_nc_ids = Config().SOFTWARE_DIR + '/database/lnc_rna/greenc/greenc_lncRNA_link.txt'
        with open(green_nc_ids, 'r') as f:
            green_dict = {line.strip().split()[0]: line.strip().split()[1] for line in f.readlines()}
        trans2link = dict()
        linc_id2link = dict()
        with open(lnc_ids, 'r') as f:
            heads = f.readline().strip().split("\t")
            for line in f:
                cols = line.strip().split("\t")
                link_dict = dict(zip(heads, cols))
                out_dict = dict()
                source = link_dict['source']
                for k, v in link_dict.items():
                    if v != "" and k.endswith("_transcript_id"):
                        db = k.split("_transcript_id")[0].upper()
                        if db in lnc_url:
                            if db != source.upper():
                                ids = []
                                links = []
                                for lncs in v.split(","):
                                    ids.append(lncs)
                                    links.append(lnc_url[db].format(lncs))
                                out_dict[db] = {"ids": ids, "links": links}
                        elif db == "GREENC":
                            if db != source.upper():
                                ids = []
                                links = []
                                for lncs in v.split(","):
                                    ids.append(lncs)
                                    links.append(green_dict[lncs])
                                out_dict[db] = {"ids": ids, "links": links}
                trans2link[cols[0]] = out_dict
                source = link_dict["source"]
                if source.upper() == "GREENC":
                    linc_id2link[cols[0]] = green_dict.get(cols[0], '')
                elif source.upper() in ["ENSEMBL", "NCBI", "NONCODE"]:
                    linc_id2link[cols[0]] = lnc_url[source.upper()].format(cols[0])
                else:
                    pass

        return trans2link, linc_id2link
    

    def gene_detail_example(self):
        # 基因详情表要存储的信息示例
        sg_gene_detail= {
            "gene_id": "gene123",
            "transcripts": [
                # ["t1", "start", "end", "length"],
                {"transcript_id": "transcript1", "start": 50, "end": 1050, "length":1050},
                {"transcript_id": "transcript2", "start": 56, "end": 1080, "length":1044},
                {"transcript_id": "transcript3", "start": 67, "end": 1070, "length":993},
            ],
            "chrom": 6,
            "start": 40,
            "end": 1200,
            "strand": '+',
            "gene_length": 1161,
            "gene_type": "protein_coding",
            "location":'40-1200',
            # "ncbi": "https://www.ncbi.nlm.nih.gov/gene/?term=50614", 需前端去sg_annotation_query表获得entrez_id后拼接
            "gene_ensemble": "http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000182870",
            # ensembl = "{}/Gene/Summary?g={}".format(species_urls, gene_id)
            "main_id": 'ObjectID(xxx)'
        }

    def seqdb_example(self):
        pass
        # seqdb 包含3张表，一张是转录本的注释表，一张是基因的序列表，还有一张是转录本的序列表
        """
        【transcript_annot 表结构】
        --------------------------------------------------------------------------------
        1. 对于新转录本，可能预测有几个cds（orf)，每条cds(orf)仅对应一个pep， 而且cds_id和pep_id相同，它包含了transcript_id信息。
        2. 对于已知转录本，cds_id就是transcript_id， pep_id是蛋白id， transcript_id和pep_id的对应关系是有biomart文件提供的。
        3. type: 对于已转录本：type是complete；对于新转录本的type有：5prime_partial/complete
        表头如下：    
        |transcript_id | cds_id | pep_id | cds_seq | pep_seq | type |
        --------------------------------------------------------------------------------
        
        【gene_seq 表结构】
        ------------------------------------------------------------
        | gene_id | gene_seq |
        ------------------------------------------------------------
        
        【trans_seq 表结构】
        ------------------------------------------------------------
        | transcript_id | transcript_seq |
        ------------------------------------------------------------
        """

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
        else:
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
                                    print('transcript id -> protein {} failed in biomart'.format(pep_id_else))
                                    trans_id = None
                                    continue
                                else:
                                    trans_id = p2t[pep_id_else]
                            else:
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
                pep_dict[trans_id] = dict(name=pep_id, sequence=pep_sequence, sequence_length=seq_len)
                pep_dict[trans_id_else] = dict(name=pep_id, sequence=pep_sequence, sequence_length=seq_len)
        if not pep_dict:
            print('提取蛋白序列信息为空')
        print("共统计出{}条转录本的蛋白序列信息".format(len(pep_dict)))
        return pep_dict

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
                    gene_info[line[3]] = {"chr": line[0], "strand": line[5], "start": str(int(line[1]) + 1),
                                          "end": line[2]}
                    if len(line) >= 10:
                        gene_info[line[3]]["exon_num"] = line[9]
                else:
                    raise Exception("{}中基因{}在多行出现".format(gene_bed_path, line[3]))
        return gene_info

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
                    # gene id和转录本id可能包含括号
                    # 取消特殊字符判断
                    if '(+)' in seq_id:
                        seq_id = seq_id.split('(+)')[0]
                    elif '(-)' in seq_id:
                        seq_id = seq_id.split('(-)')[0]
                    elif '(.)' in seq_id:
                        seq_id = seq_id.split('(.)')[0]
                    # elif '|' in seq_id:
                    #     seq_id = seq_id.split('|')[0]
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
        """
        generator for parsing sequence file
        :param seq_file: fasta sequence file
        :param match: regexp pattern for parsing line startswith '>'
        :return: seq_tag, sequence
        """
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

    def get_predict_cds_pep(self, cds_file, pep_file):
        # parse pep
        pep_parser = self.parse_seq_file(pep_file)
        pep_dict = dict()
        for pep_desc, pep_seq in pep_parser:
            pep_dict[pep_desc] = pep_seq
        # parse cds
        t2cds_pep = dict()
        match = re.compile(r'>(.*?)::(.*?)::.*type:(.*?)\s+len:\d+.*:(\d+-\d+\(.\)).*').match
        cds_parser = self.parse_seq_file(cds_file)
        for cds_desc, cds_seq in cds_parser:
            # ('TRINITY_DN1001_c0_g1', 'TRINITY_DN1001_c0_g1_i1', 'complete', '92-931(+)')
            g_id, t_id, _type, cds_pos = match(cds_desc).groups()
            pep_seq = pep_dict[cds_desc]
            seq_id = re.match(r'>(.*?)\s.*', cds_desc).groups()[0]
            t2cds_pep.setdefault(t_id, list())
            t2cds_pep[t_id].append((seq_id, _type, cds_seq, pep_seq))  # pep_id = cds_id
        return t2cds_pep

'''
    def run1(self):
        db_path = "/mnt/ilustre/users/sanger-dev/workspace/20180625/Refrna_tsg_30770/refrna_seqs.db"
        gene_bed = "/mnt/ilustre/users/sanger-dev/workspace/20180625/Refrna_tsg_30770/GeneFa/ref_new_bed"
        transcript_bed = "/mnt/ilustre/users/sanger-dev/workspace/20180625/Refrna_tsg_30770/GeneFa/ref_new_trans_bed"
        transcript_fasta = "/mnt/ilustre/users/sanger-dev/workspace/20180625/Refrna_tsg_30770/RefrnaAssemble/output/NewTranscripts/all_transcripts.fa"
        t2g = "/mnt/ilustre/users/sanger-dev/workspace/20180625/Refrna_tsg_30770/RefrnaAssemble/output/NewTranscripts/trans2gene"
        gene_fasta = "/mnt/ilustre/users/sanger-dev/workspace/20180625/Refrna_tsg_30770/GeneFa/output/gene.fa"
        species_urls = "http://www.ensembl.org/Mus_musculus/Info/Index"
        biomart_file = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/biomart/Mus_musculus.GRCm38.biomart_gene.txt"
        biomart_type = "type1"
        known_cds = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/cds/Mus_musculus.GRCm38.cds.all.fa"
        known_pep = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/cds/Mus_musculus.GRCm38.pep.all.fa"
        new_cds = "/mnt/ilustre/users/sanger-dev/workspace/20180625/Refrna_tsg_30770/AnnotOrfpfam/output/new_transcripts.fa.transdecoder.cds"
        new_pep = "/mnt/ilustre/users/sanger-dev/workspace/20180625/Refrna_tsg_30770/AnnotOrfpfam/output/new_transcripts.fa.transdecoder.pep"
        self.add_gene_detail(db_path, gene_bed, transcript_bed, species_urls, biomart_file, biomart_type, known_cds, known_pep, new_cds, new_pep, transcript_fasta, gene_fasta, t2g)

class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """
    from mbio.workflows.ref_rna_v2.refrna_test_api import RefrnaTestApiWorkflow
    from biocluster.wsheet import Sheet
    import random

    data = {
        "id": "tsg_30770",
        "project_sn": "188_5b03a60c4f528",
        "type": "workflow",
        "name": "ref_rna_v2.ref_test_api",
        "options": {
        },
    }
    wsheet = Sheet(data=data)
    wf = RefrnaTestApiWorkflow(wsheet)
    wf.IMPORT_REPORT_DATA = True
    wf.IMPORT_REPORT_AFTER_END = False
    wf.test_api = wf.api.api("ref_rna_v2.add_gene_detail")
    wf.test_api.run1()

if __name__ == '__main__':
    unittest.main()

'''

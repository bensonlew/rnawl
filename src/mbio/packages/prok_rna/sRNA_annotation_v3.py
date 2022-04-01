## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "yitong.feng"
# 20180724

import os
import re
import sys
from collections import defaultdict
from Bio.Blast import NCBIXML
from biocluster.config import Config
from mako.template import Template


class srna_annot(object):
    def __init__(self, predict_fa, evalue, rfam):
        self.predict_fa = predict_fa
        self.predix = os.path.basename(self.predict_fa)
        self.evalue = evalue
        self.rfam = rfam
        # self.blastn_path = Config().SOFTWARE_DIR + '/bioinfo/align/ncbi-blast-2.3.0+/bin/blastn'
        self.blastn_path = Config().SOFTWARE_DIR + '/bioinfo/ref_rna_v2/miniconda2/bin/blastn'

        self.db_path = Config().SOFTWARE_DIR + '/database/prok_rna/sRNA_Anno/'
        self.parafly = Config().SOFTWARE_DIR + '/program/parafly-r2013-01-21/src/ParaFly'

    def run_blast(self):
        blast_cmd = \
"""
${blastn} -task blastn-short -strand plus  -num_threads 5 -db ${db_path}/SIPHI_nc  -query ${predict_fa} -outfmt 5 -evalue ${evalue} -num_alignments 1 -out ${predix}_vs_SIPHI.xml
${blastn} -task blastn-short -strand plus  -num_threads 5 -db ${db_path}/sRNA_Gene_Seq  -query ${predict_fa} -outfmt 5 -evalue ${evalue} -num_alignments 1 -out ${predix}_vs_sRNAMap.xml
${blastn} -task blastn-short -strand plus  -num_threads 5 -db ${db_path}/sRNATarBase_V0002 -query ${predict_fa} -outfmt 5 -evalue ${evalue} -num_alignments 1 -out ${predix}_vs_sRNATarBase.xml
${blastn} -task blastn-short -strand plus  -num_threads 5 -db ${db_path}/BSRD  -query ${predict_fa} -outfmt 5 -evalue ${evalue} -num_alignments 1 -out ${predix}_vs_BSRD.xml
"""
        f = Template(blast_cmd)
        blast_cmd = f.render(blastn=self.blastn_path,
                             db_path=self.db_path,
                             predict_fa=self.predict_fa,
                             evalue=self.evalue,
                             predix=self.predix,
                             )
        with open('blast_cmd', 'w') as blast_l:
            blast_l.write(blast_cmd)
        os.system(self.parafly + " -c blast_cmd -CPU 40")

    def xml2table(self, xml_fp, table_out):
        """

            :param xml_fp: 输入xml文件路径
            :param table_out: 输出文件路径
            :param header: 选择列，列的选择必须是上面定义的 all_values中的值组成的列表
            :param anno_head: 是否写入表头
            """
        need_values = ['Score', 'E-Value', 'Hit-Len', 'Identity-%', 'Similarity-%',
                      'Query-Name', 'Q-Len', 'Q-Begin', 'Q-End', 'Hit-Name', 'Hit-Len', 'Hsp-Begin',
                      'Hsp-End', 'Hit-Description']
        header = 'Score\tE-Value\tHSP-Len\tIdentity(%)\tSimilarity(%)\tsRNA ID\tsRNA Len (nt)\tsRNA Start\tsRNA End\tHit Name\tHit Len (nt)\tHit Start\tHit End\tDescription\n'
        if not os.path.isfile(xml_fp):
            raise Exception('输入xml文件不存在:{}'.format(xml_fp))
        with open(xml_fp) as f, open(table_out, 'w') as w:
            w.write(header)
            records = NCBIXML.parse(f)
            values = {i: 'N/A' for i in need_values}
            for rec in records:
                query = re.split(' ', rec.query, maxsplit=1)[0]
                for align in rec.alignments:
                    for hsp in align.hsps:
                        one_hsp = values.copy()
                        one_hsp['Query-Name'] = query
                        hit_def = align.hit_def
                        if u'|' in hit_def:
                            one_hsp['Hit-Name'] = hit_def.split('|')[0]
                            one_hsp['Hit-Description'] = '|'.join(hit_def.split('|')[1:]).strip()
                        else:
                            one_hsp['Hit-Name'] = hit_def.split(' ')[0]
                            one_hsp['Hit-Description'] = ' '.join(hit_def.split(' ')[1:]).strip()
                        # one_hsp['Hit-Name'] = align.hit_def.split(' ')[0]
                        # one_hsp['Hit-Description'] = align.hit_def
                        one_hsp['Score'] = str(hsp.score)
                        one_hsp['E-Value'] = str(hsp.expect)
                        one_hsp['HSP-Len'] = str(hsp.align_length)
                        one_hsp['Identity'] = str(hsp.identities)
                        one_hsp['Positives'] = str(hsp.positives)
                        one_hsp['Q-Len'] = str(rec.query_length)
                        one_hsp['Q-Begin'] = str(hsp.query_start)
                        one_hsp['Q-End'] = str(hsp.query_end)
                        one_hsp['Q-Frame'] = str(hsp.frame[0])
                        one_hsp['Hit-Len'] = str(align.length)
                        one_hsp['Hsp-Begin'] = str(hsp.sbjct_start)
                        one_hsp['Hsp-End'] = str(hsp.sbjct_end)
                        one_hsp['Hsp-Frame'] = str(hsp.frame[1])
                        one_hsp['Q-Strand'] = str(hsp.strand[0])
                        one_hsp['Hsp-Strand'] = str(hsp.strand[1])
                        one_hsp['Identity-%'] = str(round(float(hsp.identities) / hsp.align_length, 3) * 100)
                        one_hsp['Similarity-%'] = str(round(float(hsp.positives) / hsp.align_length, 3) * 100)
                        one_hsp['Mismatch'] = str(
                            int(hsp.align_length) - int(len(hsp.match)))  # 此处使用len(hsp.match) 还是hsp.indentities没有具体证据
                        one_hsp['Gapopen_num'] = str(hsp.gaps)  # 此处的gaps不是 gapopen 而是 gaps总数，因为xml获取不到
                        line = list()
                        for i in need_values:
                            line.append(one_hsp[i])
                        w.write('\t'.join(line) + '\n')

    def get_list_file(self, infile, anno_type):
        with open(infile, 'r') as xls_r, open(infile.strip('xls') + 'list', 'w') as list_w:
            query_list = list()
            xls_r.readline()
            for line in xls_r:
                line = line.strip().split('\t')
                if anno_type != 'rfam' and len(line) >= 6:
                    query_list.append(line[5])
                elif anno_type == 'rfam' and len(line) >= 4:
                    if line[3] not in query_list:
                        query_list.append(line[3])
            list_w.write('\n'.join(query_list))

    def modify_rfam_description(self):
        rf2des = dict()
        with open(self.db_path + 'description.xls', 'r') as des_r:
            des_r.readline()
            for line in des_r:
                line = line.strip().split('\t')
                if len(line) > 3:
                    rf2des[line[0]] = line[3]
        rfam_r = self.rfam
        rfam_o = self.predix + '_vs_rfam.xls'
        with open(rfam_r, 'r') as r, open(rfam_o, 'w') as rfam_w:
            header = r.readline().strip().split(' ')
            header = [i for i in header if i != '']
            rfam_w.write('\t'.join(header) + '\n')
            r.readline()
            for line in r:
                info = line.strip().split(' ')
                info = [i for i in info if i != '']
                if len(info) > 3:
                    rf = info[2]
                    try:
                        info[-1] = rf2des[rf]
                    except:
                        info[-1] = '_'
                    rfam_w.write('\t'.join(info) + '\n')
                elif info[0] == '#':
                    break

    def merge_stat_rfam(self):
        rfam_info = defaultdict(dict)
        tmp = ''
        # with open(self.db_path + 'Rfam.seed', 'r') as rfam_seed:
        with open(self.db_path + 'Rfam-14.6/Rfam.seed', 'r') as rfam_seed:
            for line in rfam_seed:
                if line.startswith('#'):
                    if re.match('^#=GF\s+AC\s+(\S+)\s*$', line):
                        tmp = re.match('^#=GF\s+AC\s+(\S+)\s*$', line).group(1)
                    if re.match('^#=GF\s+DE\s+(.*)\s*$', line):
                        rfam_info[tmp]['DE'] = re.match('^#=GF\s+DE\s+(.*)\s*$', line).group(1)
                    if re.match('^#=GF\s+TP\s+(.*)\s*$', line):
                        rfam_info[tmp]['TP'] = re.match('^#=GF\s+TP\s+(.*)\s*$', line).group(1)

        with open(self.predict_fa, 'r') as fa:
            all_ = 0
            for line in fa:
                if line.startswith('>'):
                    all_ += 1

        with open(self.predix + '_vs_rfam.xls', 'r') as rfam_r, \
                open('rfam_merge', 'w') as merge_w, \
                open('rfam_stat.xls', 'w') as rfam_s:
            rfam_r.readline()
            TP_list = list()
            merge_w.write('query_name\tTP\tDE\tAC\tHit\tHspLength\tEvalue\tScore\n')
            for line in rfam_r:
                line = line.strip().split('\t')
                if len(line) >= 18:
                    query_name = line[3]
                    hit = line[1]
                    AC = line[2]
                    length = abs(int(line[9]) - int(line[10])) + 1
                    TP = DE = '-'
                    if AC in rfam_info:
                        TP = rfam_info[AC]['TP']
                        DE = rfam_info[AC]['DE']
                    if TP != '-':
                        TP_list.append(TP)
                    merge_w.write(query_name + '\t' + TP + '\t' + DE + '\t' + AC + '\t' + hit + '\t'
                                  + str(length) + '\t' + line[17] + '\t' + line[16] + '\n')

            num2DE = defaultdict(set)
            for i in TP_list:
                num2DE[TP_list.count(i)].add(i)
            rfam_s.write('All_reads' + '\t' + str(all_) + '\t' + '100%\n')
            rfam_s.write('All_matched' + '\t' + str(len(TP_list)) + '\t' +
                         str(round(float(len(TP_list))/all_*100, 2)) + '%\n')
            for i in sorted(num2DE.keys(), reverse=True):
                for j in num2DE[i]:
                    rfam_s.write(j + '\t' + str(i) + '\t' + str(round(float(i)/all_*100, 2)) + '%\n')

    def merge_annot(self, blast_xls, annot_file, merged_file):
        with open(blast_xls, 'r') as blast_r, \
            open(annot_file, 'r') as annot_r, \
            open(merged_file, 'w') as merge_w:
            annot_info = dict()
            header = 'Query-Name\tHit-Name\t' + '\t'.join(annot_r.readline().strip().split('\t')[2:])
            merge_w.write(header + '\n')
            for line in annot_r.readlines():
                line = line.strip().split('\t')
                annot_info[line[0]] = '\t'.join(line[1:])
            _ = blast_r.readline()
            for line in blast_r.readlines():
                line = line.strip().split('\t')
                if len(line) >= 10:
                    if line[9] in annot_info:
                        merge_w.write(line[5] + '\t' + line[9] + '\t' + annot_info[line[9]] + '\n')

    def stat_annot(self):
        rfam_s2hit = defaultdict(set)
        rfam_list = list()
        with open(self.predix + '_vs_rfam.xls', 'r') as rfam_r:
            rfam_r.readline()
            for line in rfam_r:
                line = line.strip().split('\t')
                if len(line) >= 4:
                    rfam_s2hit[line[3]].add(line[2])
                    rfam_list.append(line[3])

        SIPHI_s2hit = defaultdict(set)
        SIPHI_list = list()
        with open(self.predix + '_vs_SIPHI.xls', 'r') as SIPHI_r:
            _ = SIPHI_r.readline()
            for line in SIPHI_r.readlines():
                line = line.strip().split('\t')
                if len(line) >= 9:
                    SIPHI_s2hit[line[5]].add(line[9].split(';')[0])
                    SIPHI_list.append(line[5])

        sRNAMap_s2hit = defaultdict(set)
        sRNAMap_list = list()
        with open(self.predix + '_vs_sRNAMap.xls', 'r') as sRNAMap_r:
            _ = sRNAMap_r.readline()
            for line in sRNAMap_r.readlines():
                line = line.strip().split('\t')
                if len(line) >= 9:
                    sRNAMap_s2hit[line[5]].add(line[9].split(';')[0])
                    sRNAMap_list.append(line[5])

        sRNATarBase_s2hit = defaultdict(set)
        sRNATarBase_list = list()
        with open(self.predix + '_vs_sRNATarBase.xls', 'r') as sRNATarBase_r:
            _ = sRNATarBase_r.readline()
            for line in sRNATarBase_r.readlines():
                line = line.strip().split('\t')
                if len(line) >= 9:
                    sRNATarBase_s2hit[line[5]].add(line[9].split(';')[0])
                    sRNATarBase_list.append(line[5])

        # BSRD_s2hit = defaultdict(set)
        # BSRD_list = list()
        # with open(self.predix + '_vs_BSRD.xls', 'r') as BSRD_r:
        #     BSRD_r.readline()
        #     for line in BSRD_r:
        #         line = line.strip().split('\t')
        #         if len(line) >= 9:
        #             BSRD_s2hit[line[5]].add(line[9].split(';')[0])
        #             BSRD_list.append(line[5])

        all_srna = set(rfam_list + SIPHI_list + sRNAMap_list + sRNATarBase_list)    # +BSRD_list
        with open('annotation_stat.xls', 'w') as annot_stat:
            annot_stat.write('sRNA ID\tsRNATarBase\tsRNAMap\trfam\tSIPHI\tSum\n')
            stat_str = ''
            for s in sorted(all_srna):
                sum_ = 0
                if s in sRNATarBase_s2hit:
                    stat_str += s + '\t' + ';'.join(sRNATarBase_s2hit[s])
                    sum_ += 1
                else:
                    stat_str += s + '\t' + '-'
                if s in sRNAMap_s2hit:
                    stat_str += '\t' + ';'.join(sRNAMap_s2hit[s])
                    sum_ += 1
                else:
                    stat_str += '\t' + '-'
                if s in rfam_s2hit:
                    stat_str += '\t' + ';'.join(rfam_s2hit[s])
                    sum_ += 1
                else:
                    stat_str += '\t' + '-'
                if s in SIPHI_s2hit:
                    stat_str += '\t' + ';'.join(SIPHI_s2hit[s])
                    sum_ += 1
                else:
                    stat_str += '\t' + '-'
                # if s in BSRD_s2hit:
                #     stat_str += '\t' + ';'.join(BSRD_s2hit[s])
                #     sum_ += 1
                # else:
                #     stat_str += '\t' + '-'
                stat_str += '\t' + str(sum_) + '\n'
            annot_stat.write(stat_str)

    def run(self):
        self.run_blast()
        self.modify_rfam_description()
        self.xml2table(self.predix + '_vs_SIPHI.xml', self.predix + '_vs_SIPHI.xls')
        self.xml2table(self.predix + '_vs_sRNAMap.xml', self.predix + '_vs_sRNAMap.xls')
        self.xml2table(self.predix + '_vs_sRNATarBase.xml', self.predix + '_vs_sRNATarBase.xls')
        self.xml2table(self.predix + '_vs_BSRD.xml', self.predix + '_vs_BSRD.xls')
        self.get_list_file(self.predix + '_vs_rfam.xls', 'rfam')
        self.get_list_file(self.predix + '_vs_SIPHI.xls', 'SIPHI')
        self.get_list_file(self.predix + '_vs_sRNAMap.xls', 'sRNAMap')
        self.get_list_file(self.predix + '_vs_sRNATarBase.xls', 'sRNATarBase')
        self.get_list_file(self.predix + '_vs_BSRD.xls', 'BSRD')
        self.merge_stat_rfam()
        self.merge_annot(self.predix + '_vs_SIPHI.xls', self.db_path + 'all.sRNA.annot3.xls', 'SIPHI_merge')
        self.merge_annot(self.predix + '_vs_sRNATarBase.xls', self.db_path + 'sRNATarBase_V0002.annot.xls', 'sRNATarBase_merge')
        self.stat_annot()


if __name__ == '__main__':
    predict_fa = sys.argv[1]
    evalue = sys.argv[2]
    rfam = sys.argv[3]
    ANNO = srna_annot(predict_fa, evalue, rfam)
    ANNO.run()



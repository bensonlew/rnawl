# -*- coding: utf-8 -*-
# __author__ = 'sheng.he'

import os
import re
from Bio.Blast import NCBIXML
import sys
from Bio import SeqIO
from Bio import Seq


def get_xml_info(xmls):
    blast_hit_dict = {}
    matched_gene_set = set()
    for xml in xmls:
        records = NCBIXML.parse(open(xml))
        for rec in records:
            query = re.split(' ', rec.query, maxsplit=1)[0]
            if query in matched_gene_set:
                continue
            else:
                matched_gene_set.add(query)
            for align in rec.alignments:
                for hsp in align.hsps:
                    one_hsp = dict()
                    if align.hit_id.startswith("gnl|BL_ORD_ID|"):
                        # 修改rfam名称特殊情况
                        one_hsp['Hit-Name'] = align.hit_def.split(" ")[0]
                        one_hsp['Hit-Description'] = " ".join(align.hit_def.split(" ")[1:])
                    else:
                        one_hsp['Hit-Name'] = align.hit_id
                        one_hsp['Hit-Description'] = align.hit_def
                    one_hsp['Query-Name'] = query
                    one_hsp['Score'] = str(hsp.score)
                    one_hsp['Q-Frame'] = str(hsp.frame[0])
                    one_hsp['Q-Begin'] = str(hsp.query_start)
                    one_hsp['Q-End'] = str(hsp.query_end)
                    one_hsp['Hsp-Begin'] = str(hsp.sbjct_start)
                    one_hsp['Hsp-End'] = str(hsp.sbjct_end)
                    one_hsp['Hsp-Frame'] = str(hsp.frame[1])
                    if query in blast_hit_dict:
                        if hsp.score > blast_hit_dict[query]['Score']:
                            blast_hit_dict[query] = one_hsp
                        else:
                            pass
                    else:
                        blast_hit_dict[query] = one_hsp
    return blast_hit_dict

def extend_orf(pep, i, start, end):
    # 根据蛋白长度延伸比对区域
    pep_start = (start - 1 - i)/3
    pep_end = (end - 1 - i)/3
    pep_start_i  = i + 1
    pep_end_i = i + 1 + len(pep) * 3
    if pep_start >= 0:
        for n in range(pep_start, 0, -1):
            if pep[n] == "M":
                pep_start_i = i + 1 + n*3
                break
            elif pep[n] == "*":
                pep_start_i = i + 1 + (n + 1)*3
                break
            else:
                pass
    if pep_end <= len(pep):
        for n in range(pep_end - 1, len(pep)):
            if pep[n] == "*":
                pep_end_i = i + 1 + (n+1)*3
                break
            else:
                pass
    return pep_start_i, pep_end_i



def get_orf_by_blast(dna_seq, blast_dict, prefix="blast_hit"):
    record = SeqIO.parse(open(dna_seq), "fasta")
    with open(prefix + '.xml_orf.xls', 'w') as f, \
         open(prefix + '.cds.fa', 'w') as f1, \
         open(prefix + '.unhit.fa', 'w') as f3, \
         open(prefix + '.pep.fa', 'w') as f2:

        for seq in record:
            if seq.id in blast_dict:
                hsp = blast_dict[seq.id]
                if hsp['Q-Frame'] in ["3", "1", "2"]:
                    i = (int(hsp['Q-Frame']) - 1) % 3
                    frames_pep = Seq.translate(seq.seq[i:], table='Standard', stop_symbol='*')
                    pep_start, pep_end = extend_orf(frames_pep, i, int(hsp['Q-Begin']), int(hsp['Q-End']))
                    # f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(seq.id, hsp['Hit-Name'], pep_start, pep_end, hsp['Q-Frame'], frames_pep, Seq.translate(str(seq.seq[pep_start-1:pep_end-1]), table='Standard', stop_symbol='*')))
                    pep = Seq.translate(str(seq.seq[pep_start-1:pep_end-1]), table='Standard', stop_symbol='*')
                    if pep.startswith("M") and pep.endswith("*"):
                        orf_type = "complete"
                    elif pep.startswith("M"):
                        orf_type = "5prime_partial"
                    elif pep.endswith("*"):
                        orf_type = "3prime_partial"
                    else:
                        orf_type = "internal"
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(seq.id + "_orf1", hsp['Hit-Name'], pep_start, pep_end - 1, hsp['Q-Frame'], orf_type, "+"))
                    des = "ORF type:{} len:{} hit:{} {}:{}-{}({})".format(orf_type, len(pep), hsp['Hit-Name'], seq.id, pep_start, pep_end - 1, "+")
                    f1.write(">{} {}\n{}\n".format(seq.id + "_orf1", des, seq.seq[pep_start-1:pep_end-1]))
                    f2.write(">{} {}\n{}\n".format(seq.id + "_orf1", des, pep))
                elif hsp['Q-Frame'] in ["-1", "-2", "-3"]:
                    seq_rc = str(seq.seq.reverse_complement())
                    i =  (2 - int(hsp['Q-Frame'])) % 3
                    frames_pep = Seq.translate(seq_rc[i:], table='Standard', stop_symbol='*')
                    pep_start, pep_end = extend_orf(frames_pep, i, len(seq) + 2 - int(hsp['Q-End']), len(seq) + 2 - int(hsp['Q-Begin']))


                    pep = Seq.translate(seq_rc[pep_start-1:pep_end-1], table='Standard', stop_symbol='*')
                    if pep.startswith("M") and pep.endswith("*"):
                        orf_type = "complete"
                    elif pep.startswith("M"):
                        orf_type = "5prime_partial"
                    elif pep.endswith("*"):
                        orf_type = "3prime_partial"
                    else:
                        orf_type = "internal"
                    des = "ORF type:{} len:{} hit:{} {}:{}-{}({})".format(orf_type, len(pep), hsp['Hit-Name'], seq.id, len(seq) + 2 - pep_end, len(seq) + 2 - pep_start - 1, "-")

                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(seq.id + "_orf1", hsp['Hit-Name'],len(seq) + 2 -  pep_end, len(seq) + 2 - pep_start - 1, hsp['Q-Frame'], orf_type, "-"))
                    f1.write(">{} {}\n{}\n".format(seq.id + "_orf1", des,  seq_rc[pep_start-1:pep_end-1]))
                    f2.write(">{} {}\n{}\n".format(seq.id + "_orf1", des,  pep))
                else:
                    print "query frame wrong {}".format(hsp['Q-Frame'])
            else:
                f3.write(">{}\n{}\n".format(seq.id, str(seq.seq)))



if __name__ == '__main__':  # for test
    blast_files = sys.argv[1].split(",")
    transcript = sys.argv[2]
    blast_dict = get_xml_info(blast_files)
    get_orf_by_blast(transcript, blast_dict)

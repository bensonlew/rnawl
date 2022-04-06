# !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "yitong.feng"
# 20180721


import os
import sys
# from mako.template import Template
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def get_gene_range(predict_rna):
    """
    获取基因cds在坐标信息
    """
    gene_range = dict()
    with open(predict_rna, 'r') as trans_r:
        for line in trans_r.readlines()[1:]:
            cols = line.strip("\n").split("\t")
            if len(cols) >= 3 and cols[1] and cols[2]:
                gene_range[cols[6]] = (min(int(cols[1]), int(cols[2])),
                                       max(int(cols[1]), int(cols[2])),
                                       cols[5],
                                       cols[7])
    return gene_range


def choose_anti_info(anti_info, anti):
    if len(anti_info) == 1:
        choosed_aanti_info = anti_info[0]

    else:
        choosed_aanti_info = None
        for aanti_info in anti_info:
            if aanti_info[1] > anti["end"] or aanti_info[2] < anti["start"]:
                continue
            else:
                choosed_aanti_info = aanti_info
    return choosed_aanti_info

def get_seq_anti(antisense, unanti, antisense_f):
    '''
    写入反义转录本信息
    '''
    gene1_info = []
    gene2_info = []
    for anti in antisense:
        if anti["strand"] == "+":
            gene1_info = [anti["id"], anti["start"],
                          anti["end"], anti["description"]]
        else:
            gene2_info = [anti["id"], anti["start"],
                          anti["end"], anti["description"]]

        if anti["antisense_id"] in unanti:
            anti_info = unanti[anti["antisense_id"]]
            choosed_aanti_info = choose_anti_info(anti_info, anti)           

            if choosed_aanti_info:
                if choosed_aanti_info[-2] == anti["strand"]:
                    print("error : strand is same")
                elif anti["strand"] == "+":
                    gene2_info = [choosed_aanti_info[0], choosed_aanti_info[1],
                                      choosed_aanti_info[2], choosed_aanti_info[-1]]
                else:
                    gene1_info = [choosed_aanti_info[0], choosed_aanti_info[1],
                                      choosed_aanti_info[2], choosed_aanti_info[-1]]
            # print gene1_info, gene2_info

            over_laps = [int(gene1_info[1]), int(gene1_info[2]),
                         int(gene2_info[1]), int(gene2_info[2])]
            over_laps.sort()
            if int(gene1_info[1]) <= int(gene2_info[1]) and int(gene1_info[2]) >= int(gene2_info[2]):
                anti_type = "enclosed"
            elif int(gene1_info[1]) >= int(gene2_info[1]) and int(gene1_info[2]) <= int(gene2_info[2]):
                anti_type = "enclosed"
            elif int(gene1_info[1]) < int(gene2_info[1]) and int(gene1_info[2]) < int(gene2_info[2]):
                anti_type = "enclosed"
            else:
                anti_type = "divergent"

            antisense_f.write("\t".join(gene1_info + gene2_info + [
                str(over_laps[1]),
                str(over_laps[2]),
                str(over_laps[2] - over_laps[1] + 1),
                anti_type
            ]
            ) + "\n")
    return True

def parse_operons(seq, gene_range, file_dict):
    # opera 信息提取
    with open('Rockhopper_Results/' + seq.split('.')[0] + '_operons.txt', 'r') as ope_r:
        header_ope = ope_r.readline()
        for line in ope_r.readlines():
            cols = line.strip("\n").split("\t")
            gene_detail = []
            genes = cols[-1].split(", ")
            for gene_id in gene_range:
                if gene_range[gene_id][0] <= int(cols[1]) and gene_range[gene_id][1] >= int(cols[0]):
                    if gene_id in genes or gene_range[gene_id][2] in genes:
                        gene_detail.append(
                            gene_id + "(" + gene_range[gene_id][2] + "|" + gene_range[gene_id][3] + ")")
            file_dict["ope_w"].write(seq + '\t' + line.strip() + "\t" +";".join(gene_detail) + "\n")


def parse_new_genes(file_dict, new_cds_num, seq, tran_start, tran_end, type__, strand, length):
    # print("new_cds_num {}".format(new_cds_num))
    type_new = "novel{}".format(str(new_cds_num).zfill(4))
    file_dict["genome_cds_prebed"].write(
        seq + '\t' + str(int(tran_start) - 1) + '\t' + tran_end + '\t' + type_new + '\t0' + '\t' + strand + '\t' + str(int(tran_start) - 1) + '\t' + tran_end + '\t0' + '\t1' + '\t' + str(length) + '\t0' + '\n')
    des_new = blast_result[type__][1]
    predicted_id = type_new
    if strand == '-':
        file_dict["genome_cds_prexls"].write(
            seq + '\t' + tran_end + '\t' + str(int(tran_start)) + '\t' + type_new + '\t' + 'predicted_cds' + '\t' +
            strand + '\t' + des_new + '\t' + str(length) + '\n')
    else:
        file_dict["genome_cds_prexls"].write(
            seq + '\t' + str(int(tran_start)) + '\t' + tran_end + '\t' + type_new + '\t' + 'predicted_cds' + '\t' +
            strand + '\t' + des_new + '\t' + str(length) + '\n')
    file_dict["fasta_cds_prebed"].write(type__ + '\t0' + '\t' + str(length) + '\t' + type__ + '\t0' + '\t' +
                        strand + '\t0' + '\t' + str(length) + '\t0' + '\t1' + '\t' + str(length) + '\t0' + '\n')
    
    file_dict["genome_cds_xls"].write(type_new + '\t' + seq + '\t' + tran_start + '\t' + tran_end + '\t' + strand + '\t'  + str(length) + '\n')
    return predicted_id

def parse_srna(file_dict, type_, seq, tran_start, tran_end, type__, strand, length):
    des = type_
    file_dict["genome_prebed"].write(
        seq + '\t' + str(int(tran_start) - 1) + '\t' + tran_end + '\t' + type__ + '\t0' + '\t' + strand + '\t' + str(int(tran_start) - 1) + '\t' + tran_end + '\t0' + '\t1' + '\t' + str(length) + '\t0' + '\n')
    if strand == '-':
        file_dict["genome_prexls"].write(
            seq + '\t' + tran_end + '\t' + str(int(tran_start)) + '\t' + type__ + '\t' + 'predicted_RNA' + '\t' +
            strand + '\t' + des + '\t' + str(length) + '\n')
    else:
        file_dict["genome_prexls"].write(
            seq + '\t' + str(int(tran_start)) + '\t' + tran_end + '\t' + type__ + '\t' + 'predicted_RNA' + '\t' +
            strand + '\t' + des + '\t' + str(length) + '\n')
    file_dict["fasta_prebed"].write(type__ + '\t0' + '\t' + str(length) + '\t' + type__ + '\t0' + '\t' +
                    strand + '\t0' + '\t' + str(length) + '\t0' + '\t1' + '\t' + str(length) + '\t0' + '\n')

def parse_cds(file_dict, type_, seq, tran_start, tran_end, code_start, code_end, strand, name, des):
    tran_start_org = tran_start
    tran_end_org = tran_end
    if tran_start == '':
        tran_start = code_start
    if tran_end == '':
        tran_end = code_end
    tran_start_tmp, tran_end_tmp, code_start_tmp, code_end_tmp = tran_start, tran_end, code_start, code_end
    if strand == '-':
        tran_end, tran_start = tran_start, tran_end
        code_end, code_start = code_start, code_end
    length = int(tran_end) - int(tran_start) + 1
    file_dict["genome_bed"].write(seq + '\t' + str(int(tran_start)-1) + '\t' + tran_end + '\t' + type_ + '\t0' + '\t' + strand +
                    '\t' + str(int(code_start)-1) + '\t' + code_end + '\t0' + '\t1' + '\t' + str(length) + '\t0' + '\n')
    file_dict["fasta_bed"].write(type_ + '\t' + '0\t' + str(length) + '\t' + type_ + '\t0' + '\t' + strand + '\t' + str(int(code_start)-int(
        tran_start)) + '\t' + str(int(code_end)-int(tran_start)+1) + '\t0' + '\t1' + '\t' + str(length) + '\t0' + '\n')
    file_dict["tss_tts"].write(type_ + '\t' + name + '\t' + des + '\t' + seq + '\t' + strand + '\t' +
                tran_start_org + '\t' + code_start_tmp + '\t' + code_end_tmp + '\t' + tran_end_org + '\n')
    if line[0] != '' and line[0] != line[1]:
        file_dict["utr"].write(seq + '\t' + tran_start_tmp + '\t' + code_start_tmp + '\t' +
                type_ + '\t' + name + '\t' + strand + '\t' + 'UTR5\t' + des + '\n')
        if strand == '+':
            file_dict["utr5"].write(seq + '\t' + str(int(tran_start_tmp)-1) + '\t' +
                    code_start_tmp + '\t' + type_ + '_UTR5' + '\t0' + '\t' + strand + '\n')
        else:
            file_dict["utr5"].write(seq + '\t' + str(int(
                code_start_tmp) - 1) + '\t' + tran_start_tmp + '\t' + type_ + '_UTR5' + '\t0' + '\t' + strand + '\n')
    if line[3] != '' and line[2] != line[3]:
        file_dict["utr"].write(
            seq + '\t' + code_end_tmp + '\t' + tran_end_tmp + '\t' + type_ + '\t' + name + '\t' + strand + '\t' + 'UTR3\t' + des + '\n')
        if strand == '+':
            file_dict["utr3"].write(seq + '\t' + str(int(code_end_tmp)-1) + '\t' +
                    tran_end_tmp + '\t' + type_ + '_UTR3' + '\t0' + '\t' + strand + '\n')
        else:
            file_dict["utr3"].write(seq + '\t' + str(int(
                tran_end_tmp) - 1) + '\t' + code_end_tmp + '\t' + type_ + '_UTR3' + '\t0' + '\t' + strand + '\n')



def parse_rna(file_dict, type_, seq, tran_start, tran_end, strand):
    if tran_end != '' and tran_start != '':
        if strand == '-':
            tran_end, tran_start = tran_start, tran_end
        # print(tran_end, tran_start)
        length = int(tran_end) - int(tran_start) + 1
        file_dict["genome_kbed"].write(seq + '\t' + str(int(tran_start)-1) + '\t' + tran_end + '\t' + type_ + '\t0' + '\t' + strand + '\t' + str(
            int(tran_start)-1) + '\t' + tran_end + '\t0' + '\t1' + '\t' + str(length) + '\t0' + '\n')
        file_dict["fasta_kbed"].write(type_ + '\t' + '0\t' + str(length) + '\t' + type_ + '\t0' + '\t' +
                        strand + '\t0' + '\t' + str(length) + '\t0' + '\t1' + '\t' + str(length) + '\t0' + '\n')


def parse_antisense(type_, tran_start, tran_end, code_start, code_end, strand, des, des_new, predicted_id, antisense, name, unanti):
    if tran_start == "":
        tran_start = code_start
    if tran_end == "":
        tran_end = code_end
    if 'antisense: ' in des:
        antisense.append({
            "id": predicted_id,
            "start": tran_start,
            "end": tran_end,
            "strand": strand,
            "antisense_id": des.split("antisense: ")[1],
            "description": des_new

        })
    else:
        if name in unanti:
            unanti[name].append(
                [type_, tran_start, tran_end, strand, des]
            )
        else:
            unanti[name] = [
                [type_, tran_start, tran_end, strand, des]]
 

def parse_transcript(seq, gene_range, file_dict, new_num, new_cds_num, antisense, unanti):
    with open('Rockhopper_Results/' + seq.split('.')[0] + '_transcripts.txt', 'r') as trans_r:
        trans_r.readline()
        for line in trans_r.readlines():
            line = line.strip('\n').split('\t')
            if not len(line) < 6:
                tran_start = line[0]
                code_start = line[1]
                code_end = line[2]
                tran_end = line[3]
                strand = line[4]
                name = line[5]
                type_ = line[6]
                des = line[7]
                des_new = des
                predicted_id = ""

                if u'predicted' in type_:
                    if strand == '-':
                        tran_end, tran_start = tran_start, tran_end
                    length = int(tran_end) - int(tran_start) + 1
                    # id需要固定顺序以确保和blast结果一致
                    type__ = 'sRNA' + str(new_num).zfill(4)
                    predicted_id = type__
                    new_num += 1
                    if type__ in blast_result:
                        # "new gene
                        new_cds_num += 1
                        predicted_id = parse_new_genes(file_dict, new_cds_num, seq, tran_start, tran_end, type__, strand, length)
                    elif length <= 500 and length >= 30:
                        parse_srna(file_dict, type_, seq, tran_start, tran_end, type__, strand, length)

                elif code_start != '' and code_end != '':
                    parse_cds(file_dict, type_, seq, tran_start, tran_end, code_start, code_end, strand, name, des)
                else:
                    parse_rna(file_dict, type_, seq, tran_start, tran_end, strand)
                parse_antisense(type_, tran_start, tran_end, code_start, code_end, strand, des, des_new, predicted_id, antisense, name, unanti)


def parse_one_seq(seq, gene_range, file_dict, new_num, new_cds_num):
    antisense = list()
    unanti = dict()
    parse_operons(seq, gene_range, file_dict)
    parse_transcript(seq, gene_range, file_dict, new_num, new_cds_num, antisense, unanti)
    get_seq_anti(antisense, unanti, file_dict["antisense_f"])



def extract_rockhopper_info():
    with open('Rockhopper_Results/genome.gene.bed', 'w') as genome_bed, \
        open('Rockhopper_Results/genome.knownnc.bed', 'w') as genome_kbed, \
        open('Rockhopper_Results/fasta.gene.bed', 'w') as fasta_bed, \
        open('Rockhopper_Results/fasta.knownnc.bed', 'w') as fasta_kbed, \
        open('Rockhopper_Results/genome.predicted_RNA.bed', 'w') as genome_prebed, \
        open('Rockhopper_Results/genome.predicted_RNA.bed.tmp', 'w') as genome_prexls, \
        open('Rockhopper_Results/fasta.predicted_RNA.bed', 'w') as fasta_prebed, \
        open('Rockhopper_Results/genome.predicted_cds.bed', 'w') as genome_cds_prebed, \
        open('Rockhopper_Results/genome.predicted_cds.bed.xls', 'w') as genome_cds_prexls, \
        open('Rockhopper_Results/genome.predicted_cds.xls', 'w') as genome_cds_xls, \
        open('Rockhopper_Results/fasta.predicted_cds.bed', 'w') as fasta_cds_prebed, \
        open('Rockhopper_Results/antisense.xls', 'w') as antisense_f, \
        open('Rockhopper_Results/operon.xls', 'w') as ope_w, \
        open('Rockhopper_Results/TSS_and_TTS.xls', 'w') as tss_tts, \
        open('Rockhopper_Results/UTR.xls', 'w') as utr, \
        open('Rockhopper_Results/UTR5.bed', 'w') as utr5, \
        open('Rockhopper_Results/UTR3.bed', 'w') as utr3:

        #  add file head
        ope_w.write(
            "Chr_name\tStart\tStop\tStrand\tNumber_of_genes\tGenes\tGenes_detail\n")
        tss_tts.write('Id\tName\tDescription\tLocation\tStrand\tTranscription start site\tTranslation initiation site\tTranslation stop site\tTranscription terminator site\n')
        utr.write("Chr\tStart\tEnd\tGeneID\tGeneName\tStrand\tType\tGene Desc\n")
        genome_prexls.write(
            'Location\tStart\tEnd\tsRNA ID\tType\tStrand\tAntisense_Genes\tlength\n')
        genome_cds_prexls.write(
            'Location\tStart\tEnd\tpredict ID\tType\tStrand\tDescription\tlength\n')
        genome_cds_xls.write("Gene ID\tLocation\tStart\tEnd\tStrand\tLength\n")
        antisense_f.write("\t".join([
            "GeneID(+)",
            "Start(+)",
            "End(+)",
            "Description(+)",
            "GeneID(-)",
            "Start(-)",
            "End(-)",
            "Description(-)",
            "Overlap_start",
            "Overlap_end",
            "Overlap_length",
            "Type"
        ]) + "\n")


        file_dict = {
            "genome_bed": genome_bed,
            "genome_kbed": genome_kbed,
            "fasta_bed": fasta_bed,
            "fasta_kbed": fasta_kbed,
            "genome_prebed": genome_prebed,
            "genome_prexls": genome_prexls,
            "fasta_prebed": fasta_prebed,
            "genome_cds_prebed": genome_cds_prebed,
            "genome_cds_prexls": genome_cds_prexls,
            "genome_cds_xls": genome_cds_xls,
            "fasta_cds_prebed": fasta_cds_prebed,
            "antisense_f": antisense_f,
            "ope_w": ope_w,
            "tss_tts":tss_tts,
            "utr" : utr,
            "utr5": utr5,
            "utr3": utr3
        }

        new_num = 0
        new_cds_num = 0
        for seq in seqs:
            gene_range = get_gene_range('Rockhopper_Results/' + seq.split('.')[0] + '_transcripts.txt')

            parse_one_seq(seq, gene_range, file_dict, new_num, new_cds_num)

            


def extract_fasta():
    #os.system('cat Rockhopper_Results/genome.gene.bed Rockhopper_Results/genome.knownnc.bed Rockhopper_Results/genome.predicted_RNA.bed >Rockhopper_Results/genome.feature.bed && cat Rockhopper_Results/fasta.gene.bed Rockhopper_Results/fasta.knownnc.bed Rockhopper_Results/fasta.predicted_RNA.bed >Rockhopper_Results/fasta.feature.bed')
    # os.system('cat Rockhopper_Results/genome.gene.bed Rockhopper_Results/genome.predicted_RNA.bed >Rockhopper_Results/genome.feature.bed && cat Rockhopper_Results/fasta.gene.bed Rockhopper_Results/fasta.knownnc.bed Rockhopper_Results/fasta.predicted_RNA.bed >Rockhopper_Results/fasta.feature.bed')
    os.system('cat Rockhopper_Results/genome.gene.bed Rockhopper_Results/genome.predicted_RNA.bed >Rockhopper_Results/genome.feature.bed && cat Rockhopper_Results/fasta.gene.bed Rockhopper_Results/fasta.predicted_RNA.bed >Rockhopper_Results/fasta.feature.bed')
    bedtool_path = script_path + '/bioinfo/seq/bedtools-2.25.0/bin/bedtools'
    cmd = '''${bedtool_path} getfasta -fi ${fna} -bed Rockhopper_Results/genome.feature.bed -s -name -fo Rockhopper_Results/genome.feature.fa
    ${bedtool_path} getfasta -fi ${fna} -bed Rockhopper_Results/genome.gene.bed -s -name -fo Rockhopper_Results/genome.gene.fa
    ${bedtool_path} getfasta -fi ${fna} -bed Rockhopper_Results/genome.predicted_RNA.bed -s -name -fo Rockhopper_Results/genome.predicted_RNA.fa
    ${bedtool_path} getfasta -fi ${fna} -bed Rockhopper_Results/UTR5.bed -s -name -fo Rockhopper_Results/UTR5.fa
    ${bedtool_path} getfasta -fi ${fna} -bed Rockhopper_Results/UTR3.bed -s -name -fo Rockhopper_Results/UTR3.fa
    '''

    # f = Template(cmd)
    # bash_info = f.render(bedtool_path=bedtool_path,
    #                      fna=fna,
    #                      )
    with open('run_rockhopper2bed.bash', 'w') as rock_sh:
        rock_sh.write("{} getfasta -fi {} -bed Rockhopper_Results/genome.feature.bed -s -name -fo Rockhopper_Results/genome.feature.fa".format(bedtool_path, fna)+'\n' +
                    "{} getfasta -fi {} -bed Rockhopper_Results/genome.gene.bed -s -name -fo Rockhopper_Results/genome.gene.fa".format(bedtool_path, fna) + '\n' +
                    "{} getfasta -fi {} -bed Rockhopper_Results/genome.predicted_RNA.bed -s -name -fo Rockhopper_Results/genome.predicted_RNA.fa".format(bedtool_path, fna) + '\n' +
                    "{} getfasta -fi {} -bed Rockhopper_Results/genome.predicted_cds.bed -s -name -fo Rockhopper_Results/genome.predicted_cds.fa".format(bedtool_path, fna) + '\n' +
                    "{} getfasta -fi {} -bed Rockhopper_Results/UTR5.bed -s -name -fo Rockhopper_Results/UTR5.fa".format(bedtool_path, fna) + '\n' +
                    "{} getfasta -fi {} -bed Rockhopper_Results/UTR3.bed -s -name -fo Rockhopper_Results/UTR3.fa".format(
                        bedtool_path, fna)
                    )
        # rock_sh.write(bash_info)

    os.system('bash run_rockhopper2bed.bash')

    #  ------给genome.predicted_RNA.bed加入seqence那一列-----

    with open('Rockhopper_Results/genome.predicted_RNA.bed.tmp', 'r') as pre_tmp, \
            open('Rockhopper_Results/genome.predicted_RNA.fa', 'r') as pre_fa, \
            open('Rockhopper_Results/genome.predicted_RNA.bed.xls', 'w') as pre_xls:
        s2seq = dict()
        for seq in pre_fa.read().split('\n>'):
            seq = seq.strip().split('\n')
            s2seq[seq[0].strip().lstrip('>')] = ''.join(seq[1:])
        pre_xls.write(pre_tmp.readline().strip() + '\tSequence\n')
        for line in pre_tmp.readlines():
            tmp = line.strip().split('\t')
            try:
                pre_xls.write(line.strip() + '\t' + s2seq[tmp[3]] + '\n')
            except:
                pre_xls.write(line.strip() + '\t' + '-' + '\n')

    # ---------通过ptt文件生成基因的fasta文件并将其翻译成碱基序列--------

        chr_list = os.listdir('rock_index')
        for i in chr_list:
            with open('rock_index/' + i + '/' + i + '.ptt', 'r') as ptt_r, \
                    open('Rockhopper_Results/' + 'cds.bed', 'a') as bed_w, \
                    open('Rockhopper_Results/' + 'ptt.bed', 'a') as ptt_w:
                _ = ptt_r.readline()
                _ = ptt_r.readline()
                _ = ptt_r.readline()
                for line in ptt_r.readlines():
                    line = line.strip('\n').split('\t')
                    bed_w.write(i + '\t' + str(int(line[0].split('..')[0])-1) + '\t' + line[0].split(
                        '..')[1] + '\t' + line[5] + '\t0' + '\t' + line[1] + '\n')
                    ptt_w.write(i + '\t' + '\t'.join(line) + '\n')
        cmd = "{bedtool_path} getfasta -fi {fna} -bed Rockhopper_Results/cds.bed -s -name -fo Rockhopper_Results/cds.fa".format(
            bedtool_path=bedtool_path, fna=fna)
        os.system(cmd)

    with open('Rockhopper_Results/cds.fa', 'r') as fa_r, \
            open('Rockhopper_Results/cds.faa', 'w') as faa_w:
        for block in fa_r.read().split('\n>'):
            block = block.lstrip('>').split('\n')
            coding_dna = Seq(''.join(block[1:]), IUPAC.ambiguous_dna)
            protein = coding_dna.translate()
            faa_w.write('>' + block[0].strip() + '\n' + str(protein) + '\n')

    # 不添加蛋白序列可能不完整

    with open('Rockhopper_Results/genome.predicted_cds.bed.tmp', 'r') as ptt_r, \
        open('Rockhopper_Results/' + 'cds.bed', 'aw') as bed_w, \
        open('Rockhopper_Results/' + 'ptt.bed', 'aw') as ptt_w:
        ptt_r.readline()
        for line in ptt_r.readlines():
            line = line.strip('\n').split('\t')
            bed_w.write(line[0] + '\t' + line[1] + '\t' + line[2] +
                        '\t' + line[3] + '\t0' + '\t' + line[4] + '\n')
            ptt_w.write('\t'.join(line) + '\n')


    with open('Rockhopper_Results/genome.predicted_cds.fa', 'r') as fa_r, \
            open('Rockhopper_Results/cds.fa', 'aw') as fa_w:
        for line in fa_r:
            fa_w.write(line)



if __name__ == '__main__':


    fna = sys.argv[1]
    script_path = sys.argv[2]
    blast_path = sys.argv[3]
    operons = list()
    transcripts = list()
    rock_results = os.listdir('Rockhopper_Results')
    seqs_op = set()
    seqs_tr = set()

    for file in rock_results:
        if u'transcripts' in file:
            transcripts.append(file)
        if u'operons' in file:
            operons.append(file)

    for file in operons:
        seqs_op.add(file.split('_operons.txt')[0])
    for file in transcripts:
        seqs_tr.add(file.split('_transcripts.txt')[0])
    seqs_tmp = seqs_op & seqs_tr
    seqs = set()
    for i in os.listdir('rock_index'):
        if '.' in i:
            for j in seqs_tmp:
                if j in i:
                    seqs.add(i)
        else:
            for j in seqs_tmp:
                if j == i:
                    seqs.add(i)

    blast_result = dict()
    if os.path.exists(blast_path):
        with open(blast_path, 'r') as f:
            for line in f:
                seq_id = line.split('\t')[0]
                target_id = line.split('\t')[2]
                des = " ".join(line.split("\t")[6].split()[1:])
                blast_result[seq_id] = [target_id, des]

    extract_rockhopper_info()
    extract_fasta()
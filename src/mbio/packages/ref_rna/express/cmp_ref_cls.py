# -*- coding: utf-8 -*-
# __author__ = fiona
# modify by khl
# time: 2017/3/30 14:50

import re, os, Bio, argparse, sys, fileinput, urllib2, subprocess


def get_dic_from_merged_gtf(merged_file):
    d = {}
    for line in open(merged_file):
            # if re.search(r'\t+transcript\t+', line):
            cls_m = re.search(r'class_code\s+\"(\S+)\";', line)
            # ref_txpt_id_m = re.search(r'cmp_ref\s+\"(\S+)\";', line)
            assembly_txpt_id_m = re.search(r'transcript_id\s+\"(\S+)\";', line)
            assembly_gene_id_m = re.search(r'gene_id\s+\"(\S+)\";', line)
            assembly_txpt_id = ''
            assembly_gene_id = ''
            # ref_txpt_id = ''
            cls = ''
            line_id = line.strip().split("\t")
            if assembly_txpt_id_m and assembly_gene_id_m:
                assembly_txpt_id = assembly_txpt_id_m.group(1)
                assembly_gene_id = assembly_gene_id_m.group(1)
                if cls_m:
                    cls = cls_m.group(1)
                else:
                    cls = '='
                # if re.match(r'^[^u]$', cls) and not ref_txpt_id_m:
                #     raise Exception(
                #         '{} of combined gtf file has logical problem: class code 为非U得情况没有cmp_ref的值 '.format(line))
                # if re.match(r'^[u]$', cls) and ref_txpt_id_m:
                #     raise Exception(
                #         '{} of combined gtf file has logical problem: class code 为U得情况有nearest_ref的值 '.format(line))
                # if ref_txpt_id_m:
                #     ref_txpt_id = ref_txpt_id_m.group(1)
            # else:
            #     raise Exception(
            #         'line: {} in annotate gtf {} has no valid internal txpt id or class code'.format(line.strip(),merged_file))
                d[assembly_txpt_id] = {'assembly_txpt_id': assembly_txpt_id, 'cls': cls,
                                       'assembly_gene_id': assembly_gene_id}
    return d


def get_gname_gid_dic_from_ref_gtf(ref_gtf_content):
    gene_id2gene_name = dict()
    for line in ref_gtf_content:
        gene_id_m = re.search(r'gene_id\s+\"(\S+)\"', line.strip())
        gname_m = re.search(r'gene_name\s+\"(\S+)\"', line.strip())
        txpt_id_m = re.search(r'transcript_id\s+\"(\S+)\"', line.strip())
        if gene_id_m and txpt_id_m:
            gene_id = gene_id_m.group(1)
            # txpt_id = txpt_id_m.group(1)
            if gname_m:
                gname = gname_m.group(1)
                if len(gname) > 1:
                    gene_id2gene_name[gene_id] = gname
            # d[txpt_id] = {'gene_id': gene_id, 'gname': gname}
        else:
            pass
            # raise Exception('ref gtf文件不合法的第九列: {},没有完整的gene id 和gene_name 记录'.format(line.strip()))
    return gene_id2gene_name


def write_text(ref_dic, merged_dic, f):
    fw = open(f, 'w')
    head = '#{}\n'.format('\t'.join(
        ['assemble_txpt_id', 'assemble_gene_id', 'class_code', 'ref_gene_name']))
    fw.write(head)

    for internal_id in merged_dic.keys():
        cls = merged_dic[internal_id]['cls']
        assembly_gene_id = merged_dic[internal_id]['assembly_gene_id']
        # assembly_txpt_id = merged_dic[internal_id]['assembly_txpt_id']
        # chr = merged_dic[internal_id]["chr"]
        # strand = merged_dic[internal_id]["strand"]
        # start = merged_dic[internal_id]["start"]
        # end = merged_dic[internal_id]["end"]

        ref_gname = '-'
        try:
            # if internal_id in ref_dic.keys():  # add 1 line by khl 20170412
            #     # ref_gid = ref_dic[assembly_txpt_id]['gene_id']
            #     ref_gname = ref_dic[internal_id]['gname']
            # else:
            #     ref_gname = '-'
            if assembly_gene_id in ref_dic:
                ref_gname = ref_dic[assembly_gene_id]
        except Exception:
            # print assembly_txpt_id
            pass
        newline = '\t'.join([internal_id, assembly_gene_id, cls, ref_gname]) + '\n'
        fw.write(newline)

    fw.close()


def write_relation_text_for_merge_gtf(ref_gtf, merged_gtf, text_file):
    ref_tmp_cmd = """ awk -F '\t' 'NF>=9{print $9}' %s  | uniq  """ % (ref_gtf)
    ref_tmp_content = subprocess.check_output(ref_tmp_cmd, shell=True).strip().split('\n')
    ref_t_gid_gname_dic = get_gname_gid_dic_from_ref_gtf(ref_tmp_content)
    merged_tid_ref_tid_cls_dic = get_dic_from_merged_gtf(merged_gtf)
    write_text(ref_t_gid_gname_dic, merged_tid_ref_tid_cls_dic, text_file)


if __name__ == '__main__':
    # ref = '/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/ref.gtf'
    # cmp_f = '/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/cuffmerge_gffcompare/gffcmp_cuffmerge_.annotated.gtf'

    # out = '/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/cuffmerge_gffcompare/relation_table.tsv'
    # ref = "/mnt/ilustre/users/sanger-dev/sg-users/qindanhua/test_files/ref_rna/pipe/ref/ref_genome.gtf"
    # cmp_f = "/mnt/ilustre/users/sanger-dev/workspace/20170508/Single_assembly_module_tophat_stringtie_true_file/RefrnaAssemble/output/Gffcompare/cuffcmp.annotated.gtf"
    # out = "/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tofiles/newdata_class_code"
    ref = "/mnt/ilustre/users/sanger-dev/sg-users/qindanhua/test_files/ref_rna/pipe/ref/ref_genome.gtf"

    # cmp_f = "/mnt/ilustre/users/sanger-dev/workspace/20170508/Single_assembly_module_tophat_stringtie_true_file/RefrnaAssemble/output/Gffcompare/cuffcmp.annotated.gtf"
    # merged_f = "/mnt/ilustre/users/sanger-dev/workspace/20170616/Single_assembly_module_tophat_stringtie_true_file_1/RefrnaAssemble/NewTranscripts/output/change_id_merged.gtf"
    merged_f = '/mnt/ilustre/users/sanger-dev/workspace/20170508/Single_assembly_module_tophat_stringtie_true_file/RefrnaAssemble/output/StringtieMerge/merged.gtf'
    out = "/mnt/ilustre/users/sanger-dev/workspace/20170523/Single_merge_rsem_fpkm_11/MergeRsem/new_class_code"

    write_relation_text_for_merge_gtf(ref, merged_f, out)
    # ref = "/mnt/ilustre/users/sanger-dev/workspace/20170210/Refrna_refrna_test_01/FilecheckRef/Danio_rerio.GRCz10.85.gff3.gtf"
    # cmp_f = "/mnt/ilustre/users/sanger-dev/workspace/20170411/Single_assembly_module_tophat_cufflinks/Assembly/Gffcompare/output/cuffcmp.annotated.gtf"
    # out = "/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/module/cufflinks_class_code/class_code"
    # write_relation_text_for_merge_gtf("cufflinks_cmp_1000.gtf", "ref_1000.gtf", "class_code")

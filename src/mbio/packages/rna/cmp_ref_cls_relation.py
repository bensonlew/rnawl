# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/30 14:50

import re, os, Bio, argparse, sys, fileinput, urllib2, subprocess


def get_dic_from_cmp_gtf(cmp_file):
    d = {}
    for line in open(cmp_file):
        if re.search(r'\t+transcript\t+', line):
            cls_m = re.search(r'class_code\s+\"(\S+)\";', line)
            ref_txpt_id_m = re.search(r'cmp_ref\s+\"(\S+)\";', line)
            internal_txpt_id_m = re.search(r'transcript_id\s+\"(\S+)\";', line)
            internal_gene_id_m = re.search(r'gene_id\s+\"(\S+)\";',line)
            internal_txpt_id = ''
            internal_gene_id = ''
            ref_txpt_id = ''
            cls = ''
            line_id= line.strip().split("\t")
            if internal_txpt_id_m and cls_m and internal_gene_id_m:
                cls = cls_m.group(1)
                internal_txpt_id = internal_txpt_id_m.group(1)
                internal_gene_id = internal_gene_id_m.group(1)
                if re.match(r'^[^u]$', cls) and not ref_txpt_id_m:
                    raise Exception(
                        '{} of combined gtf file has logical problem: class code 为非U得情况没有cmp_ref的值 '.format(line))
                if re.match(r'^[u]$', cls) and ref_txpt_id_m:
                    raise Exception(
                        '{} of combined gtf file has logical problem: class code 为U得情况有nearest_ref的值 '.format(line))
                if ref_txpt_id_m:
                    ref_txpt_id = ref_txpt_id_m.group(1)
            else:
                raise Exception(
                    'line: {} in annotate gtf {} has no valid internal txpt id or class code'.format(line.strip(),
                                                                                                     cmp_file))
            d[internal_txpt_id] = {'ref_txpt_id': ref_txpt_id, 'cls': cls, 'internal_gene_id': internal_gene_id,"chr":line_id[0],
                                   "strand":line_id[6],"start":line_id[3],"end":line_id[4]}
    return d


def get_gname_gid_dic_from_ref_gtf(ref_gtf_content):
    d = {}
    for line in ref_gtf_content:
        gene_id_m = re.search(r'gene_id\s+\"(\S+)\"', line.strip())
        gname_m = re.search(r'gene_name\s+\"(\S+)\"', line.strip())
        txpt_id_m = re.search(r'transcript_id\s+\"(\S+)\"', line.strip())
        if gene_id_m and gname_m and txpt_id_m:
            gname = gname_m.group(1)
            gene_id = gene_id_m.group(1)
            txpt_id = txpt_id_m.group(1)
            d[txpt_id] = {'gene_id': gene_id, 'gname': gname}
        else:
            pass
            #raise Exception('ref gtf文件不合法的第九列: {},没有完整的gene id 和gene_name 记录'.format(line.strip()))
    return d


def write_text(ref_dic, cmp_dic, f):
    fw = open(f, 'w')
    head = '#{}\n'.format('\t'.join(['assemble_txpt_id', 'assemble_gene_id', 'class_code', 'ref_txpt_id', 'ref_gene_id', 'ref_gene_name',"chr","strand","start","end"]))
    fw.write(head)
    
    for internal_id in cmp_dic.keys():
        cls = cmp_dic[internal_id]['cls']
        internal_gene_id = cmp_dic[internal_id]['internal_gene_id']
        ref_tid = cmp_dic[internal_id]['ref_txpt_id']
        chr = cmp_dic[internal_id]["chr"]
        strand = cmp_dic[internal_id]["strand"]
        start = cmp_dic[internal_id]["start"]
        end = cmp_dic[internal_id]["end"]

        if ref_tid:
            try:
               if ref_tid in ref_dic.keys():   #add 1 line by khl 20170412
                   ref_gid = ref_dic[ref_tid]['gene_id']
                   ref_gname = ref_dic[ref_tid]['gname']
               else:
                   ref_gid = ref_gname = '-'
            except Exception:
                print ref_tid
        else:
               ref_tid = ref_gname = ref_gid = '-'

        newline = '\t'.join([internal_id, internal_gene_id, cls, ref_tid, ref_gid, ref_gname,chr,strand,start,end]) + '\n'
        fw.write(newline)
    
    fw.close()


def write_relation_text_for_merge_gtf(ref_gtf, cmp_gtf, text_file):
    ref_tmp_cmd = """ awk -F '\t' 'NF>=9{print $9}' %s  | uniq  """ % (ref_gtf)
    ref_tmp_content = subprocess.check_output(ref_tmp_cmd, shell=True).strip().split('\n')
    ref_t_gid_gname_dic = get_gname_gid_dic_from_ref_gtf(ref_tmp_content)
    cmp_tid_ref_tid_cls_dic = get_dic_from_cmp_gtf(cmp_gtf)
    write_text(ref_t_gid_gname_dic, cmp_tid_ref_tid_cls_dic, text_file)


if __name__ == '__main__':
    #ref = '/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/ref.gtf'
    #cmp_f = '/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/cuffmerge_gffcompare/gffcmp_cuffmerge_.annotated.gtf'
    
    #out = '/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/cuffmerge_gffcompare/relation_table.tsv'
    #ref = "/mnt/ilustre/users/sanger-dev/sg-users/qindanhua/test_files/ref_rna/pipe/ref/ref_genome.gtf" 
    #cmp_f = "/mnt/ilustre/users/sanger-dev/workspace/20170508/Single_assembly_module_tophat_stringtie_true_file/RefrnaAssemble/output/Gffcompare/cuffcmp.annotated.gtf"
    #out = "/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tofiles/newdata_class_code"
    ref = "/mnt/ilustre/users/sanger-dev/sg-users/qindanhua/test_files/ref_rna/pipe/ref/ref_genome.gtf"
  
    cmp_f= "/mnt/ilustre/users/sanger-dev/workspace/20170508/Single_assembly_module_tophat_stringtie_true_file/RefrnaAssemble/output/Gffcompare/cuffcmp.annotated.gtf"
    out = "/mnt/ilustre/users/sanger-dev/workspace/20170523/Single_merge_rsem_fpkm_11/MergeRsem/new_class_code"
    
    write_relation_text_for_merge_gtf(ref, cmp_f, out)
    #ref = "/mnt/ilustre/users/sanger-dev/workspace/20170210/Refrna_refrna_test_01/FilecheckRef/Danio_rerio.GRCz10.85.gff3.gtf"
    #cmp_f = "/mnt/ilustre/users/sanger-dev/workspace/20170411/Single_assembly_module_tophat_cufflinks/Assembly/Gffcompare/output/cuffcmp.annotated.gtf"
    #out = "/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/module/cufflinks_class_code/class_code"
    #write_relation_text_for_merge_gtf("cufflinks_cmp_1000.gtf", "ref_1000.gtf", "class_code")

# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import os
import re
import datetime
from bson.son import SON
import types
import gridfs
import json
import unittest
from biocluster.config import Config
import pandas as pd


class RefAnnotation():
    def __init__(self):
        self.result_dir = ''
        self.result_file = {}
        self.trans_gene = {}
        self.trans_isgene = {}
        self.has_new = True
        self.anno_type = 'origin'
        self.new_gene_set = set()
        self.new_trans_set = set()
        self.known_gene_set = set()
        self.known_trans_set = set()
        self.annot_dir = ""
        self.out_dir = ""


    def get_trans2gene(self, trans2gene, trans2gene_ref=None):
        if trans2gene and os.path.exists(trans2gene):
            pass
        elif trans2gene_ref and os.path.exists(trans2gene_ref):
            pass
        else:
            raise Exception('转录本基因对应的结果文件{}不存在，请检查'.format(trans2gene))
        print("读入基因转录本对应关系文件 {}".format(trans2gene))
        if trans2gene:
            with open(trans2gene, 'rb') as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip().split("\t")
                    self.trans_gene[line[0]] = line[1]
                    self.new_trans_set.add(line[0])
                    if line[2] == "yes":
                        self.trans_isgene[line[0]] = True
                        self.new_gene_set.add(line[1])
                    else:
                        self.trans_isgene[line[0]] = False
        if trans2gene_ref:
            with open(trans2gene_ref, 'rb') as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip().split("\t")
                    self.trans_gene[line[0]] = line[1]
                    self.known_trans_set.add(line[0])
                    if line[2] == "yes":
                        self.known_gene_set.add(line[1])
                        self.trans_isgene[line[0]] = True
                    else:
                        self.trans_isgene[line[0]] = False

        print("读入基因转录本对应关系文件结束")

    def get_exp_set(self, gene_exp, trans_exp):
        exp_g_set = set()
        with open(gene_exp, 'r') as gene_exp_f:
            gene_exp_f.readline()
            for line in gene_exp_f:
                exps = line.strip().split("\t")[1:]
                if sum(map(float, exps)) != 0:
                    exp_g_set.add(line.strip().split("\t")[0])
                else:
                    pass
        exp_t_set = set()
        with open(trans_exp, 'r') as trans_exp_f:
            trans_exp_f.readline()
            for line in trans_exp_f:
                exps = line.strip().split("\t")[1:]
                if sum(map(float, exps)) != 0:
                    exp_t_set.add(line.strip().split("\t")[0])
                else:
                    pass
        self.all_exp_g = exp_g_set
        self.all_exp_t = exp_t_set


    def run(self, annotation_mudule_dir, trans2gene, trans2gene_ref, exp_level='transcript', gene_exp = None, trans_exp = None, output_dir = ""):
        self.get_trans2gene(trans2gene, trans2gene_ref)
        self.get_exp_set(gene_exp, trans_exp)
        self.annot_dir = annotation_mudule_dir
        self.out_dir = output_dir
        for i in ['AnnotDetail', 'AnnotStatistics']:
            if os.path.exists(os.path.join(output_dir, i)):
                pass
            else:
                os.makedirs(os.path.join(output_dir, i))
        self.get_detail_result()
        self.get_stat_result()


    def get_detail_result(self):
        self.change_annot_result(os.path.join(self.annot_dir, "refannot_class/all_annot.xls"), os.path.join(self.out_dir, 'AnnotDetail/ref_gene_annot_detail.xls'), os.path.join(self.out_dir, 'AnnotDetail/ref_transcript_annot_detail.xls'))
        if os.path.exists(os.path.join(self.annot_dir, "newannot_class/all_annot.xls")):
            self.change_annot_result(os.path.join(self.annot_dir, "newannot_class/all_annot.xls"), os.path.join(self.out_dir, 'AnnotDetail/new_gene_annot_detail.xls'), os.path.join(self.out_dir, 'AnnotDetail/new_transcript_annot_detail.xls'))
            self.change_annot_result(os.path.join(self.annot_dir, "allannot_class/all_annot.xls"), os.path.join(self.out_dir, 'AnnotDetail/all_gene_annot_detail.xls'), os.path.join(self.out_dir, 'AnnotDetail/all_transcript_annot_detail.xls'))

    def change_annot_result(self, annot, gene_out, trans_out):
        '''
        根据页面展示修改详情表
        '''
        annot_df = pd.read_table(annot, header=0)
        annot_df = annot_df.fillna("")
        go_dict = {
            'cellular_component': "CC",
            'molecular_function': "MF",
            'biological_process': "BP"
        }
        annot_df['GO_id'] =annot_df['go'].map(lambda x: ";".join([a.split("(")[0] for a in x.split("; ")]))
        annot_df['GO_term'] =annot_df['go'].map(lambda x: ";".join([a.split("(")[1].split(":")[0]
                                                               for a in x.split("; ") if len(a.split("(")) > 1]))
        annot_df['Go_description'] =annot_df['go'].map(lambda x: ";".join([go_dict.get(a.split("(")[1].split(":")[0], 'UN') + ":" +
                                                               a.split("(")[1].split(":")[1].rstrip(")")
                                                               for a in x.split("; ") if len(a.split("(")) > 1]))
        annot_df['COG_id'] = annot_df['cog'].map(lambda x: ";".join([a.split("(")[0] for a in x.split("; ")]))
        annot_df['COG_Functional_Categories'] = annot_df['cog'].map(lambda x: ";".join([a.split("(")[1].rstrip(')') for a in x.split("; ") if len(a.split("(")) > 1]))
        annot_df['NR_hit-name'] = annot_df['nr'].map(lambda x: x.split("(")[0])
        annot_df['NR_description'] = annot_df['nr'].map(lambda x: x.split("(")[1].rstrip(")") if len(x.split("(")) > 1 else "")
        annot_df['Swiss-Prot_hit-name'] = annot_df['swissprot'].map(lambda x: x.split("(")[0])
        annot_df['Swiss-Prot_description'] = annot_df['swissprot'].map(lambda x: x.split("(")[1].rstrip(")") if len(x.split("(")) > 1 else "")

        annot_df['Pathway_id'] =annot_df['paths'].map(lambda x: ";".join([a.split("(")[0] for a in x.split("; ")]))
        annot_df['Pathway_definition'] =annot_df['paths'].map(lambda x: ";".join([a.split("(")[1].rstrip(")")
                                                               for a in x.split("; ") if len(a.split("(")) > 1]))
        annot_df['Pfam_id'] =annot_df['pfam'].map(lambda x: ";".join([a.split("(")[0] for a in x.split("; ")]))
        annot_df['Domain'] =annot_df['pfam'].map(lambda x: ";".join([a.split("(")[1].split(":")[0]
                                                               for a in x.split("; ") if len(a.split("(")) > 1]))
        annot_df['Domain_description'] =annot_df['pfam'].map(lambda x: ";".join([
                                                               a.split("(")[1].split(":")[1].rstrip(")")
                                                               for a in x.split("; ") if len(a.split("(")) > 1]))
        annot_df = annot_df.drop(columns=['pfam', 'go', 'cog', 'swissprot', 'nr', 'swissprot'])
        annot_df = annot_df.rename(index=str, columns={"gene_id": "Gene ID", "transcript": "Transcript ID"})
        annot_df_gene = annot_df.loc[annot_df['is_gene'] == "yes"]
        annot_df = annot_df.reindex(["Transcript ID", "Gene ID" ,'length','GO_id', "GO_term", 'Go_description',
                                     "KO_id", "KO_name", "Pathway_id", "Pathway_definition",
                                     'COG_id', 'COG_Functional_Categories',
                                     'NR_hit-name',  'NR_description', 'Swiss-Prot_hit-name', 'Swiss-Prot_description',
                                     'Pfam_id', 'Domain', 'Domain_description'
                                     ], axis=1)
        annot_df_gene = annot_df_gene.reindex(["Gene ID",'length', 'GO_id', "GO_term", 'Go_description',
                                     "KO_id", "KO_name", "Pathway_id", "Pathway_definition",
                                     'COG_id', 'COG_Functional_Categories',
                                     'NR_hit-name',  'NR_description', 'Swiss-Prot_hit-name', 'Swiss-Prot_description',
                                     'Pfam_id', 'Domain', 'Domain_description'
                                     ], axis=1)

        annot_df.to_csv(trans_out, sep="\t", header=True, index=False)
        annot_df_gene.to_csv(gene_out, sep="\t", header=True, index=False)


    def get_stat_result(self):
        self.add_annotation_stat_detail(stat_path=self.annot_dir + '/refannot_class/all_stat.xls',
                                        venn_path=self.annot_dir + '/refannot_class/',
                                        seq_type = "ref",
                                        exp_level = "transcript",
                                        out_file = self.out_dir + '/AnnotStatistics/ref_annot_statistics.xls')
        if os.path.exists(self.annot_dir + '/newannot_class/all_stat.xls'):
            self.add_annotation_stat_detail(stat_path=self.annot_dir + '/newannot_class/all_stat.xls',
                                            venn_path=self.annot_dir + '/newannot_class/',
                                            seq_type = "new",
                                            exp_level = "transcript",
                                            out_file = self.out_dir + '/AnnotStatistics/new_annot_statistics.xls')
            self.add_annotation_stat_detail(stat_path=self.annot_dir + '/allannot_class/all_stat.xls',
                                        venn_path=self.annot_dir + '/allannot_class/',
                                        seq_type = "all",
                                        exp_level = "transcript",
                                        out_file = self.out_dir + '/AnnotStatistics/all_annot_statistics.xls')


    def add_annotation_stat_detail(self, stat_path, venn_path, seq_type, exp_level, out_file):
        """
        database: 进行统计的数据库
        stat_path: all_annotation_statistics.xls
        venn_path: venn图目录
        """
        if not os.path.exists(stat_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(stat_path))
        if not os.path.exists(venn_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(venn_path))


        if seq_type == "ref":
            exp_t_set = self.all_exp_t & self.known_trans_set
            exp_g_set = self.all_exp_g & self.known_gene_set
        elif seq_type == "new":
            exp_t_set = self.all_exp_t & self.new_trans_set
            exp_g_set = self.all_exp_g & self.new_gene_set
        else:
            exp_t_set = self.all_exp_t
            exp_g_set = self.all_exp_g


        gene_dict = dict()
        trans_dict = dict()

        data_list = []
        tail = "gene"

        database_venn = {
            'NR': 'nr/nr_venn_{}.txt'.format(tail),
            'Swiss-Prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Swiss-prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Pfam': 'pfam/pfam_venn_{}.txt'.format(tail),
            'KEGG': 'kegg/kegg_venn_{}.txt'.format(tail),
            'GO': 'go/go_venn_{}.txt'.format(tail),
            'COG': 'cog/cog_venn_{}.txt'.format(tail),
        }
        tail = "tran"
        database_venn_tran = {
            'NR': 'nr/nr_venn_{}.txt'.format(tail),
            'Swiss-Prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Swiss-prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Pfam': 'pfam/pfam_venn_{}.txt'.format(tail),
            'KEGG': 'kegg/kegg_venn_{}.txt'.format(tail),
            'GO': 'go/go_venn_{}.txt'.format(tail),
            'COG': 'cog/cog_venn_{}.txt'.format(tail),
        }
        v3_db2db = {
            "pfam": 'Pfam',
            "kegg": "KEGG",
            "swissprot": "Swiss-Prot",
            "cog": "COG",
            "go": "GO",
            "nr": "NR",
            "annotation": "Total_anno",
            "total": "Total",
            "Pfam": 'Pfam',
            "KEGG": "KEGG",
            "Swiss-Prot": "Swiss-Prot",
            "COG": "COG",
            "GO": "GO",
            "NR": "NR",
            "Total_anno": "Total_anno",
            "Total": "Total"
        }
        with open(stat_path, 'r') as f, open(out_file, 'w') as f_o:
            lines = f.readlines()
            anno_gene_set = set()
            anno_tran_set = set()

            f_o.write("\tExpre_Gene number(percent)\tExpre_Transcript number(percent)\tAll_Gene number(percent)\tAll_Transcript number(percent)\n")
            for line in lines[1:]:
                line = line.strip().split('\t')
                if exp_level.lower() == "gene":
                    data = [
                        str(line[2]) + "({})".format(round(float(line[4]), 4))
                    ]
                else:
                    data = [
                        str(line[2]) + "({})".format(round(float(line[4]), 4)),
                        str(line[1]) + "({})".format(round(float(line[3]), 4)),
                    ]
                venn_list, gene_venn_list = None, None
                database = ["nr", "swissprot", "pfam", "kegg", "go", "string", "cog"]

                if database_venn.has_key(v3_db2db[line[0]]):
                    db = v3_db2db[line[0]]
                    tran_venn = venn_path + "/" + database_venn_tran[v3_db2db[line[0]]]
                    gene_venn = venn_path + "/" + database_venn[v3_db2db[line[0]]]
                    if os.path.exists(tran_venn) and os.path.exists(gene_venn):
                        with open(tran_venn, "rb") as f:
                            tran_venn_set = set([tran.strip() for tran in f])
                        with open(gene_venn, "rb") as f:
                            gene_venn_set = set([gene.strip() for gene in f ])
                        tran_set = exp_t_set & tran_venn_set
                        gene_set = exp_g_set & gene_venn_set
                        anno_gene_set = anno_gene_set | gene_set
                        #print gene_set
                        #print "1"
                        # print anno_gene_set
                        anno_tran_set = anno_tran_set | tran_set
                        data = [
                            str(len(gene_set)) + "({})".format(round(float(len(gene_set))/float(len(exp_g_set)), 4)),
                            str(len(tran_set)) + "({})".format(round(float(len(tran_set))/float(len(exp_t_set)), 4))
                        ] + data
                    else:
                        raise Exception("{}对应的venn.txt文件{} {}不存在".format(line[0], venn, gene_venn))
                elif v3_db2db[line[0]] == "Total_anno":
                    data = [
                        str(len(anno_gene_set)) + "({})".format(round(float(len(anno_gene_set))/float(len(exp_g_set)), 4)),
                        str(len(anno_tran_set)) + "({})".format(round(float(len(anno_tran_set))/float(len(exp_t_set)), 4))
                    ] + data
                elif v3_db2db[line[0]] == "Total":
                    data = [
                        str(len(exp_g_set)) + "({})".format(1),
                        str(len(exp_t_set)) + "({})".format(1)
                    ] + data
                f_o.write(v3_db2db[line[0]]+ "\t" + "\t".join(data) + "\n")


if __name__ == '__main__':
    from mbio.packages.denovo_rna_v2.extract_annotation_result import RefAnnotation
    annot_file = RefAnnotation()
    merge_annot_v3 = "/mnt/ilustre/users/sanger-dev/workspace/20190702/Refrna_tsg_34677/AnnotMerge/output"
    gene_exp = "/mnt/ilustre/users/sanger-dev/workspace/20190702/Refrna_tsg_34677/Quant/gene.tpm.matrix"
    trans_exp = "/mnt/ilustre/users/sanger-dev/workspace/20190702/Refrna_tsg_34677/Quant/transcript.tpm.matrix"
    gene2trans = "/mnt/ilustre/users/sanger-dev/workspace/20190702/Refrna_tsg_34677/AnnotMerge/output/newannot_class/all_tran2gene.txt"
    gene2trans_ref = "/mnt/ilustre/users/sanger-dev/workspace/20190702/Refrna_tsg_34677/AnnotMerge/output/refannot_class/all_tran2gene.txt"
    pp1="/mnt/ilustre/users/sanger-dev/workspace/20191016/Denovorna_tsg_35822/AnnotClassBeta/output/all_annot.xls"
    out_gene = "/mnt/ilustre/users/sanger-dev/workspace/20191016/Denovorna_tsg_35822/AnnotClassBeta/output/test/test_gene"
    out_trans="/mnt/ilustre/users/sanger-dev/workspace/20191016/Denovorna_tsg_35822/AnnotClassBeta/output/test/test_trans"
    annot_file.change_annot_result(pp1,out_gene,out_trans)
    # annot_file.run(merge_annot_v3, gene2trans, gene2trans_ref, gene_exp=gene_exp, trans_exp=trans_exp, output_dir=out_dir)

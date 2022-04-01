# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20180107

import sys
import argparse
import os
import pandas as pd


#from biocluster.config import Config


class OverviewMerge(object):
    def __init__(self):
        self.nr_level = ["#Query", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        self.cog_level = ["#Query", "NOG", "NOG_description", "Function", "Fun_description", "Category"]
        self.kegg_level = ["#Query", "Gene", "KO", "Definition", "Pathway", "Level3", "Enzyme", "Module", "Level1",
                           "Level2"]
        self.cazy_level = ["#Query", "Family", "Class", "Class_description"]
        self.card_level = ["#Query", "ARO", "ARO_name", "ARO_description", "AMR_Gene_Family", "Drug_Class",
                           "Antibiotic_class", "Resistance_Mechanism"]
        self.card_level2 = ["#Query", "ARO", "ARO_name", "ARO_description", "ARO_category", "Class"]
        self.ardb_level = ["#Query", "ARG", "Type", "Antibiotic_type", "Class", "Class_description", "Antibiotic_class"]
        self.ardb_level2 = ["#Query", "ARG", "Type", "Antibiotic_type", "Class", "Class_description"]
        self.vfdb_level = ["#Query", "VFs", "VF_Function", "Species", "Level1", "Level2"]
        self.level_dic = self.get_level_dic()  #guanqing.zou
        self.new_columns_dic = self.get_new_columns_dic()

    def run_merge(self, gene_length, outfile, nr=None, cog= None, kegg=None, cazy=None, ardb = None, card=None,
                  vfdb= None):
        gene_length_table = pd.read_table(gene_length, sep='\t', header=0)
        gene_length_table.index = gene_length_table["GeneID"]
        merge_file = gene_length_table
        merge_file.columns = ["#Query", "Length"]
        if nr != None:
            nr_table = pd.read_table(nr, sep='\t', header=0)
            nr_table = nr_table[self.nr_level]
            #nr_table = nr_table[nr_table["#Query"].isin(gene_length_table.index)]
            merge_file = merge_file.merge(nr_table, left_on='#Query', right_on='#Query', how='outer')
        if cog != None:
            cog_table = pd.read_table(cog, sep='\t', header=0)
            cog_table = cog_table[self.cog_level]
            #cog_table = cog_table[cog_table["#Query"].isin(gene_length_table.index)]
            cog_table.columns = ["#Query", "NOG", "COG_description", "COG_Function", "COG_Fun_description",
                                 "COG_Category"]
            merge_file = merge_file.merge(cog_table, left_on='#Query', right_on='#Query', how='outer')
        if kegg != None:
            kegg_table = pd.read_table(kegg, sep='\t', header=0)
            kegg_table = kegg_table[self.kegg_level]
            #kegg_table = kegg_table[kegg_table["#Query"].isin(gene_length_table.index)]
            kegg_table.columns = ["#Query", "KEGG_Gene", "KO", "KO_description", "KEGG_Pathway",
                                "KEGG_Path_description",
                                "KEGG_Enzyme", "KEGG_Modules", "KEGG_Level1", "KEGG_Level2"]
            enzyme_profile = "/".join(kegg.split("/")[0:len(kegg.split("/")) - 1]) + "/kegg_enzyme_profile.xls"
            module_profile = "/".join(kegg.split("/")[0:len(kegg.split("/")) - 1]) + "/kegg_module_profile.xls"
            if os.path.exists(enzyme_profile):
                enzyme_des = pd.read_table(enzyme_profile, sep='\t', header=0)
                enzyme_des = enzyme_des[["#Enzyme", "Description"]]
                enzyme_des.columns = ["KEGG_Enzyme", "KEGG_En_description"]
                kegg_table = kegg_table.merge(enzyme_des, left_on='KEGG_Enzyme', right_on='KEGG_Enzyme', how='outer')
                #print kegg_table.head()
            else:
                kegg_dir = "/".join(kegg.split("/")[0:len(cazy.split("/") - 1)])
                raise Exception("{}目录下缺少文件kegg_enzyme_profile.xls!".format(kegg_dir))
            if os.path.exists(module_profile):
                module_des = pd.read_table(module_profile, sep='\t', header=0)
                module_des = module_des[["#Module", "Description"]]
                module_des.columns = ["KEGG_Modules", "KEGG_Mo_description"]
                kegg_table = kegg_table.merge(module_des, left_on='KEGG_Modules', right_on='KEGG_Modules', how='outer')
            else:
                kegg_dir = "/".join(cazy.split("/")[0:len(cazy.split("/") - 1)])
                raise Exception("{}目录下缺少文件kegg_enzyme_profile.xls!".format(kegg_dir))
            if "KEGG name" in kegg_table.columns:
                kegg_table = kegg_table[["#Query", "KEGG_Gene", "KO", "KEGG name", "KO_description", "KEGG_Pathway", "KEGG_Path_description", \
                    "KEGG_Enzyme", "KEGG_En_description", "KEGG_Modules", "KEGG_Mo_description", "KEGG_Level1", "KEGG_Level2"]]
            else:
                kegg_table = kegg_table[["#Query", "KEGG_Gene", "KO", "KO_description", "KEGG_Pathway", "KEGG_Path_description", \
                    "KEGG_Enzyme", "KEGG_En_description", "KEGG_Modules", "KEGG_Mo_description", "KEGG_Level1", "KEGG_Level2"]]
            merge_file = merge_file.merge(kegg_table, left_on='#Query', right_on='#Query', how='outer')
        if cazy != None:
            cazy_table = pd.read_table(cazy, sep='\t', header=0)
            cazy_table = cazy_table[self.cazy_level]
            #cazy_table = cazy_table[cazy_table["#Query"].isin(gene_length_table.index)]
            cazy_table.columns = ["#Query", "CAZY_Family", "CAZY_Class", "CAZY_Cl_description"]
            fa_profile = "/".join(cazy.split("/")[0:len(cazy.split("/")) - 1]) + "/cazy_family_profile.xls"
            if os.path.exists(fa_profile):
                fa_des = pd.read_table(fa_profile, sep='\t', header=0)
                fa_des = fa_des[["#Family", "Description"]]
                fa_des.columns = ["CAZY_Family", "CAZY_Fa_description"]
                cazy_table = cazy_table.merge(fa_des, left_on='CAZY_Family', right_on='CAZY_Family', how='outer')
            else:
                cazy_dir = "/".join(cazy.split("/")[0:len(cazy.split("/"))- 1])
                raise Exception("{}目录下缺少文件cazy_family_profile.xls!".format(cazy_dir))
            merge_file = merge_file.merge(cazy_table, left_on='#Query', right_on='#Query', how='outer')
        if ardb != None:
            ardb_table = pd.read_table(ardb, sep='\t', header=0)
            try: ##先做新的数据库
                ardb_table = ardb_table[self.ardb_level]
                #ardb_table = ardb_table[ardb_table["#Query"].isin(gene_length_table.index)]
                ardb_table.columns = ["#Query", "ARDB_ARG", "ARDB_Type", "ARDB_Antibiotic_type", "ARDB_Class","ARDB_Cl_description", "ARDB_Antibiotic_class"]
            except:## 然后兼容老的数据库版本 fix by qingchen.zhang@20210203
                ardb_table = ardb_table[self.ardb_level2]
                ardb_table.columns = ["#Query", "ARDB_ARG", "ARDB_Type", "ARDB_Antibiotic_type", "ARDB_Class","ARDB_Cl_description"]
            merge_file = merge_file.merge(ardb_table, left_on='#Query', right_on='#Query', how='outer')
        if card != None:
            card_table = pd.read_table(card, sep='\t', header=0)
            try:## 先做新的数据库
                card_table = card_table[self.card_level]
                #card_table = card_table[card_table["#Query"].isin(gene_length_table.index)]
                card_table.columns = ["#Query", "CARD_ARO", "CARD_ARO_name", "CARD_ARO_description", "CARD_AMR_Gene_Family",
                                  "CARD_Drug_Class", "CARD_Antibiotic_class", "CARD_Resistance_Mechanism"]
            except:## 兼容老的数据库版本 fix by qingchen.zhang@20210203
                card_table = card_table[self.card_level2]
                card_table.columns = ["#Query", "CARD_ARO", "CARD_ARO_name", "CARD_ARO_description", "CARD_ARO_category",
                                  "CARD_Class"]
            merge_file = merge_file.merge(card_table, left_on='#Query', right_on='#Query', how='outer')
        if vfdb != None:
            vfdb_table = pd.read_table(vfdb, sep='\t', header=0)
            vfdb_table = vfdb_table[self.vfdb_level]
            #vfdb_table = vfdb_table[vfdb_table["#Query"].isin(gene_length_table.index)]
            vfdb_table.columns = ["#Query", "VFDB_VFs", "VFDB_VF_Function", "VFDB_Species", "VFDB_Level1",
                                  "VFDB_Level2"]
            merge_file = merge_file.merge(vfdb_table, left_on='#Query', right_on='#Query', how='outer')
        merge_file = merge_file[merge_file["#Query"].isin(gene_length_table.index)].drop_duplicates('#Query')
        merge_file["count"] = merge_file.T.count()
        merge_file = merge_file.sort_values(["count"], ascending=False)
        merge_file = merge_file.drop('count', 1)
        merge_file = merge_file.fillna("-")
        merge_file.to_csv(outfile, sep="\t", index=False)

    #guanqing
    def common_fun(self, infile, nead_columns, new_columns, result_table):
        table = pd.read_table(infile, sep='\t', header=0)
        table_sub = table[nead_columns]
        table_sub.columns = new_columns
        result_table = result_table.merge(table_sub, left_on = '#Query', right_on = '#Query', how = 'outer')
        return result_table

    def get_gene_length(self, gene_length):
        gene_length_table = pd.read_table(gene_length, sep='\t', header=0)
        gene_length_table.index = gene_length_table["GeneID"]
        merge_file = gene_length_table
        merge_file.columns = ["#Query", "Length"]
        return merge_file

    def get_level_dic(self):
        levels = {'mvirdb':['#Query', 'status','Short Description','Virulence Factor Type', 'Database Source','Virulence Factor ID'],
                  'tcdb' : ['#Query', 'TCDB Family', 'TCDB Subclass', 'TCDB Class', 'TCDB ID'],
                  'pfam' : ['#Query', 'Domain', 'Type', 'Clan ID', 'Pfam ID'], #'Pfam Accession']
                  'p450' : ['#Query', 'homo_family', 'super_family'],
                  'phi' : ['#Query', 'Pathogen Species', 'Phenotype', 'Host Description', 'Function', 'PHI ID', 'protein', 'Hit Gene Name', 'Experimental Host Species'],
                  'qs' : ['#Query', 'Class'],
                  'probio' : ['#Query', 'Probiotic_name', 'Genus', 'Use_in', 'Disease_class','Commercial_Development_Stage','Probiotic_Effect'],
                  'go' : ['#Query',  'GO (Lev1)', 'GO Term (Lev2)', 'GO Term (Lev3)', 'GO Term (Lev4)'],
                  #'sec' : ['#Query','Gram neg', 'Gram pos', 'Fungi'],
                  'nr_lca' : ["#Query", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"],
                  'nr_de_unclassified' : ["#Query", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"],
                  't3ss' : ['gene_id','score'],
                  'sec_Gram_neg' : ['Gene ID', 'Smean'],
                  'sec_Gram_pos' : ['Gene ID', 'Smean'],
                  'sec_Euk' : ['Gene ID', 'Smean'],
        }
        return levels

    def get_new_columns_dic(self):
        new_columns = {
            'mvirdb':['#Query', 'Mvirdb_status','Mvirdb_Short_Description','Mvirdb_Virulence_Factor_Type', 'Mvirdb_Database_Source','Mvirdb_Virulence_Factor_ID'],
            'tcdb' : ['#Query', 'TCDB_Family', 'TCDB_Subclass', 'TCDB_Class', 'TCDB_ID'],
            'pfam' : ['#Query', 'pfam_Domain', 'pfam_Type', 'pfam_Clan_ID', 'pfam_Pfam_ID'],
            'p450' : ['#Query', 'p450_homo_family', 'p450_super_family'],
            'phi' : ['#Query', 'phi_Pathogen_Species', 'phi_Phenotype', 'phi_HOST', 'phi_Function', 'phi_PHI_ID', 'phi_protein', 'phi_Gene_Name', 'phi_Experimental_Host_Species'],
            'qs' : ['#Query', 'qs_Class'],
            'probio' : ['#Query', 'probio_Probiotic_name', 'probio_Genus', 'probio_Use_in', 'probio_Disease_class','probio_Commercial_Development_Stage','probio_Probiotic_Effect'],
            'go' : ['#Query',  'go_GO_(Lev1)', 'go_GO_Term_(Lev2)', 'go_GO_Term_(Lev3)', 'go_GO_Term_(Lev4)'],
            #'sec' : ['#Query','sec_Gram_neg', 'sec_Gram_pos', 'sec_Fungi'],
            'nr_lca' : ["#Query", "lca_Domain", "lca_Kingdom", "lca_Phylum", "lca_Class", "lca_Order", "lca_Family", "lca_Genus", "lca_Species"],
            'nr_de_unclassified' : ["#Query", "duc_Domain", "duc_Kingdom", "duc_Phylum", "duc_Class", "duc_Order", "duc_Family", "duc_Genus", "duc_Species"],
            't3ss' : ['#Query','t3ss_score'],
            'sec_Gram_neg' : ['#Query', 'sec_Gram_neg'],
            'sec_Gram_pos' : ['#Query', 'sec_Gram_pos'],
            'sec_Euk' : ['#Query', 'sec_Euk'],
        }
        return new_columns


    def run_merge_personal(self,types,files, result_table):
        database_list = types.split(',')
        files_list = files.split(',')
        num = len(database_list)
        if num != len(files_list):
            print('len(database_list) != len(files_list)')
            exit()
        for id in range(num):
            database = database_list[id]
            result_table = self.common_fun(files_list[id],self.level_dic[database], self.new_columns_dic[database], result_table)
        if 't3ss_score' in result_table.columns.values:
            result_table['t3ss_score']=result_table['t3ss_score'].apply(lambda x : True if x > 0 else '')

        for i in ['sec_Gram_neg','sec_Gram_pos','sec_Euk']:
            if i  in result_table.columns.values:
                result_table[i] = result_table[i].apply(lambda x: True if x > 0 else '')
        return result_table

    def write_file(self, merge_file, gene_length_table, outfile):
        merge_file = merge_file[merge_file["#Query"].isin(gene_length_table.index)]
        merge_file["count"] = merge_file.T.count()
        merge_file = merge_file.sort_values(["count"], ascending=False)
        merge_file = merge_file.drop('count', 1)
        merge_file = merge_file.fillna("-")
        merge_file.to_csv(outfile, sep="\t", index=False)

    def merge_pre_sum(self, data, sum_pre, outfile):
        table_sub = pd.read_table(sum_pre, sep='\t', header=0)
        #data = data.drop(['Length'],axis = 1)
        #del data['Length']
        table_sub_columns = table_sub.columns.values
        data_columns = data.columns.values
        #data_columns.remove('#Query')

        for c in data_columns:
            if c == '#Query':
                continue
            if c in table_sub_columns:
                del table_sub[c]
        result_table = data.merge(table_sub, left_on = '#Query', right_on = '#Query', how = 'outer')
        result_table["count"] = result_table.T.count()
        result_table = result_table.sort_values(["count"], ascending=False)
        result_table = result_table.drop('count', 1)
        result_table = result_table.fillna("-")
        result_table.to_csv(outfile, sep="\t", index=False)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-gl', metavar='[gene_length_table]', required=True, help='Input gene_length_table')
    parser.add_argument('-o', metavar='[output file]', required=True, help='output file name')
    parser.add_argument('-nr', metavar='[gene_nr_anno]', help='Input gene_nr_anno.xls file')
    parser.add_argument('-cog', metavar='[gene_cog_anno]', help='Input gene_cog_anno.xls file')
    parser.add_argument('-kegg', metavar='[gene_kegg_anno]', help='Input gene_kegg_anno.xls file')
    parser.add_argument('-cazy', metavar='[gene_cazy_anno]', help='Input gene_cazy_anno.xls file')
    parser.add_argument('-ardb', metavar='[gene_ardb_anno]', help='Input gene_ardb_anno.xls file')
    parser.add_argument('-card', metavar='[gene_card_anno]', help='Input gene_card_anno.xls file')
    parser.add_argument('-vfdb', metavar='[gene_vfdb_anno]', help='Input gene_vfdb_anno.xls file')
    ##zouguanqing >>
    parser.add_argument('-other')
    parser.add_argument('-otherf')
    parser.add_argument('-add')
    #parser.add_argument('-addo')
    ###
    args = parser.parse_args()
    gene_length = args.gl
    outfile = args.o
    run_overview_merge = OverviewMerge()
    nr = cog = kegg = card = cazy = ardb = vfdb = None
    if args.nr:
        nr = args.nr
    if args.cog:
        cog = args.cog
    if args.kegg:
        kegg = args.kegg
    if args.cazy:
        cazy = args.cazy
    if args.ardb:
        ardb = args.ardb
    if args.card:
        card = args.card
    if args.vfdb:
        vfdb = args.vfdb
    #zouguaniqing >>>
    if not args.other and not args.add:    #标准注释表汇总
        run_overview_merge.run_merge(gene_length, outfile, nr=nr, cog=cog, kegg=kegg, cazy=cazy, ardb=ardb, card=card,
                                 vfdb=vfdb)
    else:    #个性化注释汇总
        if args.other == 'None':
            pass
        else:
            length_table = run_overview_merge.get_gene_length(gene_length)
            merge_data = run_overview_merge.run_merge_personal(args.other, args.otherf, length_table)
            merge_data["count"] = merge_data.T.count()
            merge_data = merge_data.sort_values(["count"], ascending=False)
            merge_data = merge_data.drop('count', 1)
            merge_data = merge_data.fillna("-")
            run_overview_merge.write_file(merge_data, length_table, outfile)
            if args.add :          #2张汇总表汇总
                #run_overview_merge.merge_pre_sum(merge_data, args.add, args.addo)
                run_overview_merge.merge_pre_sum(merge_data, args.add,outfile)




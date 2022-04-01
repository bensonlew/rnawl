# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
#20190916

import os,re
import pandas as pd
import argparse


class CommonAnno():

    def add_merge(self, dir, type, out, sample_add=None, add_del=None):
        """
        将注释目录下的文件相同type的合并
        :param dir: 注释总目录
        :param type: 数据库类型，对应就是不同数据库合并，相同数据库就给定不同的字段
        :param out: 输出文件
        :param sample_add: 是否在文件中加入样品名称
        :return:
        """
        files = os.listdir(dir)
        n = 1
        if os.path.exists(out):
            os.remove(out)
        for file in files:
            aa = os.listdir(dir + "/" + file)
            for a in aa:
                if re.search("{}".format(type), a):
                    a = pd.read_table(dir + "/" + file + "/" + a, sep='\t', header=0, dtype={'fullVisitorId': 'str'})
                    if sample_add in ['true', 'True', "TRUE"]:
                        a['sample'] = file
                    elif sample_add in ['false', 'False', "FALSE"]:
                        pass
                    if add_del:
                        a = a.drop([add_del], axis=1)
                    if n == 1:
                        a.to_csv(out, mode='a', sep='\t', header=True, index=False)
                    elif n > 1:
                        a.to_csv(out, mode='a', sep='\t', header=0, index=False)
                    n += 1


    def anno_abund(self, file, func, gene_id, out, add=None):
        """
        计算某一水平下基因个数
        :param file: 所有样品对应的注释文件
        :param func: 功能水平对应列名
        :param gene_id: 基因的列名，主要针对不同的注释结果名称
        :param out: 结果所有样品对应不同功能的基因个数，没有填充为0
        :return:
        """
        a = pd.read_table(file, sep='\t', header=0, dtype={'fullVisitorId': 'str'})
        a = a.groupby([func, 'sample'])[gene_id].count().unstack()
        a = a.fillna("0.0")
        #if add in ["Yes", "yes"]:
            #a = a.drop(["-"], axis=0)
        a.to_csv(out, sep='\t', header=True, index=True)

    def anno_genelist(self, file, fun, gene_id, func_levels, out):
        """
        得到某一个数据库最小功能水平的gene_list，不存在的“-”
        :param file: 所有样品对应的注释文件
        :param fun: 最低功能水平的列名
        :param gene_id: 基因的列名，主要针对不同的注释结果名称
        :param func_levels: 最小功能水平对应上级的功能
        :param out: 输出的所有样品对应功能的gene_list文件
        :return:
        """
        anno = pd.read_table(file, sep='\t', header=0, dtype={'fullVisitorId': 'str'})
        samples = set(anno.ix[:, 'sample'])
        aa = dict(list(anno.groupby([fun, 'sample'])))
        anno = {}
        for key in aa.keys():
            if len(list(aa[key][gene_id])) >= 1:
                des = ','.join(list(aa[key][gene_id]))
            else:
                des = "-"
            if key[0] not in anno:
                anno[key[0]] = {key[1]: des}
            else:
                anno[key[0]][key[1]] = des
        list2 = (func_levels).split(",")
        def anno_des(file, list, func):
            dict = {}
            with open(file, "r") as f:
                lines = f.readlines()
                names = lines[0].strip('\n').split("\t")
                for line in lines[1:]:
                    lin = line.strip().split("\t")
                    arry = []
                    for i in list:
                        arry.append(lin[names.index(i)])
                    fun = lin[names.index(func)]
                    dict[fun] = "\t".join(arry)
            return dict

        ss = ''
        func = anno_des(file, list2, fun)
        for key in sorted(anno):
            ss += key + "\t" + func[key]
            for sample in samples:
                if sample in anno[key]:
                    ss += "\t" + anno[key][sample]
                else:
                    ss += "\t-"
            ss += "\n"
        if os.path.exists(out):
            os.remove(out)
        with open(out, "w") as f:
            f.write(fun + "\t" + "\t".join(list2) + "\t" + '\t'.join(samples) + "\n")
            f.write(ss)

    def anno_ano_genelist(self, file, fun, gene_id, out):
        """
        得到某一个数据库最小功能水平的gene_list，不存在的“-”
        :param file: 所有样品对应的注释文件
        :param fun: 最低功能水平的列名
        :param gene_id: 基因的列名，主要针对不同的注释结果名称
        :param out: 输出的所有样品对应功能的gene_list文件
        :return:
        """
        anno = pd.read_table(file, sep='\t', header=0, dtype={'fullVisitorId': 'str'})
        samples = set(anno.ix[:, 'sample'])
        aa = dict(list(anno.groupby([fun, 'sample'])))
        anno = {}
        for key in aa.keys():
            if len(list(aa[key][gene_id])) >= 1:
                des = ','.join(list(aa[key][gene_id]))
            else:
                des = "-"
            if key[0] not in anno:
                anno[key[0]] = {key[1]: des}
            else:
                anno[key[0]][key[1]] = des
        ss= ''
        for key in sorted(anno):
            ss += key
            for sample in samples:
                if sample in anno[key]:
                    ss += "\t" + anno[key][sample]
                else:
                    ss += "\t-"
            ss += "\n"
        if os.path.exists(out):
            os.remove(out)
        with open(out, "w") as f:
            f.write(fun + "\t" + '\t'.join(samples) + "\n")
            f.write(ss)

    def anno_chuli(self, file, func, out):
        with open(file, "r") as f, open(out, "w") as w:
            lines = f.readlines()
            w.write(lines[0])
            names = lines[0].strip('\n').split("\t")
            num = names.index(func)
            for line in lines[1:]:
                lin = line.strip().split("\t")
                if re.search(';', lin[num]):
                    func = lin[num].split(";")
                    for i in func:
                        lin[num] = i
                        w.write("\t".join(lin) + "\n")
                else:
                    w.write("\t".join(lin) + "\n")

    def fill_data(self, file, sample_list, out, data, fun):
        table = pd.read_table(file, sep='\t', header=0)
        name =list(table.columns)
        add_name =list(set(sample_list).difference(set(name)))
        table = pd.concat([table, pd.DataFrame(columns=add_name)])
        table = table.fillna(data)
        last_data = table.set_index(fun)
        last_data.to_csv(out, sep='\t', header=True, index=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='different function datbase in your table ')
    parser.add_argument('-dir', help="annotation dir")
    parser.add_argument('-t', help="the type of database or different of type the same database")
    parser.add_argument('-o', help="out of dir")
    args = parser.parse_args()
    anno =CommonAnno()
    if args.t == "cog":
        anno.add_merge(args.dir, "cog_anno", args.o + "all.cog_anno.xls", sample_add="true")
        anno.anno_chuli(args.o + "all.cog_anno.xls", "Function", args.o + "all.cog_Function_anno.xls")
        anno.anno_abund(args.o + "all.cog_anno.xls", "NOG", "Gene ID", args.o + "all.cog_abundance.xls")
        anno.anno_abund(args.o + "all.cog_Function_anno.xls", "Function", "Gene ID", args.o + "all.Function_abundance.xls")
        anno.anno_genelist(args.o + "all.cog_anno.xls", "NOG", "Gene ID", "NOG_description,Function,Fun_description,Category",args.o +  "all.cog_genelist.xls")
        os.remove(args.o + "all.cog_Function_anno.xls")
    elif args.t == "kegg":
        anno.add_merge(args.dir, "kegg_anno", args.o + "all.kegg_anno.xls", sample_add="true")
        anno.anno_chuli(args.o + "all.kegg_anno.xls", "Level3", args.o + "all.kegg_level3_anno.xls")
        anno.anno_abund(args.o + "all.kegg_anno.xls", "KO", "Gene ID", args.o + "all.KO_abundance.xls")
        anno.anno_abund(args.o + "all.kegg_level3_anno.xls", "Level3", "Gene ID", args.o + "all.level3_abundance.xls", add="yes")
        anno.anno_genelist(args.o + "all.kegg_anno.xls", "KO", "Gene ID", "Gene Name,Definition,Pathway,Enzyme,Module,Hyperlink,Level1,Level2,Level3",args.o +  "all.KEGG_genelist.xls")
        os.remove(args.o + "all.kegg_level3_anno.xls")
    elif args.t == "cazy":
        anno.add_merge(args.dir, "cazy_anno", args.o + "all.cazy_anno.xls", sample_add="true")
        anno.anno_genelist(args.o + "all.cazy_anno.xls", "Family", "Gene ID",
                           "Class,Class_description,Family_description",
                           args.o + "all.cazy_genelist.xls")
    elif args.t == "vfdb":
        anno.add_merge(args.dir, "vfdb_core_anno", args.o + "all.vfdb_core_anno.xls", sample_add="true")
        anno.anno_genelist(args.o + "all.vfdb_core_anno.xls", "VFDB ID", "#Query",
                           "Gi_num,VfdbGene,Gene_description,Species,VFs,VFs_Description,Level2,Level1",
                           args.o + "all.vfdb_core_genelist.xls")
        anno.add_merge(args.dir, "vfdb_predict_anno", args.o + "all.vfdb_predict_anno.xls", sample_add="true")
        anno.anno_genelist(args.o + "all.vfdb_predict_anno.xls", "VFDB ID", "#Query",
                           "Gi_num,VfdbGene,Gene_description,Species,VFs,VFs_Description,Level2,Level1",
                           args.o + "all.vfdb_predict_genelist.xls")
        anno.add_merge(args.dir, "vfdb_anno", args.o + "all.vfdb_anno.xls", sample_add="true")
        anno.anno_genelist(args.o + "all.vfdb_anno.xls", "VFDB ID", "Gene ID",
                           "Gi_num,VfdbGene,Gene_description,Species,VFs,VFs_Description,Level2,Level1",
                           args.o + "all.vfdb_all_genelist.xls")
    elif args.t == "phi":
        anno.add_merge(args.dir, "phi_anno", args.o + "all.phi_anno.xls", sample_add="true")
        anno.anno_genelist(args.o + "all.phi_anno.xls", "PHI ID", "Gene ID",
                           "Protein ID,Hit Gene Name,NCBI Tax ID,Pathogen Species,Phenotype,Host Descripton,Host NCBI ID,Host Species,Gene Function,Disease",
                           args.o + "all.phi_genelist.xls")
    elif args.t == "card":
        anno.add_merge(args.dir, "card_anno", args.o + "all.card_anno.xls", sample_add="true")
        anno.anno_genelist(args.o + "all.card_anno.xls", "ARO_Accession", "Gene ID",
                           "ARO_name,ARO_description,ARO_category,Drug_class,Resistance_mechanism",
                           args.o + "all.card_genelist.xls")
    elif args.t == "tcdb":
        anno.add_merge(args.dir, "tcdb_anno", args.o + "all.tcdb_anno.xls", sample_add="true")
        anno.anno_genelist(args.o + "all.tcdb_anno.xls", "TCDB ID", "Gene ID",
                           "TCDB Description,TCDB Family,TCDB Subclass,TCDB Class",
                           args.o + "all.tcdb_genelist.xls")
    elif args.t == "secretory":
        anno.add_merge(args.dir, "secretory_anno", args.o + "all.secretory_anno.xls", sample_add="true")
        anno.anno_genelist(args.o + "all.secretory_anno.xls", "Type", "Gene ID", "Gene,KO,des", args.o + "all.secretory_genelist.xls")
    elif args.t == "tmhmm":
        anno.add_merge(args.dir, "tmhmm_anno", args.o + "all.tmhmm_anno.xls", sample_add="true")
    elif args.t == "signalp":
        anno.add_merge(args.dir, "Gram--_SignalP", args.o + "all.Gram+_SignalP.xls")
        anno.add_merge(args.dir, "Gram-_SignalP", args.o + "all.Gram-_SignalP.xls")
    elif args.t ==  "all":
        anno.add_merge(args.dir, "anno_summary", args.o + "all.anno_summary.xls")




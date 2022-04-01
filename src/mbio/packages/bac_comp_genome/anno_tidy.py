# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
#20190916

import os
import pandas as pd
import argparse


def add_gene_info(sample, gff, cog, kegg, cazy, card, phi, vfdb, signalp1, signalp2, tmhmm, tcdb, secretory, out):
    gene = pd.read_table(gff, sep='\t', header=0)
    gene = gene.ix[:, ["Gene ID", "Sequence id", "Strand", "Start", "End", "Gene Length(bp)", "Type"]]
    gene["Location"] = gene["Sequence id"].str.split("_ORF", expand=True)[0]
    del gene['Sequence id']

    def cent_fun(values):
        n = 0
        all = 0
        for v in values:
            if v == "-":
                n += 1
            all += 1
        return float(n)/all*100

    gene2 = gene.ix[:, ["Gene ID", "Location"]]
    cog = pd.read_table(cog, sep='\t', header=0)
    cog["Gene ID"] = cog["#Query"]
    del cog["#Query"]
    cog2 = gene2.merge(cog, on='Gene ID', how='left')
    cog2 = cog2.fillna("-")
    num = float(len(cog2.columns) -2)/len(cog2.columns)*100
    cog2 = cog2[cog2.apply(lambda x: cent_fun(x) < num, axis=1)]
    cog2.index = cog2["Gene ID"]
    cog2.to_csv(out + ".cog_anno.xls", sep='\t', header=True, index=False)
    cog = cog.ix[:, ["Gene ID", "NOG", "NOG_description", "Function", "Fun_description"]]
    cog.columns = ["Gene ID", "COG ID", "COG Description", "Function", "Fun Description"]
    kegg = pd.read_table(kegg, sep='\t', header=0)
    kegg["Gene ID"] = kegg["#Query"]
    del kegg["#Query"]
    kegg2 = gene2.merge(kegg, on='Gene ID', how='left')
    kegg2 = kegg2.fillna("-")
    num = float(len(kegg2.columns) - 2) / len(kegg2.columns) * 100
    kegg2 = kegg2[kegg2.apply(lambda x: cent_fun(x) < num, axis=1)]
    kegg_anno2 = dict(list(kegg2.groupby(["Gene ID"])))
    list_kegg2 = []
    kkg2 = pd.DataFrame(columns=["Gene ID","Location","Gene Name","KO","Definition","Pathway","Enzyme","Module","Hyperlink","Level1","Level2","Level3","Identity(%)","Align_len"])
    for key in kegg_anno2.keys():
        if len(kegg_anno2[key]) == 1:
            list_kegg2.append(kegg_anno2[key])
        elif len(kegg_anno2[key]) > 1:
            data = {}
            for j in ["Gene ID", "Location", "Gene Name", "KO", "Definition", "Enzyme", "Module", "Hyperlink","Identity(%)", "Align_len"]:
                data[j] = [";".join(set(kegg_anno2[key][j].astype('str')))]
            dd = []
            for s in set(kegg_anno2[key]["Pathway"].astype('str')):
                dd.append(list(kegg_anno2[key]["Pathway"].astype('str')).index(s))
            for j in ["Pathway", "Level1", "Level2", "Level3"]:
                lists = []
                for d in dd:
                    lists.append(list(kegg_anno2[key][j].astype('str'))[d])
                data[j] = [";".join(lists)]
            kkg2 = kkg2.append(pd.DataFrame(data), ignore_index=True)
    kegg2 = kkg2.append(list_kegg2)
    kegg2.index = kegg2["Gene ID"]
    kegg2.to_csv(out + ".kegg_anno.xls", columns=["Gene ID","Location","Gene Name","KO","Definition","Pathway","Enzyme","Module","Hyperlink","Level1","Level2","Level3","Identity(%)","Align_len"], sep='\t', header=True, index=False)
    kegg = kegg.ix[:, ["Gene ID", "Gene Name", "KO", "Definition"]]
    kegg.columns = ["Gene ID", "Gene Name", "KO ID", "KO Description"]
    kegg_anno= dict(list(kegg.groupby(["Gene ID"])))
    list_kegg = []
    kkg = pd.DataFrame(columns=["Gene ID", "Gene Name", "KO ID", "KO Description"])
    for key in kegg_anno.keys():
        if len(kegg_anno[key]) ==1:
            list_kegg.append(kegg_anno[key])
        elif len(kegg_anno[key]) >1:
            data ={}
            for j in ["Gene ID","Gene Name", "KO ID", "KO Description"]:
                data[j] = [";".join(set(kegg_anno[key][j]))]
            kkg = kkg.append(pd.DataFrame(data), ignore_index=True)
    kegg = kkg.append(list_kegg)
    cazy = pd.read_table(cazy, sep='\t', header=0)
    cazy["Gene ID"] = cazy["#Query"]
    del cazy["#Query"]
    cazy2 = gene2.merge(cazy, on='Gene ID', how='left')
    cazy2 = cazy2.fillna("-")
    num = float(len(cazy2.columns) - 2)/ len(cazy2.columns) * 100
    cazy2 = cazy2[cazy2.apply(lambda x: cent_fun(x) < num, axis=1)]
    cazy2.index = cazy2["Gene ID"]
    cazy2.to_csv(out + ".cazy_anno.xls", sep='\t', header=True, index=False)
    cazy = cazy.ix[:, ["Gene ID", "Family","Family_description", "Class", "Class_description"]]
    cazy.columns = ["Gene ID", "Family", "Family Description", "Class", "Class Description"]
    card = pd.read_table(card, sep='\t', header=0)
    card2 = gene2.merge(card, on='Gene ID', how='left')
    card2 = card2.fillna("-")
    num = float(len(card2.columns) - 2) / len(card2.columns) * 100
    card2 = card2[card2.apply(lambda x: cent_fun(x) < num, axis=1)]
    card2.index = card2["Gene ID"]
    card2.to_csv(out + ".card_anno.xls", sep='\t', header=True, index=False)
    card = card.ix[:, ["Gene ID", "ARO_Accession", "ARO_description", "ARO_category"]]
    card.columns = ["Gene ID", "ARO Accession", "ARO Description", "ARO Category"]
    phi = pd.read_table(phi, sep='\t', header=0)
    phi2 = gene2.merge(phi, on='Gene ID', how='left')
    phi2 = phi2.fillna("-")
    num = float(len(phi2.columns) - 2) / len(phi2.columns) * 100
    phi2 = phi2[phi2.apply(lambda x: cent_fun(x) < num, axis=1)]
    phi2.index = phi2["Gene ID"]
    phi2.to_csv(out + ".phi_anno.xls", sep='\t', header=True, index=False)
    phi = phi.ix[:, ["Gene ID", "PHI ID", "Gene Function", "Pathogen Species", "Host Species"]]
    phi.columns = ["Gene ID", "PHI ID", "Gene Function", "Pathogen Species", "Host Species"]
    vfdb = pd.read_table(vfdb, sep='\t', header=0)
    vfdb["Gene ID"] = vfdb["#Query"]
    del vfdb["#Query"]
    vfdb2 = gene2.merge(vfdb, on='Gene ID', how='left')
    vfdb2 = vfdb2.fillna("-")
    num = float(len(vfdb2.columns) - 2) / len(vfdb2.columns) * 100
    vfdb2 = vfdb2[vfdb2.apply(lambda x: cent_fun(x) <= num, axis=1)]
    vfdb2.index = vfdb2["Gene ID"]
    vfdb2.to_csv(out + ".vfdb_anno.xls", sep='\t', header=True, index=False)
    vfdb = vfdb.ix[:, ["Gene ID", "VFDB ID", "VFs", "VFs_Description"]]
    vfdb.columns = ["Gene ID", "VFDB ID", "VFs", "VFs Description"]
    signalp = pd.read_table(signalp1, sep='\t', header=0)
    del signalp["Location"]
    signalp5 = gene2.merge(signalp, on='Gene ID', how='left')
    signalp5 = signalp5.fillna("-")
    num = float(len(signalp5.columns) - 2) / len(signalp5.columns) * 100
    signalp5 = signalp5[signalp5.apply(lambda x: cent_fun(x) < num, axis=1)]
    signalp5.index = signalp5["Gene ID"]
    signalp5.to_csv(out + "_Gram+_SignalP.txt", sep='\t', header=True, index=False)
    signalp['secretory protein(G+)']="Yes"
    signalp = signalp.ix[:, ["Gene ID", "secretory protein(G+)"]]
    signalp.columns = ["Gene ID", "secretory protein(G+)"]
    signalp4 = pd.read_table(signalp2, sep='\t', header=0)
    del signalp4["Location"]
    signalp3 = gene2.merge(signalp4, on='Gene ID', how='left')
    signalp3 = signalp3.fillna("-")
    num = float(len(signalp3.columns) - 2) / len(signalp3.columns) * 100
    signalp3 = signalp3[signalp3.apply(lambda x: cent_fun(x) < num, axis=1)]
    signalp3.index = signalp3["Gene ID"]
    signalp3.to_csv(out + "_Gram-_SignalP.txt", sep='\t', header=True, index=False)
    signalp4['secretory protein(G-)'] = "Yes"
    signalp4 = signalp4.ix[:, ["Gene ID", "secretory protein(G-)"]]
    signalp4.columns = ["Gene ID", "secretory protein(G-)"]
    tmhmm = pd.read_table(tmhmm, sep='\t', header=0)
    tmhmm2 = gene2.merge(tmhmm, on='Gene ID', how='left')
    tmhmm2 = tmhmm2.fillna("-")
    num = float(len(tmhmm2.columns) - 2) / len(tmhmm2.columns) * 100
    tmhmm2 = tmhmm2[tmhmm2.apply(lambda x: cent_fun(x) < num, axis=1)]
    tmhmm2.index = tmhmm2["Gene ID"]
    tmhmm2.to_csv(out + ".tmhmm_anno.xls", sep='\t', header=True, index=False)
    tmhmm = tmhmm.ix[:, ["Gene ID", "Number of predicted TMHs"]]
    tmhmm.columns = ["Gene ID", "TMH No"]
    tcdb = pd.read_table(tcdb, sep='\t', header=0)
    tcdb2 = gene2.merge(tcdb, on='Gene ID', how='left')
    tcdb2 = tcdb2.fillna("-")
    num = float(len(tcdb2.columns) - 2) / len(tcdb2.columns) * 100
    tcdb2 = tcdb2[tcdb2.apply(lambda x: cent_fun(x) < num, axis=1)]
    tcdb2.index = tcdb2["Gene ID"]
    tcdb2.to_csv(out + ".tcdb_anno.xls", sep='\t', header=True, index=False)
    tcdb = tcdb.ix[:, ["Gene ID", "TCDB Description"]]
    tcdb.columns = ["Gene ID", "TCDB Description"]
    secretory = pd.read_table(secretory, sep='\t', header=0)
    secretory["Gene ID"] = secretory["#Query"]
    del secretory["#Query"]
    secretory2 = gene2.merge(secretory, on='Gene ID', how='left')
    secretory2 = secretory2.fillna("-")
    num = float(len(secretory2.columns) - 2) / len(secretory2.columns) * 100
    secretory2 = secretory2[secretory2.apply(lambda x: cent_fun(x) < num, axis=1)]
    secretory2.index = secretory2["Gene ID"]
    secretory2.to_csv(out + ".secretory_anno.xls", sep='\t', header=True, index=False)
    secretory = secretory.ix[:, ["Gene ID", "Type"]]
    secretory.columns = ["Gene ID", "Secretion Type"]
    a = gene.merge(cog, on='Gene ID', how='left')
    a = a.merge(kegg, on='Gene ID', how='left')
    a = a.merge(cazy, on='Gene ID', how='left')
    a = a.merge(card, on='Gene ID', how='left')
    a = a.merge(phi, on='Gene ID', how='left')
    a = a.merge(vfdb, on='Gene ID', how='left')
    a = a.merge(signalp, on='Gene ID', how='left')
    a = a.merge(signalp4, on='Gene ID', how='left')
    a = a.merge(tmhmm, on='Gene ID', how='left')
    a = a.merge(tcdb, on='Gene ID', how='left')
    a = a.merge(secretory, on='Gene ID', how='left')
    a = a.fillna("-")
    a = a.set_index('Gene ID', drop=False)
    a = a.ix[:,["Gene ID", "Location", "Strand", "Start", "End", "Gene Length(bp)", "Type","COG ID", "COG Description", "Function", "Fun Description", "Gene Name", "KO ID", "KO Description", "Family", "Family Description", "Class", "Class Description", "ARO Accession", "ARO Description", "ARO Category", "PHI ID", "Gene Function", "Pathogen Species", "Host Species", "VFDB ID", "VFs", "VFs Description", "secretory protein(G+)", "secretory protein(G-)", "TMH No", "TCDB Description", "Secretion Type"]]
    a["sample"] = sample
    a.to_csv(out + ".anno_summary.xls", sep='\t', header=True, index=False)

def _main():
    parser = argparse.ArgumentParser(description='add function datbase in your table ')
    parser.add_argument('-n', help="sample name")
    parser.add_argument('-gff', help="file")
    parser.add_argument('-cog', help="cog")
    parser.add_argument('-kegg', help="kegg")
    parser.add_argument('-cazy',  help="cazy")
    parser.add_argument('-card',  help="card")
    parser.add_argument('-phi', help="phi")
    parser.add_argument('-vfdb',  help="vfdb")
    parser.add_argument('-signalp', help="signalp(G+)")
    parser.add_argument('-signalps',help="signalp(G-)")
    parser.add_argument('-tmhmm', help="tmhmm")
    parser.add_argument('-tcdb', help="tcdb")
    parser.add_argument('-secretory', help="secretory")
    parser.add_argument('-o', help="output file")
    args = parser.parse_args()
    add_gene_info(args.n, args.gff, args.cog, args.kegg, args.cazy, args.card, args.phi, args.vfdb, args.signalp, args.signalps, args.tmhmm, args.tcdb, args.secretory, args.o)


if __name__ == "__main__":
    _main()
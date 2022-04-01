#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from collections import defaultdict 
import pickle 


## 读入公共数据库：KEGG， HMDB， LIPIDMAPS， EPA，SMPDB； 
##返回几个字典： dic_name2cid, dic_name2cpd, dic_name2cas, dic_name2formula
def get_stardard(name):
    dict_al2al = {
    u'α':u'alpha',
    u'β':u'beta',
    u'γ':u'gamma',
    u'δ':u'delta',
    u'ε':u'epsilon',
    u'ζ':u'zeta',
    u'η':u'eta',
    u'θ':u'theta',
    u'ι':u'iota',
    u'κ':u'kappa',
    u'λ':u'lambda',
    u'μ':u'mu',
    u'ν':u'nu',
    u'ξ':u'xi',
    u'ο':u'omicron',
    u'π':u'pi',
    u'ρ':u'rho',
    u'σ':u'sigma',
    u'τ':u'tau',
    u'υ':u'upsilon',
    u'φ':u'phi',
    u'χ':u'chi',
    u'ψ':u'psi',
    u'ω':u'omega',
    u'&plusmn;': u'+/-',
    u'&':u'',
    u'′':u"'",
    u'<WBR>':u'',
    u'±': u'+/-',
    u'Δ': u'delta',
    u'- ': u'-',
    }
    try:
        for alpha in dict_al2al:
            if alpha in name:
                name = name.replace(alpha, dict_al2al[alpha])
    except:
        pass
    if ' / ' in name:
        name = name.split(' / ')[0]
    name_len = len(name)
    if name_len > 0:
        name_head = name[0]
        name_tail = name[1:name_len + 1]
        name = ''.join([name_head.upper(), name_tail])
    return name
    
#处理每个字典：对键 值去重；
def short_dict(mydict):
    ##需要为列表式字典
    new_dict= defaultdict(list)
    for key ,values in mydict.items():
        new_values = list(set(values))
        for val in new_values:
            new_dict[key].append(val) 
    return new_dict

#把多个字典合并，
def combine_many_dict(mydicts):
    combine_dict = defaultdict(list)
    for my_dt in mydicts:
        new_my_dt = short_dict(my_dt)
        for key,value in new_my_dt.items():
            combine_dict[key].extend(value)
    new_combine_dict = short_dict(combine_dict) 
    return new_combine_dict



def get_public_anno_info(out_file,kegg_file,hmdb_file,lmdb_file,smpdb_file):
    ##1.开始处理kegg_file；
    #kegg_df = pd.read_csv(kegg_file)
    # kegg_file = 'Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\KEGG_compound.20210526_addHMDB.xls'
    # hmdb_file = 'Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\hmdb_metabolites.detail.xls'
    # lmdb_file = 'Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\LIPIDMAPS.detail.xls'
    # smpdb_file = 'Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\smpdb_metabolites.20180914.xls'
    name_cid_dict_kegg = defaultdict(list)
    name_formula_dict_kegg = defaultdict(list) 
    name_cas_dict_kegg = defaultdict(list)  
    kegg_df = pd.read_csv(kegg_file,sep="\t",quotechar="'") 
    # print(kegg_df.columns.tolist())
    for row in kegg_df.itertuples(): 
        cid = getattr(row, 'ENTRY')
        names = getattr(row, 'NAME')
        cpd = getattr(row, "HMDBID")
        cas = getattr(row, "KEGGCAS")
        formula = getattr(row,"FORMULA")
        for name in names.split(";"):
            name = name.lower()
            name_cid_dict_kegg[name].append(cid)
            name_cid_dict_kegg[get_stardard(name)].append(cid)
            name_formula_dict_kegg[name].append(formula)
            name_formula_dict_kegg[get_stardard(name)].append(formula)
            if cas != "_":
                name_cas_dict_kegg[name].append(cas)
                name_cas_dict_kegg[get_stardard(name)].append(cas)
    #2.处理HMDBID；
    hmdb_df = pd.read_csv(hmdb_file,sep="\t",quotechar='"') 
    # print(hmdb_df.columns.tolist())
    name_cid_dict_hmdb = defaultdict(list)
    name_cpd_dict_hmdb = defaultdict(list)
    name_formula_dict_hmdb = defaultdict(list) 
    name_cas_dict_hmdb = defaultdict(list) 
    
    for row in hmdb_df.itertuples():
        cpd = getattr(row,"accession")
        names =';'.join([getattr(row,"name"),getattr(row,"iupac_name"),getattr(row,"traditional_iupac"),getattr(row, "synonyms")]) 
        cid = getattr(row,"kegg_id")
        formula = getattr(row, "chemical_formula")
        cas = getattr(row,"cas_registry_number")
        for name in names.split(";"):
            if name != "_":
                name = name.lower()
                name_cpd_dict_hmdb[name].append(cpd)
                name_cpd_dict_hmdb[get_stardard(name)].append(cpd)
                name_formula_dict_hmdb[name].append(formula)
                name_formula_dict_hmdb[get_stardard(name)].append(formula)
                if cas != "_":
                    name_cas_dict_hmdb[name].append(cas)
                    name_cas_dict_hmdb[get_stardard(name)].append(cas)
                if cid != '_':
                    name_cid_dict_hmdb[name].append(cid)
                    name_cid_dict_hmdb[get_stardard(name)].append(cid)

    ##3.开始处理LIPIDMAPS
    lpdb_df= pd.read_csv(lmdb_file,sep="\t",quotechar="'") 
    # print(lpdb_df.columns.tolist())
    
    name_cid_dict_lpdb = defaultdict(list)
    name_cpd_dict_lpdb = defaultdict(list)
    name_formula_dict_lpdb = defaultdict(list) 
    name_cas_dict_lpdb = defaultdict(list) 
    
    for row in lpdb_df.itertuples():
        cpd = getattr(row,"LM_ID")
        names =';'.join([getattr(row,"NAME"),getattr(row,"SYSTEMATIC_NAME"),getattr(row,"SYNONYMS")]) 
        cid = getattr(row,"KEGG_ID")
        formula = getattr(row, "FORMULA")
        for name in names.split(";"):
            if name != "_":
                name = name.lower()
                name_cpd_dict_lpdb[name].append(cpd)
                name_cpd_dict_lpdb[get_stardard(name)].append(cpd)
                name_formula_dict_lpdb[name].append(formula)
                name_formula_dict_lpdb[get_stardard(name)].append(formula)
                if cid != '_':
                    name_cid_dict_lpdb[name].append(cid)
                    name_cid_dict_lpdb[get_stardard(name)].append(cid)

    ##4.开始处理SMPDB的数据库：
    smpdb_df= pd.read_csv(smpdb_file, sep="\t",quotechar="'") 
    # print(smpdb_df.columns.tolist())
    
    name_cid_dict_smpdb = defaultdict(list)
    name_cpd_dict_smpdb = defaultdict(list)
    name_formula_dict_smpdb = defaultdict(list) 
    name_cas_dict_smpdb = defaultdict(list) 
    
    for row in smpdb_df.itertuples():
        #print(row)
        cpd = getattr(row,"MetaboliteID")
        #names =';'.join([getattr(row,"Metabolite Name"),getattr(row,"SYSTEMATIC_NAME"),getattr(row,"SYNONYMS")])
        name= getattr(row,"MetaboliteName").lower()
        cid = getattr(row,"KEGGID")
        formula = getattr(row, "Formula")
        cas = getattr(row,"CAS")
        name_cpd_dict_smpdb[name].append(cpd)
        name_cpd_dict_smpdb[get_stardard(name)].append(cpd)
        name_formula_dict_smpdb[name].append(formula)
        name_formula_dict_smpdb[get_stardard(name)].append(formula)
        if cid != '_':
            name_cid_dict_smpdb[name].append(cid)
            name_cid_dict_smpdb[get_stardard(name)].append(cid)
        if cas != "_":
            name_cas_dict_smpdb[name].append(cas)
            name_cas_dict_smpdb[get_stardard(name)].append(cas) 
                    
##处理全部数据库；
    dic_name2cas = combine_many_dict([name_cas_dict_kegg,name_cas_dict_hmdb,name_cas_dict_smpdb])
    dic_name2cpd = combine_many_dict([name_cpd_dict_hmdb,name_cpd_dict_lpdb,name_cas_dict_smpdb])
    dic_name2cid  = combine_many_dict([name_cid_dict_kegg,name_cid_dict_hmdb,name_cid_dict_lpdb,name_cid_dict_smpdb])
    dic_name2formula = combine_many_dict([name_formula_dict_kegg,name_formula_dict_hmdb,name_formula_dict_lpdb,name_formula_dict_smpdb])
    with open(out_file, 'wb') as ow:
            pickle.dump([dic_name2cid, dic_name2cpd, dic_name2cas, dic_name2formula], ow)
    return dic_name2cid, dic_name2cpd, dic_name2cas, dic_name2formula 


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description="处理publicMetabolite db to get kegg id")
#     parser.add_argument("-kegg",default='Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\KEGG_compound.20210526_addHMDB.xls',help="input file as KEGG_compound.20210526_addHMDB.xls")
#     parser.add_argument("-hmdb",default='Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\hmdb_metabolites.detail.xls',help="input file as class.hmdb_metabolites.xls")
#     parser.add_argument("-lipids",default='Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\LIPIDMAPS.detail.xls',help="input file as LIPIDMAPS.detail.xls")
#     parser.add_argument("-smpdb",default='Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\smpdb_metabolites.20180914.xls',help="input file as SMPDB_metabolite.xls")
#     args = parser.parse_args()
#     dic_name2cid, dic_name2cpd, dic_name2cas, dic_name2formula = get_public_anno_info(args.kegg, args.hmdb, args.lipids, args.smpdb)
    
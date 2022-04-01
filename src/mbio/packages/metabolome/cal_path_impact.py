# -*- coding: UTF-8 -*-
import os
import re
import pandas as pd
from collections import defaultdict
import json
from scipy.stats import fisher_exact, hypergeom
import sys
import argparse

parser = argparse.ArgumentParser(description="Pathway Topology Analysis")
#parser.add_argument("-p_class", type=str, required=True, help="ko.pathway.xls")
parser.add_argument("-cid", type=str, required=True, help="name2cid.xls")
parser.add_argument("-org", type=str, required=True, help="hsa. is subfolder of db_json")
parser.add_argument("-method", type=str, default='rbc', help="rbc or rod")
parser.add_argument("-out", type=str, required=True, help="outfile for topology.xls")
parser.add_argument("-db_json", type=str, required=True, help="database dir path.include imp.json,node.json and *pathway.xls")
parser.add_argument("-backgroup", type=str, help='backgroup file. must contain ID column')

args = parser.parse_args()

primary_list = ['01100', '01110', '01120', '01130', '01200', '01210', '01212', '01230', '01220']

if args.backgroup:
    backgroup_file = args.backgroup
    backgroup_data = pd.read_table(backgroup_file,sep='\t',header=0)
    back_map_ids =[i[-5:] for i in backgroup_data['ID'].tolist()]
    need_back = True
else:
    need_back = False

CIDs = []
cid_map_metab_id = {}  #add 20200414
with open(args.cid, "r") as cpd_r:
    for line in cpd_r.readlines():
        items = line.strip().split("\t")
        CPD = items[0]
        CID = items[1]
        if ';' in CID.strip():
            CID_ = CID.split(';')
            CIDs.extend(CID_)
            for c_id in CID_:
                if c_id not in cid_map_metab_id:
                    cid_map_metab_id[c_id] = []
                cid_map_metab_id[c_id].append(CPD)
        else:
            CIDs.append(CID)
            if CID not in cid_map_metab_id:
                cid_map_metab_id[CID] = []
            cid_map_metab_id[CID].append(CPD)

query_num = len(set(CIDs))

#db_dir = '/mnt/ilustre/users/ting.kuang/ITRAQ/db/KEGG/'
#db_dir = '/mnt/ilustre/users/sanger-dev/app/database/metabolome/topo_json/'
db_dir = args.db_json
org = args.org
#org_ = org.split('.')[1]
org_ = org

method = args.method
with open(db_dir + org + "/imp.json", 'r') as load_imp, open(db_dir + org + "/nodes.json", 'r') as load_node, \
        open(db_dir + org + "/%s.pathway.xls" %org_, 'r') as path_cpd_r, open('%s_tmp' % args.out, 'w') as out_w:
    load_imp = json.load(load_imp)
    load_node = json.load(load_node)
    allcid_in_path = sum([x for x in load_node.values()], [])
    allcid_in_path_num = len(set(allcid_in_path))
    allcid_in_path_new = [x.lstrip('cpd:') for x in allcid_in_path]
    query_in_path = [x for x in set(CIDs) if x in allcid_in_path_new]
    query_in_path_num = len(query_in_path)
    All_CIDs = []
    ko_2ko_name = {}
    for line in path_cpd_r.readlines():
        items = line.strip().split("\t")
        if len(items) < 6:
            continue
        koid = items[2].replace(':', '')
        ko_name = items[3]
        ko_2ko_name[koid] = ko_name
    out_w.write('ko\tpathway_name\tmatch_status\tp_value\timpact_value\tcompound_ids\tmetab_ids\n')
    for ko in load_node.keys():
        if ko[-5:] in primary_list:
            continue
        if need_back:
            if ko[-5:] not in back_map_ids:
                continue

        cds_in_path = load_node[ko]
        if len(cds_in_path) == 0:
            continue
        cds_in_path_new = [x.lstrip('cpd:') for x in cds_in_path]
        All_CIDs.append(cds_in_path_new)
        cds_in_path_num = len(set(cds_in_path_new))
        hit_cds = [cid for cid in CIDs if cid in cds_in_path_new]
        if hit_cds:
            hit_cds_set = set(hit_cds)
            hit_num = len(hit_cds_set)
            a = hit_num
            c = query_in_path_num - hit_num
            b = cds_in_path_num - hit_num
            d = allcid_in_path_num - cds_in_path_num - query_in_path_num + hit_num
            p_value = hypergeom.pmf(a, a + b + c + d, a + b, c)
            fish_p_value = fisher_exact([[a, b], [c, d]])[1]
            imp_info = load_imp[ko]
            impact_value = 0
            for x in imp_info:
                c_ = re.sub(r'cpd:', '', x[0])
                if method == "rbc":
                    if c_ in set(hit_cds):
                        impact_value += x[1]
                if method == "rod":
                    if c_ in set(hit_cds):
                        impact_value += x[2]

            metab_ids = []
            for c_id in hit_cds_set:
                metab_ids.append(';'.join(set(cid_map_metab_id[c_id])))

            out_w.write('%s\t%s\t%s|%s\t%s\t%s\t%s\t%s\n' % (
                    ko, ko_2ko_name[ko], hit_num, cds_in_path_num, p_value, impact_value,'|'.join(hit_cds_set), '|'.join(metab_ids)))

data = pd.read_csv('%s_tmp' % args.out, index_col=0,sep="\t")
data = data.sort_values(by="p_value")
uncorrected_p = data['p_value']
len_p = len(uncorrected_p)
p_rank = uncorrected_p.rank()
p_info = pd.concat([uncorrected_p, p_rank], axis=1, ignore_index=False)
p_info.columns = ['uncorrected_p', 'p_rank']


def cal_fdr(x, y):
    return len_p * x / (y + 1)


p_info['FDR'] = p_info.apply(lambda row: cal_fdr(row[0], row[1]), axis=1)

p_info = p_info.drop(['uncorrected_p', 'p_rank'], axis=1)
new_result = pd.concat([data, p_info], axis=1, ignore_index=False)
col_name = data.columns.tolist()
col_name.insert(col_name.index('p_value') + 1, 'FDR')
new_result = new_result.reindex(columns=col_name)
new_result.to_csv(args.out, sep="\t", header=True, index=True)

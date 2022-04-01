# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190326

"""
对注释结果进行筛选和统计
region_select: {"chr1": "0-10000"}
"""

import os
import json
import argparse


def region_filter(item, region_info):
    """
    区域筛选
    """
    if item[4] not in region_info.keys():
        return False
    for start_end in region_info[item[4]]:
        start_end = start_end.split("-")
        if start_end[0] and start_end[1]:
            if int(start_end[0]) <= int(item[5]) and int(item[6]) <= int(start_end[1]):
                return True
        elif start_end[0]:
            if int(item[5]) >= int(start_end[0]):
                return True
        elif start_end[1]:
            if int(item[6]) <= int(start_end[1]):
                return True
        else:
            return True
    return False

def region_anno_stat(region_select, pop_summary, output_dir):
    if region_select != "all":
        region_select = json.loads(region_select)
    region_summary = os.path.join(output_dir, "region.summary")
    kegg_stat_path = os.path.join(output_dir, "region.kegg.stat")
    go_stat_path = os.path.join(output_dir, "region.go.stat")
    eggnog_stat_path = os.path.join(output_dir, "region.eggnog.stat")
    pfam_stat_path = os.path.join(output_dir, "region.pfam.stat")
    kegg_detail, go_detail, egg_detail, pfam_detail = {}, {}, {}, {}
    kegg_stat, go_stat, egg_stat, pfam_stat = {}, {}, {}, {}
    all_stat = {"kegg": {"total": 0, "enrich": 0}, "go": {"total": 0, "enrich": 0}, "egg": {"total": 0, "enrich": 0}, "pfam": {"total": 0, "enrich": 0}}
    ko_list, go_list, egg_list, pfam_list = [], [], [], []
    with open(pop_summary, "rb") as f, open(region_summary, "wb") as w:
        for line in f:
            if line.startswith("#"):
                w.write(line)
                continue
            item = line.strip().split("\t")
            if region_select != "all":
                status = region_filter(item, region_select)
                if not status:
                    continue
            w.write(line)
            tran_id = item[2].split("|")[-2]
            eff = int(item[7]) + int(item[8])
            if item[15] != "--":
                all_stat["kegg"]["total"] += 1
                if eff > 0:
                    all_stat["kegg"]["enrich"] += 1
                ko_ids = item[15].split(",")
                kegg_annos = item[16].split(":")
                for i in range(len(ko_ids)):
                    if ko_ids[i] == "--":
                        continue
                    if ko_ids[i] not in ko_list:
                        ko_list.append(ko_ids[i])
                        kegg_stat[ko_ids[i]] = {"eff": 0, "total": 0, "gene": []}
                    kegg_stat[ko_ids[i]]["total"] += 1
                    kegg_stat[ko_ids[i]]["gene"].append(tran_id)
                    if eff > 0:
                        kegg_stat[ko_ids[i]]["eff"] += 1
                    kegg_detail[ko_ids[i]] = kegg_annos[i]
            if item[17] != "--":
                all_stat["go"]["total"] += 1
                if eff > 0:
                    all_stat["go"]["enrich"] += 1
                go_ids = item[17].split(",")
                go_annos = item[18].split(":")
                for i in range(len(go_ids)):
                    if go_ids[i] == "--":
                        continue
                    if go_ids[i] not in go_list:
                        go_list.append(go_ids[i])
                        go_stat[go_ids[i]] = {"eff": 0, "total": 0, "gene": []}
                    go_stat[go_ids[i]]["total"] += 1
                    go_stat[go_ids[i]]["gene"].append(tran_id)
                    if eff > 0:
                        go_stat[go_ids[i]]["eff"] += 1
                    go_detail[go_ids[i]] = go_annos[i]
            if item[19] != "--":
                all_stat["egg"]["total"] += 1
                if eff > 0:
                    all_stat["egg"]["enrich"] += 1
                egg_ids = item[19].split(":")[0].split(",")
                egg_annos = item[19].split(":")[1].split(";")
                for i in range(len(egg_ids)):
                    if egg_ids[i] == "--":
                        continue
                    if egg_ids[i] not in egg_list:
                        egg_list.append(egg_ids[i])
                        egg_stat[egg_ids[i]] = {"eff": 0, "total": 0, "gene": []}
                    egg_stat[egg_ids[i]]["total"] += 1
                    egg_stat[egg_ids[i]]["gene"].append(tran_id)
                    if eff > 0:
                        egg_stat[egg_ids[i]]["eff"] += 1
                    egg_detail[egg_ids[i]] = egg_annos[i]
            if item[21] != "--":
                all_stat["pfam"]["total"] += 1
                if eff > 0:
                    all_stat["pfam"]["enrich"] += 1
                pfam_ids = item[21].split(",")
                pfam_annos = item[22].split(":")
                for i in range(len(pfam_ids)):
                    if pfam_ids[i] == "--":
                        continue
                    if pfam_ids[i] not in pfam_list:
                        pfam_list.append(pfam_ids[i])
                        pfam_stat[pfam_ids[i]] = {"eff": 0, "total": 0, "gene": []}
                    pfam_stat[pfam_ids[i]]["total"] += 1
                    pfam_stat[pfam_ids[i]]["gene"].append(tran_id)
                    if eff > 0:
                        pfam_stat[pfam_ids[i]]["eff"] += 1
                    pfam_detail[pfam_ids[i]] = pfam_annos[i]
    with open(go_stat_path, "wb") as w:
        w.write("#GO ID\tDescription\tEffection\tTotal\tTotal Gene\n")
        for go in go_list:
            w.write(go + "\t" + go_detail[go] + "\t" + str(go_stat[go]["eff"]) + "\t" + str(go_stat[go]["total"]) + "\t")
            w.write(str(all_stat["go"]["enrich"]) + "\t" + str(all_stat["go"]["total"]) + "\t" + ";".join(go_stat[go]["gene"]) + "\n")
    with open(kegg_stat_path, "wb") as w:
        w.write("#ko ID\tDescription\tEffection\tTotal\tTotal Gene\n")
        for ko in ko_list:
            w.write(ko + "\t" + kegg_detail[ko] + "\t" + str(kegg_stat[ko]["eff"]) + "\t" + str(kegg_stat[ko]["total"]) + "\t")
            w.write(str(all_stat["kegg"]["enrich"]) + "\t" + str(all_stat["kegg"]["total"]) + "\t" + ";".join(kegg_stat[ko]["gene"]) + "\n")
    with open(eggnog_stat_path, "wb") as w:
        w.write("#Functional Categories\tEffection\tTotal\tTotal Gene\n")
        for egg in egg_list:
            w.write(egg + "\t" + egg_detail[egg] + "\t" + str(egg_stat[egg]["eff"]) + "\t" + str(egg_stat[egg]["total"]) + "\t")
            w.write(str(all_stat["egg"]["enrich"]) + "\t" + str(all_stat["egg"]["total"]) + ";".join(egg_stat[egg]["gene"]) + "\n")
    with open(pfam_stat_path, "wb") as w:
        w.write("#Pfam Accession\tPfam Annotation\tEffection\tTotal\tTotal Gene\n")
        for pfam in pfam_list:
            w.write(pfam + "\t" + pfam_detail[pfam] + "\t" + str(pfam_stat[pfam]["eff"]) + "\t" + str(pfam_stat[pfam]["total"]) + "\t")
            w.write(str(all_stat["pfam"]["enrich"]) + "\t" + str(all_stat["pfam"]["total"]) + ";".join(pfam_stat[pfam]["gene"]) + "\n")


parser = argparse.ArgumentParser(description="基因功能注释区域筛选及统计")
parser.add_argument("-i", "--pop_summary", required=True)
parser.add_argument("-r", "--region_select", required=True)
parser.add_argument("-o", "--output_dir", required=True)

args = vars(parser.parse_args())

region_anno_stat(args["region_select"], args["pop_summary"], args["output_dir"])

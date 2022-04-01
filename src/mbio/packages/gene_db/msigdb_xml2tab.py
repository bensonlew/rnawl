# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys
import xml.etree.ElementTree as ET

global class2name
class2name ={
    "H": "hallmark gene sets",
    "C1": "positional gene sets" ,
    "C2": "curated gene sets",
    "C3": "motif gene sets" ,
    "C4": "computational gene sets" ,
    "C5": "GO gene sets" ,
    "C6": "oncogenic signatures" ,
    "C7": "immunologic signatures" ,
    "CGP": "chemical and genetic perturbations",
    "BIOCARTA": "BioCarta gene sets" ,
    "REACTOME": "Reactome gene sets" ,
    "KEGG": "KEGG gene sets" ,
    "CP": "Canonical pathways" ,
    "PID": "PID gene sets" ,
    "MIR": "microRNA targets" ,
    "MIR_Legacy": "Legacy microRNA targets",
    "MIRDB": "MIRDB microRNA targets",
    "TFT": "transcription factor targets",
    "TFT_Legacy": "Legacy transcription factor target",
    "GTRD": "GTRD transcription factor targets",
    "CM": "cancer modules" ,
    "CGN": "cancer gene neighborhoods" ,
    "CC": "GO cellular component" ,
    "BP": "GO biological process" ,
    "MF": "GO molecular function"
}

def get_gmt(gmt_path):
    g_dic = dict()
    with open(gmt_path) as f:
        for line in f:
            cols = line.strip().split("\t")
            g_dic[cols[0]] = {
                'url': cols[1],
                'symbols': " ".join(cols[2:])
            }
    return g_dic

def get_xml(xml_path, gmt_path):
    g_dic = get_gmt(gmt_path)
    global class2name
    msig = ET.parse(xml_path)
    root = msig.getroot()
    a = root.getchildren()
    # b = a[0].getchildren()
    print "name\tc1\tc2\tc3\tbrief_description\turl\torganism\tgene_symbols"
    for ele in a:
        ele_dic = ele.attrib
        name = ele_dic.get("STANDARD_NAME", "")
        c1 = ele_dic.get("CATEGORY_CODE", "")
        if not c1 in class2name:
            continue
        c1 = "{}: {}".format(c1, class2name[c1])
        c2 = ele_dic.get("SUB_CATEGORY_CODE", "").split(":")[0]
        try:
            c3 = ele_dic.get("SUB_CATEGORY_CODE", "").split(":")[1]
        except:
            c3 = ""
        if c3 == "":
            if c2 != "":
                c2 = "{}: {}".format(c2, class2name[c2])
        else:
            if c3 != "":
                c3 = "{}: {}".format(c3, class2name[c3])

        brief_description = ele_dic.get("DESCRIPTION_BRIEF", "")
        url = g_dic[name]['url']
        organism = ele_dic.get("ORGANISM", "")
        gene_symbols = g_dic[name]['symbols']



        print "\t".join([name,c1,c2,c3,brief_description,url,organism, gene_symbols])

if __name__ == "__main__":
    get_xml(sys.argv[1], sys.argv[2])

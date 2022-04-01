# -*- coding: utf-8 -*-
"""
author: liubinxu
lastmodified:20180919
"""

import os
import re
import sys
import pygraphviz as pgv
from pygraphviz import *
from colour import Color
from numpy import linspace
from math import log10
import matplotlib
matplotlib.use('agg')
import pylab as pl
from Bio import Phylo
from cStringIO import StringIO
from mbio.packages.rna.accession2taxon import taxon


class Taxon(object):
    def __init__(self):
        self.children = []
        self.gene_num = 0
        self.name = ""

taxons_dict  = dict()
gene_num = 0
ncbi_taxon = taxon()


taxon_root_obj =  Taxon()
taxon_root_obj.gene_num = 0
taxon_root_obj.name = "root"
taxons_dict["root"] = taxon_root_obj

# 获取所有分类
taxon_file = sys.argv[1]
with open(taxon_file, 'r') as taxon_f:
    for line in taxon_f.readlines():
        gene_num += 1
        taxon_info = line.strip().split("\t")[2]
        taxons = taxon_info.split(";")

        parents = "root"
        for taxon in taxons:
            taxon_name = taxon.split("{")[0]
            taxon_id = taxon_name.replace(" ", "_").replace("'", "_")
            now_taxon = parents + "_" + taxon_id
            if now_taxon in taxons_dict:
                taxons_dict[now_taxon].gene_num += 1
            else:
                taxon_obj =  Taxon()
                taxon_obj.gene_num = 1
                taxon_obj.name = taxon_name
                taxon_obj.parent = parents
                taxons_dict[parents].children.append(now_taxon)
                taxons_dict[now_taxon] = taxon_obj
            parents = now_taxon

taxons_dict["root"].gene_num = gene_num
# 过滤分类中gene数量低的分类
leaf_num = 0
def filter_children():
    global taxons_dict
    for taxon in taxons_dict.keys():
        if taxon == "root":
            pass
        else:
            if taxons_dict[taxon].gene_num < gene_num * 0.0002:
                parent = taxons_dict[taxon].parent
                taxons_dict[parent].children.remove(taxon)
                if parent + 'other' in taxons_dict[parent].children:
                    taxons_dict[parent + 'other'].gene_num += taxons_dict[taxon].gene_num
                else:
                    taxon_obj =  Taxon()
                    taxon_obj.gene_num = taxons_dict[taxon].gene_num
                    taxon_obj.name = "other"
                    taxons_dict[parent].children.append(parent + 'other')
                    taxons_dict[parent + 'other'] = taxon_obj
            else:
                pass

'''
def filter_children(taxon_name):
    other = 0
    global leaf_num
    global taxons_dict
    if len(taxons_dict[taxon_name].children) == 0:
        leaf_num += 1
    else:
        pass
    childrens = taxons_dict[taxon_name].children
    for taxon in childrens:
        if taxons_dict[taxon].gene_num < gene_num * 0.001:
            print taxon + ":" + str(taxons_dict[taxon].gene_num)
            other += taxons_dict[taxon].gene_num
            #print taxons_dict[taxon_name].children
            taxons_dict[taxon_name].children.remove(taxon)
            #print taxons_dict[taxon_name].children
            leaf_num += 1
        else:
            print "tree: " + taxon + ":" + str(taxons_dict[taxon].gene_num)
            filter_children(taxon)
    if other > 0:
        taxon_obj =  Taxon()
        taxon_obj.gene_num = other
        taxon_obj.name = "other"
        taxons_dict[taxon_name].children.append(taxon_name + 'other')
        taxons_dict[taxon_name + 'other'] = taxon_obj
'''

filter_children()

def get_tree_string(taxon_name):
    global leaf_num
    name_list = []
    for taxon_child in taxons_dict[taxon_name].children:
        if len(taxons_dict[taxon_child].children) > 0:
            name_list.append(taxon_child + get_tree_string(taxon_child))
        else:
            leaf_num += 1
            name_list.append(taxon_child)
    name = "(" + ", ".join(name_list) + ")"
    return name

def get_show(taxon_name):
    if str(taxon_name) == "":
        return ""
    elif not str(taxon_name) in taxons_dict:
        print "taxon_name {} not exists".format(str(taxon_name))
    else:
        taxon_name = str(taxon_name)
        pct =  float(taxons_dict[taxon_name].gene_num)/float(gene_num) * 100
        return taxons_dict[taxon_name].name + "," + \
            ncbi_taxon.get_taxon_from_name(taxons_dict[taxon_name].name) + \
            "\nnum:" + str(taxons_dict[taxon_name].gene_num) + "(" + '%.3f' % pct  + "%)"

tree_string = get_tree_string("root")
print tree_string
tree = Phylo.read(StringIO(tree_string),  "newick")

height = leaf_num/10
pl.figure(figsize=(30,height))
pl.rcParams['figure.figsize'] = [15,leaf_num/3]
pl.rcParams['font.size'] = 6

try:
    Phylo.draw(tree, label_func=lambda x:get_show(x.name), axes=None)
except:
    print "物种太少，无法建树"

pl.axis('off')

pl.savefig(sys.argv[1] + ".pdf")

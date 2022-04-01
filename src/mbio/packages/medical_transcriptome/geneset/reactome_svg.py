# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys
import xml.etree.ElementTree as ET
from multiprocessing import Pool

class ReactomeSvg(object):
    def __init__(self, colors=['#FF0000', '#0000FF'], complex_child=dict()):
        self.colors = colors
        self.complex_child = complex_child

    def rec_get_child(self, vertex_id, child_list=[]):
        '''
        递归查找所有子节点
        '''
        # self.child_list = list()
        if vertex_id in self.complex_child:
            for child_id in self.complex_child[vertex_id]:
                #通路互相调用
                if child_id in child_list:
                    pass
                else:
                    child_list.append(child_id)
                child_list = self.rec_get_child(child_id, child_list)
            return child_list
        else:
            return child_list


    def svg_mark(self, svg_path, geneset2gene=None, rgene2gene=None, path_name="", output_dir=None):
        # flog = open(path_name + ".log", 'w')
        # f_log.write("{}\n".format(locals()))
        # flog.write("{}\n".format(self.complex_child))

        spe_abr ="-".join(path_name.split("-")[:2])

        ET.register_namespace("","http://www.w3.org/2000/svg")
        svg = ET.parse(svg_path)
        root = svg.getroot()

        layer_a = root.findall("{http://www.w3.org/2000/svg}g")
        b = layer_a[0]
        c = b.findall("{http://www.w3.org/2000/svg}g")
        # d = c[-1]
        mark_file = open(output_dir + '/{}.mark'.format(path_name), 'w')
        marked_set = set()
        for d in c:
            for e in d:
                tag =  e.tag.split("}")[1]
                attr = e.attrib

                # f_log.write("{}\n".format(attr))

                if "sgbn_id" in attr:
                    sgbn_id = e.get("sgbn_id", "")

                    # 跳过冗余


                    if sgbn_id.split("_")[0] == "arc":
                        continue
                    sgbn_vertex_id = sgbn_id.split("_")[1]
                    rgene_id = spe_abr + "-" + sgbn_id.split("_")[1]
                    # f_log.write("rgene is {}\n".format(rgene_id))

                    genes = list()
                    rgenes = list()
                    if rgene_id in rgene2gene:
                        rgenes = [rgene_id]
                        genes = rgene2gene[rgene_id]
                    elif sgbn_vertex_id in self.complex_child:
                        gene_list = []
                        child_list = []
                        child_list = set(self.rec_get_child(sgbn_vertex_id, child_list))
                        for child_gene in child_list:
                            if spe_abr + "-" + child_gene in rgene2gene:
                                rgenes.append(spe_abr + "-" + child_gene)
                            gene_list.extend(rgene2gene.get(spe_abr + "-" + child_gene, []))
                        genes = list(set(gene_list))
                        # flog.write("sgbn_vertex_id is  {}\n".format(sgbn_vertex_id))
                        # flog.write("rgene is  {}\n".format(rgenes))
                        # flog.write("child  {}\n".format(child_list))
                    if len(genes) == 0:
                        continue

                    genesets = geneset2gene.keys()
                    if len(genesets) == 1:
                        genes_set1 = set(genes).intersection(set(geneset2gene[genesets[0]]))
                        if len(genes_set1) > 0:
                            e.set("stroke-width", "3")
                            e.set("stroke", self.colors[0])
                            genes_list = list(genes_set1)
                            geneset_list = [genesets[0] for gene in genes_set1]
                            title_list = ["{}({}); ".format(gene, genesets[0]) for gene in genes_set1]
                            if sgbn_id in marked_set:
                                pass
                            else:
                                mark_file.write("\t".join([
                                    path_name,
                                    ";".join(rgenes),
                                    sgbn_id,
                                    self.colors[0],
                                    ";".join(genes_list),
                                    ";".join(geneset_list),
                                    "\\n".join(title_list)
                                ]) + "\n")
                                marked_set.add(sgbn_id)
                    else:
                        genes_set1 = set(genes).intersection(set(geneset2gene[genesets[0]]))
                        genes_set2 = set(genes).intersection(set(geneset2gene[genesets[1]]))
                        color = ""
                        if len(genes_set1) > 0 and len(genes_set2) > 0:
                            e.set("stroke-width", "3")
                            e.set("stroke", "#FFFF00")
                            color = "#FFFF00"
                        elif len(genes_set1) > 0:
                            e.set("stroke-width", "3")
                            e.set("stroke", self.colors[0])
                            color = self.colors[0]
                        elif len(genes_set2) > 0:
                            e.set("stroke-width", "3")
                            e.set("stroke", self.colors[1])
                            color = self.colors[1]
                        if len(genes_set1) > 0 or len(genes_set2) > 0:
                            genes_12 = genes_set1.intersection(genes_set2)
                            genes_1 = genes_set1 - genes_set2
                            genes_2 = genes_set2 - genes_set1
                            genes_list = list(genes_12) + list(genes_1) + list(genes_2)
                            geneset_list = [genesets[0] + "," + genesets[1] for gene in genes_12] + \
                                           [genesets[0] for gene in genes_1] + \
                                           [genesets[1] for gene in genes_2]
                            title_list = ["{}({}); ".format(gene, genesets[0] + "," + genesets[1]) for gene in genes_12] + \
                                         ["{}({}); ".format(gene, genesets[0]) for gene in genes_1] + \
                                         ["{}({}); ".format(gene, genesets[1]) for gene in genes_2]
                            if sgbn_id in marked_set:
                                continue
                            else:
                                mark_file.write("\t".join([
                                    path_name,
                                    ";".join(rgenes),
                                    sgbn_id,
                                    color,
                                    ";".join(genes_list),
                                    ";".join(geneset_list),
                                    "\\n".join(title_list)
                                ])+ "\n")
                                marked_set.add(sgbn_id)
                else:
                    pass

        mark_file.close
        # flog.close
        svg.write(output_dir + "/" + path_name + ".svg")


def svg_mark(para_list):
    [path_name, colors, complex_child, path_reac_gene, svg_path1, geneset2gene, output_dir] = para_list
    # reac_svg = ReactomeSvg(colors = colors, complex_child=complex_child)
    # print "***"
    # reac_svg.svg_mark(svg_path1, geneset2gene=geneset2gene, rgene2gene=path_reac_gene, path_name=path_name, output_dir=output_dir)
    # print "###"
    try:
        reac_svg = ReactomeSvg(colors = colors, complex_child=complex_child)
        reac_svg.svg_mark(svg_path1, geneset2gene=geneset2gene, rgene2gene=path_reac_gene, path_name=path_name, output_dir=output_dir)
    except Exception as e:
        with open(output_dir + '/{}.log'.format(path_name), 'aw') as f:
            f.write("{}".format(e))
            # f.write("{}".format(para_list))

def mutiple_svg_convert(paras=None, thread=1):
    p = Pool(thread)
    p.map(svg_mark, paras)
    p.close()
    p.join()
    '''
    for para in paras:
        svg_mark(para)
    '''


if __name__ == "__main__":
    reac = ReactomeSvg()

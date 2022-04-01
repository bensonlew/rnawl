import pandas as pd
import matplotlib.pyplot as plt
import os
from collections import defaultdict
import json
import re
import pandas as pd
import lxml
from lxml import etree
import matplotlib
matplotlib.use("Agg")
def set_default(obj):
     if isinstance(obj, set):
        return list(obj)
     raise TypeError

def extract_bind_id_to_genes(mapping_file):
    mapping_file_name = os.path.splitext(os.path.basename(mapping_file))[0]
    b2g = defaultdict(set)
    g2b = defaultdict(set)
    a=pd.read_table(mapping_file)
    for i, row in a.iterrows():
        bind_code = row["BINCODE"].split("'")[1]
        gene_id = row["IDENTIFIER"].split("'")[1]
        g2b[gene_id].add(bind_code)
        b2g[bind_code].add(gene_id)
    return mapping_file_name,g2b,b2g

def extract_bind_id_from_xml(xml_file):
    bind_ids=[]
    id2pos=defaultdict(dict)
    pathway_name = os.path.splitext(os.path.basename(xml_file))[0]
    parser = etree.XMLParser(encoding="utf-8")
    doc = etree.parse(xml_file,parser=parser)
    a = doc.xpath('DataArea')
    for i in a:
        x_coord = i.get("x")
        y_coord = i.get("y")
        pos = (x_coord,y_coord)
        child = i.getchildren()[0]
        id = child.get("id")
        id2pos[id] = pos
        bind_ids.append(id)
    return pathway_name,bind_ids,id2pos

def change_svg_scale(svg_path,x_coord,y_coord,insert_width,insert_height):
    target_text=""
    target_text += '<svg x="{x_coord}" y="{y_coord}" width="{insert_width}" height="{insert_height}" fill-opacity="0.4" >\n'.format(
        x_coord=str(x_coord),y_coord=str(y_coord),insert_width=str(insert_width),insert_height=str(insert_height))
    with open(svg_path,"r") as svg:
        text_list = svg.readlines()
        for n,line in enumerate(text_list):
            if line.startswith("<svg"):
                starts_line = n
                break
        svg_text = "".join(text_list[starts_line:])
    # svg_text = open(svg_path,"r").read()
    svg_tex1 = re.sub(r'height="(.*?)"', 'height="{}"'.format(insert_height), svg_text, count=1)
    svg_texf = re.sub(r'width="(.*?)"','width="{}"'.format(insert_width),svg_tex1, count=1)
    target_text += svg_texf
    target_text += '\n</svg>\n'
    return target_text

def insert_chart_svg_to_basic(basic_svg="",draw_svg_infos={},target_file_path=""):
    with open(basic_svg,"r") as raw:
        basic_text = raw.read()
        head = "</svg>".join(basic_text.split("</svg>")[:-1])
        tail = "</svg>"+basic_text.split("</svg>")[-1]
    for draw_svg in draw_svg_infos:
        draw_svg_path = draw_svg["file_path"]
        x_coord = draw_svg["x_coord"]
        y_coord = draw_svg["y_coord"]
        insert_width =50
        insert_height = 50
        insert_text = change_svg_scale(draw_svg_path,x_coord,y_coord,insert_width,insert_height)
        head += insert_text
    head+= tail
    with open(target_file_path,"w") as w:
        w.write(head)



def draw_line_svgs(work_dir="",data_df =None,svg_name = ""):
    # samples = data_df.columns[1:]
    # x_data = list(samples)
    # for index, row in data_df.iterrows():
    #     y_data = list(row[1:])
    #     plt.plot(x_data, y_data)
    # plt.savefig(os.path.join(work_dir,svg_name), format="svg")

    samples = data_df.columns[:]
    x_data = list(samples)
    for index, row in data_df.iterrows():
        y_data = list(row[:])
        plt.plot(x_data, y_data)
    plt.savefig(os.path.join(work_dir, svg_name), format="svg")














#-*- coding: utf-8 -*-
from collections import defaultdict
from lxml import etree
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from reportlab.lib import colors
from copy import deepcopy
import os
import re


home = os.path.abspath(__file__).split("wpm2")[0]
db_path = os.path.join(home, "app/database/Annotation/all/KEGG/v2021.09.18/html")
image_magick = os.path.join(home, "app/program/ImageMagick/bin/convert")


def mark_kegg(anno2color, pathway, outdir, db_path=db_path):
    """
    anno2color: { cmp/k id : color}
    """
    print(">>>")
    print(pathway)
    print(anno2color)
    print("<<<")
    m_kegg = MarkKEGG(anno2color, pathway, outdir, db_path)
    m_kegg.mark_to()
    return m_kegg.pdf_map


class MarkKEGG(object):
    def __init__(self, anno2color, pathway, outdir, db_path):
        self.marks = anno2color
        self.db_path = db_path
        self.pathway = re.sub(r'^\D*', 'map', pathway)
        self.outdir = outdir
        self.html = None
        self.pdf_map = ''

    def save_pdf(self):
        pathway = self.pathway
        outdir = self.outdir
        map_kgml = os.path.join(self.db_path, pathway + ".kgml")
        map_png = pathway + ".png"
        if os.path.exists(map_png):
            os.remove(map_png)

        os.link(os.path.join(self.db_path, map_png), map_png)
        kgml = KGML_parser.read(open(map_kgml))
        kgml.image = map_png
        for entry in kgml.entries.values():
            all_names = entry.link.split("/")[-1].split('?')[-1].split('+')
            hits = set(all_names) & set(self.marks.keys())
            for g in entry.graphics:
                g.bgcolor = colors.Color(alpha=0)
                if hits:
                    c = self.marks[list(hits)[0]]
                    g.bgcolor = colors.HexColor(c + '33', hasAlpha=True)
                # coords = ''
                # if g.type in ["rectangle", "roundrectangle"]:
                #     x1 = int((2 * float(g.x) - float(g.width)) / 2)
                #     x2 = int((2 * float(g.x) + float(g.width)) / 2)
                #     y1 = int((2 * float(g.y) - float(g.height)) / 2)
                #     y2 = int((2 * float(g.y) + float(g.height)) / 2)
                #     coords = "{},{},{},{}".format(x1, y1, x2, y2)
                # elif g.type == "circle":
                #     coords = "{},{},{}".format(g.x, g.y, int(g.width) / 2)
                #

        canvas = KGMLCanvas(kgml, import_imagemap=True, label_compounds=False, label_orthologs=False,
                            label_reaction_entries=False, label_maps=True)
        pdf_out = os.path.join(outdir, pathway + '.pdf')
        png_out = os.path.join(outdir, pathway + '.png')
        canvas.draw(pdf_out)
        os.system("{} -flatten -quality 100 -density 130 -background white {} {}".format(image_magick, pdf_out, png_out))
        self.pdf_map = pathway

    def save_mark(self):
        pathway = self.pathway
        map_html = os.path.join(self.db_path, pathway + ".html")
        mark_out = open(os.path.join(self.outdir, pathway + ".html.mark"), 'w')
        mark_out.write('ko\tquery\tshape\tcoords\ttitle\thref\tcolor\n')
        with open(map_html, 'r') as r:
            html = etree.HTML(r.read())
        for area in html.xpath("//area"):
            if "shape" not in area.attrib or area.attrib["shape"] not in ["rect", "circle"]:
                continue
            all_names = area.attrib["href"].split("/")[-1].split('?')[-1].split('+')
            hits = set(all_names) & set(self.marks)
            if not hits:
                continue
            mark_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(pathway, ','.join(hits), area.attrib["shape"],
                                                                 area.attrib["coords"], area.attrib["title"],
                                                                 area.attrib["href"], self.marks[list(hits)[0]]))

    def mark_to(self):
        for post in [".html", ".kgml", ".png"]:
            f_path = os.path.join(self.db_path, self.pathway + post)
            if not os.path.exists(f_path):
                print("{} 没有图片".format(self.pathway))
                return
        try:
            self.save_pdf()
            self.save_mark()
        except:
            print("{} 没有图片".format(self.pathway))
            pass

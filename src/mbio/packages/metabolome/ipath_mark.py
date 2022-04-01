#-*- coding: utf-8 -*-
from collections import defaultdict
from lxml import etree
from copy import deepcopy
import os
import re


home = os.path.abspath(__file__).split("wpm2")[0]
db_path = os.path.join(home, "app/database/IPATH3_metab")
db_type = ["Antibiotics", "Metabolism", "Microbial_metabolism", "Secondary_metabolites"]


def id_to_anno_type(annot, db_path=db_path):
    """
    annot: {annot: set(ids)}
    """
    rets = defaultdict(set)
    for t in db_type:
        with open(os.path.join(db_path, t + ".mapping"), 'r') as r:
            for l in r:
                anno = l.split('\t')[0]
                if anno in annot:
                    rets[t].update(annot[anno])
    return rets


def mark_ipath(marks, outdir, db_path=db_path):
    """
    marks: {color: [name, c or k set]}
           键为颜色，值为 两个元素的list，第一个为名称，第二个为set，所有的注释id( K id, compound id)
    """
    print(marks)
    m_ipath = MarkIpath(marks, db_path)
    for t in db_type:
        outfile = os.path.join(outdir, t + '.svg')
        m_ipath.mark_to(t, outfile)


class MarkIpath(object):
    def __init__(self, mark_dict, db_path):
        self.ori_mark = mark_dict
        self.marks = deepcopy(self.ori_mark)
        self.db_path = db_path
        self.svg = None
        self.anchor = None

    def get_svg(self, ipath_type):
        self.marks = deepcopy(self.ori_mark)
        ipath_svg = os.path.join(self.db_path, ipath_type + ".raw.svg")
        ipath_maping = os.path.join(self.db_path, ipath_type + ".mapping")
        ipath = {}
        with open(ipath_maping, 'r') as r:
            for l in r:
                l = l.strip().split('\t')
                if l[0] not in ipath:
                    ipath[l[0]] = []
                ipath[l[0]].extend(l[1].split(','))
        for color in self.marks:
            self.marks[color][1] = set(sum([ipath[i] for i in self.marks[color][1] if i in ipath], []))

        print(self.marks)
        with open(ipath_svg, 'r') as r:
            self.svg = etree.HTML(r.read())
        self.anchor = self.svg.xpath("//*[@id='legend']")[0]

    def add_legend(self):
        leg_mark, leg_text = self.anchor.xpath('*')
        #self.anchor.clear()
        mark_start = int(re.findall(r'M(\d+)', leg_mark.attrib['d'])[0])
        for color in self.marks:
            name = self.marks[color][0]

            leg_m = deepcopy(leg_mark)
            leg_m.attrib["stroke"] = color
            leg_m.attrib['d'] = "M{},75L{},75".format(mark_start, mark_start + 75)

            leg_t = deepcopy(leg_text)
            leg_t.text = name
            leg_t.attrib['x'] = str(mark_start + 85)
            leg_t.attrib["fill"] = "#000000"

            mark_start += 100 + 36 * len(name)
            self.anchor.extend([leg_m, leg_t])

        leg_mark.clear()
        leg_text.clear()

    def ele_mark(self):
        for color in self.marks:
            eles = self.marks[color][1]
            for ele in eles:
                print(ele)
                print(self.svg)
                ele = self.svg.xpath("//*[@id='{}']".format(ele))
                if not ele:
                    continue
                new_ele = deepcopy(ele[0])
                if new_ele.tag == "ellipse":
                    new_ele.attrib.update({
                        "fill": color,
                        "rx": str(int(new_ele.attrib['rx']) * 3),
                        "ry": str(int(new_ele.attrib['ry']) * 3),
                        "style": "opacity: 0.6; cursor: pointer;"
                    })
                else:
                    style = new_ele.attrib["style"]
                    style = dict([re.split(r':\s?', i) for i in re.split(r';\s?', style) if i])
                    style['opacity'] = '0.6'
                    style['stroke-width'] = str(int(style["stroke-width"]) * 5)
                    new_ele.attrib.update({
                        "stroke": color,
                        "style": ";".join(["{}: {}".format(k, v) for k, v in style.items()])
                    })
                new_ele.attrib["marked"] = "1"
                ele[0].addnext(new_ele)

    def mark_to(self, path_type, outfile):
        print(path_type)
        self.get_svg(path_type)
        self.add_legend()
        self.ele_mark()
        self.svg.getroottree().write(outfile)
        with open(outfile, 'r') as r:
            content = r.read()
        with open(outfile, 'w') as w:
            content = re.sub(r"(?ms).*(\<svg[\S\s\n]+svg\>).*$", r"\1", content)
            w.write(content)


# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys
import xml.etree.ElementTree as ET

class Reactome(object):
    def __init__(self):
        self.sbgn_ele = dict()
        self.guess_shapes = list()
        self.eles = list()

    def svg_convert(self, svg_path, test=False):
        ET.register_namespace("","http://www.w3.org/2000/svg")
        svg = ET.parse(svg_path)
        root = svg.getroot()

        root.set("width", root.attrib['viewBox'].split(" ")[2])
        root.set("height", root.attrib['viewBox'].split(" ")[3])


        svg.write(svg_path.replace(".svg", ".changed.svg"))

    def svg_convert2(self, svg_path, test=False):
        ET.register_namespace("","http://www.w3.org/2000/svg")
        svg = ET.parse(svg_path)
        root = svg.getroot()

        layer_a = root.findall("{http://www.w3.org/2000/svg}g")
        b = layer_a[0]
        c = b.findall("{http://www.w3.org/2000/svg}g")
        # d = c[-1]
        for d in c:
            if "fill" in d.attrib:
                if "fill-opacity" in d.attrib:
                    d.attrib["fill"] = "rgba({},{})".format(d.attrib["fill"][4:-1], d.attrib["fill-opacity"])
                    del d.attrib["fill-opacity"]
                else:
                    d.attrib["fill"] = "rgba({},{})".format(d.attrib["fill"][4:-1], "1")

            if "stroke" in d.attrib:
                if "stroke-opacity" in d.attrib:
                    d.attrib["stroke"] = "rgba({},{})".format(d.attrib["stroke"][4:-1], d.attrib["stroke-opacity"])
                    del d.attrib["stroke-opacity"]
                else:
                    d.attrib["stroke"] = "rgba({},{})".format(d.attrib["stroke"][4:-1], "1")

        svg.write(svg_path.replace(".svg", ".changed.svg"))


if __name__ == "__main__":
    reac = Reactome()

    for p in sys.argv[1:]:
        print p
        reac = Reactome()
        reac.svg_convert2(p)

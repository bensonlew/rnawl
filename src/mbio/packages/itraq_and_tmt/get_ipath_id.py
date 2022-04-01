# -*- coding: utf-8 -*-
# __author__ = '刘彬旭'
from biocluster.iofile import File
import os
import xml.etree.ElementTree as ET
import re
import sys


def get_ipath_id(svg_files):
    """
    根据ipath结果提取id 属性对应关系
    """
    xml = ET.parse(svg_files)
    root = xml.getroot()
    layer_a = root.findall("{http://www.w3.org/2000/svg}g")
    for one_of_a in layer_a:
        layer_b = one_of_a.findall("{http://www.w3.org/2000/svg}g")
        for one_of_b in layer_b:
            layer_c = one_of_b.findall("{http://www.w3.org/2000/svg}g")
            for one_of_c in layer_c:
                m_id = "no_id"
                if one_of_c.attrib.has_key('id'):
                    m_id = one_of_c.attrib['id']
                else:
                    pass
                layer_d = one_of_c.findall("{http://www.w3.org/2000/svg}path")
                layer_d2 = one_of_c.findall("{http://www.w3.org/2000/svg}ellipse")
                for one_of_d in layer_d:
                    print m_id,
                    attr_list = ['style', 'stroke-width' , 'stroke','fill', 'd']
                    for each_attr in attr_list:
                        if one_of_d.attrib.has_key(each_attr):
                            print "\t" + one_of_d.attrib[each_attr],
                        else:
                            print "\t" + one_of_d.attrib[each_attr],
                    print "\n",
                for one_of_d in layer_d2:
                    print m_id,
                    attr_list = ['style', 'stroke-width' , 'stroke', 'fill', 'cx', 'cy']
                    for each_attr in attr_list:
                        if one_of_d.attrib.has_key(each_attr):
                            print "\t" + one_of_d.attrib[each_attr],
                        else:
                            print "\t" + one_of_d.attrib[each_attr],
                    print "\n",
if __name__ == '__main__':  # for test
    get_ipath_id(sys.argv[1])

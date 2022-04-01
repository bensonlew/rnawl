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

    def bbox2shape(self, shape_type=None, bbox=None, center=None):
        if shape_type == "rect":
            shape_dict = {
                "x": float(bbox['x']),
                "y": float(bbox['y']),
                "width": float(bbox['w']),
                "height": float(bbox['h']),
            }
        if shape_type == "rect_r":
            shape_dict = {
                "x": float(bbox['x']) + 4.0,
                "y": float(bbox['y']) + 4.0,
                "width": float(bbox['w']) - 8.0,
                "height": float(bbox['h']) - 8.0,
            }
        elif shape_type == "ellipse":
            shape_dict = {
                "cx": float(bbox['x']) + float(bbox['w'])/2,
                "cy": float(bbox['y']) + float(bbox['h'])/2,

                "rx": float(bbox['w'])/2,
                "ry": float(bbox['h'])/2
            }
        elif shape_type == "ellipse_path":
            shape_dict = {
                "x": float(bbox['x']),
                "y": float(bbox['y']),
                "color": "rgb(165,215,145)"
            }
        elif shape_type == "complex":
            shape_dict = {
                "ele2": float(bbox['y']),
                "ele3": float(bbox['x']),
                "color": "rgb(171,209,227)"
            }
        elif shape_type == "circle":
            shape_dict = {
                "cx": float(center[0]),
                "cy": float(center[1])
            }
        elif shape_type == "port":
            shape_dict = {
                "x1": float(bbox['x']),
                "y1": float(bbox['y']),
                "x2": float(center[0]),
                "y2": float(center[1])
            }
        else:
            pass

        return shape_dict

    def line_start_arrow(self, start):
        shape = {
            "ele3": float(start['x']),
            "ele4": float(start['y'])
        }
        return shape

    def line_end_circle(self, start):
        shape = {
            "cx": float(start['x']),
            "cy": float(start['y']),
            "range": 4,
        }
        return shape

    def port2line(self, port_type="port", start=None, end=None):
        shape = {
            "x1": float(start[0]),
            "x2": float(end[0]),
            "y1": float(start[0]),
            "y2": float(end[0])
        }
        return shape


    def point_to_line(self, point_list):
        shapes = list()
        for i in range(0, len(point_list) -1):
            for j in range(i+1, len(point_list)):
                shape = {
                    "x1": float(point_list[i][0]),
                    "x2": float(point_list[j][0]),
                    "y1": float(point_list[i][1]),
                    "y2": float(point_list[j][1])
                }
                shapes.append(shape)
        return shapes


    def get_sbgn(self, sgbn_path):
        sgbn = ET.parse(sgbn_path)
        root =  sgbn.getroot()
        a = root.getchildren()
        b = a[0].getchildren()
        sgbn_eles = dict()
        shapes = list()
        for ele in b:
            ele_label = ""
            ele_bbox = dict()
            ele_port = list()
            ele_type = ele.tag.split("}")[1]
            ele_dic = ele.attrib

            if ele_type == "glyph":
                bbox_center = (None, None)
                for child in ele.getchildren():
                    child_type = child.tag.split("}")[1]
                    child_dict = child.attrib
                    if child_type == 'label':
                        ele_label = child_dict.get("text", "")
                    if child_type == 'bbox':
                        ele_bbox = child_dict
                        if ele_dic['class'] == "macromolecule":
                            shape = ('rect', ele_dic['id'], self.bbox2shape('rect', child_dict), ele_dic['id'])
                            self.guess_shapes.append(shape)
                        if ele_dic['class'] == "submap":
                            shape = ('rect', ele_dic['id'], self.bbox2shape('rect', child_dict), ele_dic['id'])
                            self.guess_shapes.append(shape)
                        if ele_dic['class'] == "simple chemical":
                            shape = ('ellipse_path',  ele_dic['id'], self.bbox2shape('ellipse_path', child_dict), ele_dic['id'])
                            self.guess_shapes.append(shape)
                            shape = ('ellipse',  ele_dic['id'], self.bbox2shape('ellipse', child_dict), ele_dic['id'])
                            self.guess_shapes.append(shape)
                        if ele_dic['class'] in ["process", "uncertain process"]:
                            shape = ('rect',  ele_dic['id'], self.bbox2shape('rect', child_dict), ele_dic['id'])
                            self.guess_shapes.append(shape)
                            bbox_center = (float(child_dict['x']) + float(child_dict['w'])/2, float(child_dict['y']) + float(child_dict['h'])/2)
                        if ele_dic['class'] in ["dissociation", "association"]:
                            bbox_center = (float(child_dict['x']) + float(child_dict['w'])/2, float(child_dict['y']) + float(child_dict['h'])/2)
                            shape = ('circle',  ele_dic['id'], self.bbox2shape('circle', center=bbox_center), ele_dic['id'])
                            self.guess_shapes.append(shape)
                        if ele_dic['class'] == "complex":
                            shape = ('complex', ele_dic['id'],  self.bbox2shape('complex', child_dict), ele_dic['id'])
                            self.guess_shapes.append(shape)
                        if ele_dic['class'] == "unspecified entity":
                            shape = ('rect',  ele_dic['id'], self.bbox2shape('rect', child_dict), ele_dic['id'])
                            self.guess_shapes.append(shape)
                            shape = ('rect',  ele_dic['id'], self.bbox2shape('rect_r', child_dict), ele_dic['id'])
                            self.guess_shapes.append(shape)
                    if child_type == 'port':
                        ele_port.append(child_dict)
                        if ele_dic['class'] in ["process", "uncertain process", "dissociation", "association"]:
                            shape = ('line',
                                     ele_dic['id'],
                                     self.bbox2shape('port', child_dict, bbox_center),
                                     ele_dic['id'])
                            self.guess_shapes.append(shape)
                self.eles.append((ele_type, ele_dic['id'], ele_label))

            elif ele_type == "arc":
                point_list = list()

                if ele_dic['class'] in ["catalysis", "inhibition", "consumption", "stimulation"]:
                    show_id = ele_dic['target'].split(".")[0]
                elif ele_dic['class'] in ["production"]:
                    show_id = ele_dic['source'].split(".")[0]

                for child in ele.getchildren():
                    child_type = child.tag.split("}")[1]
                    child_dict = child.attrib
                    if child_type == 'start':
                        point_list.append((child_dict['x'], child_dict['y']))
                        shape = ("arrow", ele_dic['id'], self.line_start_arrow(child_dict), show_id)
                        self.guess_shapes.append(shape)

                    if child_type == 'next':
                        point_list.append((child_dict['x'], child_dict['y']))
                    if child_type == 'end':
                        point_list.append((child_dict['x'], child_dict['y']))
                        shape = ("circle", ele_dic['id'], self.line_end_circle(child_dict), show_id)
                        self.guess_shapes.append(shape)
                        shape = ("arrow", ele_dic['id'], self.line_start_arrow(child_dict), show_id)
                        self.guess_shapes.append(shape)
                shapes = self.point_to_line(point_list)
                for line in shapes:
                    shape = ('line', ele_dic['id'], line, show_id)
                    self.guess_shapes.append(shape)

    def guess_id(self, tag, attr_dict):
        '''
        根据svg坐标猜测对应的id
        '''
        # print attr_dict
        if tag == "rect":
            for shape in self.guess_shapes:
                shape_dict = shape[2]
                # print shape_dict
                if shape[0] == "rect" and float(attr_dict['x']) == shape_dict['x'] and float(attr_dict['y']) == shape_dict['y'] and float(attr_dict['width']) == shape_dict['width'] and float(attr_dict['height']) == shape_dict['height']:
                    return shape, "match"
                else:
                    pass
            return False, "match"

        elif tag == "path":
            ds = attr_dict['d'].split(" ")
            for shape in self.guess_shapes:
                shape_dict = shape[2]
                if shape[0] == 'ellipse_path':
                    if ds[2].startswith("C"):
                        if float(ds[4]) == shape_dict['x'] and float(ds[3]) == shape_dict['y']:
                            return shape, "match"
                elif shape[0] == 'complex':
                    if len(ds) >= 3 and ds[2].startswith("L"):
                        if float(ds[2].lstrip("L")) == shape_dict['ele3'] and float(ds[1]) == shape_dict['ele2']:
                            return shape, "match"
                else:
                    pass
            return False, "match"
        elif tag == "polygon":
            points = attr_dict['points'].strip().split(" ")
            for shape in self.guess_shapes:
                shape_dict = shape[2]
                if shape[0] == 'arrow':
                    if float(points[2]) == shape_dict['ele3'] and float(points[3]) == shape_dict['ele4']:
                        return shape, "match"
                elif shape[0] == 'complex':
                    if len(points) == 18:
                        if float(points[-4]) == shape_dict['ele3'] and float(points[-1]) == shape_dict['ele2']:
                            return shape, "match"
            return False, "match"
        elif tag == "line":
            for shape in self.guess_shapes:
                shape_dict = shape[2]
                if shape[0] == 'line':
                    if float(attr_dict['x1']) == shape_dict['x1'] and float(attr_dict['x2']) == shape_dict['x2'] and float(attr_dict['y1']) == shape_dict['y1'] and float(attr_dict['y2']) == shape_dict['y2']:
                        return shape, "match"
                    elif float(attr_dict['x1']) == shape_dict['x2'] and float(attr_dict['x2']) == shape_dict['x1'] and float(attr_dict['y1']) == shape_dict['y2'] and float(attr_dict['y2']) == shape_dict['y1']:
                        return shape, "match"
                    elif float(attr_dict['x1']) == shape_dict['x1'] and float(attr_dict['y1']) == shape_dict['y1']:
                        print "line not in sgbn {}".format(attr_dict)
                        ## 添加线的箭头
                        shape_add_dict = {
                            "ele3": float(attr_dict['x2']),
                            "ele4": float(attr_dict['y2'])
                        }
                        shape_add = ("arrow", shape[1], shape_add_dict, shape[3])
                        self.guess_shapes.append(shape_add)
                        return shape, "guess"
                else:
                    pass
            return False, "match"
        elif tag == "circle":
            for shape in self.guess_shapes:
                shape_dict = shape[2]
                if shape[0] == 'circle':
                    if abs(shape_dict['cx'] - float(attr_dict['cx'])) <= 4.0 and abs(shape_dict['cy'] - float(attr_dict['cy'])) <= 4.0:
                        return shape, "match"
            return False, "match"
        elif tag == "ellipse":
            for shape in self.guess_shapes:
                shape_dict = shape[2]
                if shape[0] == 'ellipse':
                    if shape_dict['cx'] == float(attr_dict['cx']) and shape_dict['cy'] == float(attr_dict['cy']):
                        return shape, "match"
            return False, "match"

        else:
            print "tag {} not in sgbn".format(tag)
            return False, "match"

    def svg_convert(self, svg_path, test=False):
        ET.register_namespace("","http://www.w3.org/2000/svg")
        svg = ET.parse(svg_path)
        root = svg.getroot()

        layer_a = root.findall("{http://www.w3.org/2000/svg}g")
        b = layer_a[0]
        c = b.findall("{http://www.w3.org/2000/svg}g")
        # d = c[-1]
        for d in c:
            for e in d:
                tag =  e.tag.split("}")[1]
                attr = e.attrib
                # print attr
                # print self.guess_id(tag, attr)

                shape, guess_type = self.guess_id(tag, attr)
                if shape:
                    e.set('sgbn_id', shape[1])
                    e.set('show_id', shape[3])
                    if test:
                        if guess_type == "match":
                            e.set('stroke', "rgb(255,0,0)")
                        if guess_type == "guess":
                            e.set('stroke', "rgb(255,255,0)")
                        e.set('stroke-width', "5")
                else:
                    print "{} pos not found {}".format(tag, attr)

        svg.write(svg_path.replace(".svg", ".changed.svg"))


if __name__ == "__main__":
    reac = Reactome()
    reac.get_sbgn(sys.argv[1])

    with open(sys.argv[1] + ".shapes", 'w') as fo:
        for shape in reac.guess_shapes:
            fo.write("{}\n".format(shape))

    with open(sys.argv[1] + ".description", 'w') as fo:
        for ele in reac.eles:
            fo.write("{}\t{}\t{}\n".format(*ele))


    reac.svg_convert(sys.argv[2])

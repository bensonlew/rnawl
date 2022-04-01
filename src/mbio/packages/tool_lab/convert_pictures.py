#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2020/3/30 10:30
@file    : convert_pictures.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""

import fitz
from skimage import io
import cairosvg
import sys



def trim(img):
    img2 = img.sum(axis=2)
    img2_r = img2.sum(axis=1)
    img2_c = img2.sum(axis=0)
    (row,col) = img2.shape
    tempr0=0
    tempr1=0
    tempc0=0
    tempc1=0
    for r in range(0,row):
        if img2_r[r]!=765*col:
            tempr0=r
            break
    for r in range(row-1,0,-1):
        if img2_r[r]!=765*col:
            tempr1=r
            break
    for c in range(0,col):
        if img2_c[c]!=765*row:
            tempc0=c
            break
    for c in range(col-1,0,-1):
        if img2_c[c]!=765*row:
            tempc1=c
            break
    # print((tempr0, tempr1, tempc0, tempc1))
    new_img=img[tempr0:tempr1+3,tempc0:tempc1+3,0:3]
    return new_img


def convert(pic, *args, **kargs):
    if pic.endswith('.svg'):
        pdf = pic.split('.svg')[0] + '.pdf'
        cairosvg.svg2pdf(url=pic, write_to=pdf)
        pic = pdf
    pdfDoc = fitz.open(pic)
    zoom_x = 2
    zoom_y = 2
    mat = fitz.Matrix(zoom_x, zoom_y)
    pix = pdfDoc[0].getPixmap(matrix = mat)
    png = pic.split('.pdf')[0] + '.png'
    if 'out' in kargs:
        png = kargs['out']
    pix.writePNG(png)
    trim_ = True
    try:
        trim_ = kargs['trim']
    except:
        pass
    if trim_:
        im = io.imread(png)
        im = trim(im)
        io.imsave(png, im)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        _, pic = sys.argv
        png = pic.split('.pdf')[0] + '.png'
        convert(pic, trim=False)
        convert(pic, out=png + '_trim.png')
    elif len(sys.argv) == 3:
        _, pic, out = sys.argv
        convert(pic, out=out, trim=False)
        convert(pic, out=out + '_trim.png', trim=True)
    else:
        exit("USAGE:python $s picture out")



import fitz.fitz

import re
import os
from PIL import Image
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-pdf', required=True, help='PDF file')
parser.add_argument('-prex', default="page", help="the prex of output filename")
parser.add_argument('-type', default="png", help="output picture format")
parser.add_argument('-dpi', default="300", help="output file quality")
parser.add_argument('-weight', default=0, help="picture weight")
parser.add_argument('-output_dir', default="", help="output file directory")
args = parser.parse_args()
args.pdf = str(args.pdf)
args.prex = str(args.prex)
args.type = str(args.type)
args.dpi = int(args.dpi)
dirFile = str(args.output_dir)
weight = int(int(args.weight) * float(args.dpi / 25.4))

doc = fitz.open(args.pdf)  # open PDF file
i = 1
# dirFile = (args.pdf).split(".")[0]
# if not os.path.exists(dirFile):
#     os.mkdir(dirFile)
for page in doc:
    #    mat = fitz.Matrix(2, 2)  # zoom factor 2 in each direction
    pix = page.getPixmap()
    img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
    #    print(pix.width, pix.height)
    if args.weight != 0:
        # print(args.weight,pix.width,pix.height)
        rH = (float(weight) / pix.width) * pix.height
        img = img.resize((weight, int(rH)), Image.ANTIALIAS)

    picName = args.prex + "_" + str(i) + '.' + args.type
    img.save(dirFile + '/' + picName, quality=args.dpi)
    #    pix.writePNG(dirFile + '/' +picName+".png")
    i = i + 1

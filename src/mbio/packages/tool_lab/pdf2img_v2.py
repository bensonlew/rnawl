import fitz.fitz

import re
import os
import glob
from PIL import Image
import argparse

parser = argparse.ArgumentParser()
parser.add_argument( '-pdf' , required=True, help='PDF file')
parser.add_argument( '-type' , default= "png",help="output picture format")
parser.add_argument( '-dpi' ,  default = "300",help="output file quality")
#add by fwy
parser.add_argument('-output_dir', default="", help="output file directory")
# parser.add_argument( '-weight' , default =0,help="picture weight" )

args = parser.parse_args()
args.pdf = str(args.pdf)
args.type = str(args.type)
args.dpi = str(args.dpi)
dirFile = str(args.output_dir)

filename=os.path.splitext(os.path.basename(args.pdf))[0]
# filename = (args.pdf).split(".")[0]
# if not os.path.exists(dirFile):
#     os.mkdir(dirFile)

if (args.type == "png"):
    cmd = "pdftoppm -" + args.type + ' -r ' + args.dpi +' '+args.pdf + ' ' + dirFile + "/" + filename
elif (args.type == "jpg"):
    cmd = "pdftoppm -" + 'jpeg -r ' + args.dpi +' '+args.pdf + ' ' + dirFile + "/" + filename
elif (args.type == "tif"):
    cmd = "pdftoppm -" + 'tiff -r ' + args.dpi +' '+args.pdf + ' ' + dirFile + "/" + filename
else:
    print("only support png/jpg/tif\n")
    exit()

print(cmd)
os.system(cmd)



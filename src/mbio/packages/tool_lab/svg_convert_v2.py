import cairosvg
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '-svg' , required=True, help='SVG file')
parser.add_argument( '-type' , default= "png",help="output picture format")
parser.add_argument( '-dpi' ,  default = "300",help="output file quality")
parser.add_argument( '-out' ,help="output dir")

args = parser.parse_args()
args.dpi = str(args.dpi)
base=os.path.basename(args.svg)
svgprex = os.path.splitext(base)[0]
if (args.out is None):
    if not os.path.exists(svgprex):
        os.mkdir(svgprex)
    args.out = svgprex

# pdffile = args.out + "/" + svgprex+".pdf"
pdffile = svgprex+".pdf"

cairosvg.svg2pdf(url=args.svg, write_to=pdffile)
if args.type =="pdf":
    exit()
elif (args.type =="png"):
    cmd = "pdftoppm -" + 'png -r ' + args.dpi +' '+ pdffile + ' ' + args.out+"/"+ svgprex
elif (args.type == "jpg"):
    cmd = "pdftoppm -" + 'jpeg -r ' + args.dpi +' '+ pdffile + ' ' +args.out+"/"+ svgprex
elif (args.type == "tif"):
    cmd = "pdftoppm -" + 'tiff -r ' + args.dpi +' '+pdffile + ' ' +args.out+"/"+ svgprex
else:
    print("only support pdf/png/jpg/tif\n")
    exit()

print(cmd)
os.system(cmd)



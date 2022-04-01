import cairosvg
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '-svg' , required=True, help='SVG file')
parser.add_argument( '-type' , default= "png",help="output picture format")
parser.add_argument( '-dpi' ,  default = "300",help="output file quality")
parser.add_argument( '-weight' , default =0,help="picture weight" )
args = parser.parse_args()
args.dpi = int(args.dpi)
base=os.path.basename(args.svg)
svgprex = os.path.splitext(base)[0]

if args.type =="pdf":
    pdffile = svgprex+".pdf"
    cairosvg.svg2pdf(url=args.svg, write_to=pdffile)
if args.type is not "pdf":
    pngfile = svgprex + ".png"   
    cairosvg.svg2png(url=args.svg,dpi = args.dpi, write_to=pngfile)
    if args.weight != 0:
        os.system("convert -resize "+args.weight +" " +pngfile+" "+pngfile)
    if args.type =="png":
        exit    
    else:
        othertype = svgprex+"."+args.type
        os.system("convert "+pngfile+" "+othertype)




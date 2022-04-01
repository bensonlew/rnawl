import sys
import os
from PIL import Image
def convertPixel(r, g, b, a=1):
 color = "#%02X%02X%02X" % (r, g, b)
 opacity = a
 return (color, opacity)
for r in sys.argv[1:]:
 root, ext = os.path.splitext(r)
 image = Image.open(r)
 image = image.convert("RGB")
 mode = image.mode
 pixels = image.load()
 width, height = image.size
 print(image.mode)
 if "RGB" in mode:
  output = '<svg width="%d" height="%d" viewBox="0 0 %d %d" xmlns="http://www.w3.org/2000/svg">\n' % (width, height, width, height)
  for r in range(height):
   for c in range(width):
    color, opacity = convertPixel(*pixels[c, r])
    output += '<rect x="%d" y="%d" width="1" height="1" fill="%s" fill-opacity="%s"/>' % (c, r, color, opacity)
  output += "\n</svg>"
  with open(root + ".svg", "w") as f:
   f.write(output)
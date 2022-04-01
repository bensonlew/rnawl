import urllib
import sys
file_in = sys.argv[1]
with open(file_in, 'r') as f:
    for line in f:
        if '%' in line:
            print urllib.unquote(line.strip())
        else:
            print line

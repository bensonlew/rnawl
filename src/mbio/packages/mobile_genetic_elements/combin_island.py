# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re,os

def change_result(dir, name, output):
    island = {}
    files = os.listdir(dir)
    for file in files:
        if os.path.getsize(dir + "/" + file) >0:
            if re.search('dimob', file):
                with open(dir + "/" + file, "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        lin = line.strip().split("\t")
                        if lin[1]+"\t"+lin[2]+"\t"+lin[3] not in island.keys():
                            island[lin[1]+"\t"+lin[2]+"\t"+lin[3]] = 'Island_dimob'
                        else:
                            island[lin[1] + "\t" + lin[2] + "\t" + lin[3]] +=";" + 'Island_dimob'
            elif re.search('islander', file):
                with open(dir + "/" + file, "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        lin = line.strip().split("\t")
                        location = lin[0].split(" ")[0]
                        if lin[1]+"\t"+lin[2]+"\t"+lin[3] not in island.keys():
                            island[location+"\t"+lin[2]+"\t"+lin[3]] = 'Island_islander'
                        else:
                            island[location+"\t"+lin[2]+"\t"+lin[3]] +=";" + 'Island_islander'
            elif re.search('gihunter', file):
                with open(dir + "/" + file, "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        lin = line.strip().split("\t")
                        if lin[1]+"\t"+lin[2]+"\t"+lin[3] not in island.keys():
                            island[lin[0]+"\t"+lin[2]+"\t"+lin[3]] = 'Gihunter'
                        else:
                            island[lin[0] + "\t" + lin[2] + "\t" + lin[3]] +=";" + 'Gihunter'
    with open(output, "w") as g:
        g.write("Sample\tLocation\tStart\tEnd\tSofware\n")
        for key, value in island.items():
            g.write("{}\t{}\t{}\n".format(name, key, value))

def main():
    parser = OptionParser()
    parser.add_option('--n', dest='name', metavar='[sample name]')
    parser.add_option('--d', dest='dir', metavar='[files dir]')
    parser.add_option("--o", dest="output", metavar="[out file]")
    (options,args) = parser.parse_args()
    if not options.name or not options.dir or not options.output:
        print "python combin_island.py --n name --d dir --o output"
        return
    change_result(options.dir, options.name, options.output)

if __name__=='__main__':
    main()
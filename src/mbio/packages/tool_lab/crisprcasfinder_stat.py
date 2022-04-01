# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
import json
from pandas.core.frame import DataFrame
import pandas as pd
import shutil


def chang_result(data, sample, prefix):
    with open(data, "r") as f, open(prefix+".CRISPR_stat.xls", "w") as k, open(prefix+".CRISPR_information.xls","w") as d, open (prefix+".CRISPR_sequence.xls","w") as s:
        s.write("CRISPR ID\tType\tStart\tEnd\tSequence\n")
        data = json.load(f, encoding='utf-8')
        lista = []
        for aa in data['Sequences']:
            location = aa['Version']
            if len(aa['Crisprs']) >0:
                for bb in aa['Crisprs']:
                    list2 = [location]
                    len1 = int(bb['End']) - int(bb['Start'])
                    list2.append(len1)
                    dr_num = int(bb['Spacers'])+1
                    list2.append(dr_num)
                    cc = int(bb['DR_Length'])
                    spa_len = len1 - dr_num*cc
                    list2.append(spa_len)
                    for key, value in bb.items():
                       if key =="Name":
                          list2.append(value)
                       elif key == "End":
                          list2.append(value)
                       elif key == "Start":
                          list2.append(value)
                       elif key == "Spacers":
                          list2.append(value)
                       elif key == "DR_Length":
                          list2.append(value)
                       if key =="Regions":
                          for dd in value:
                             s.write(bb['Name']+"\t"+dd['Type']+"\t"+str(dd["Start"])+"\t"+str(dd["End"])+"\t"+dd["Sequence"]+"\n")
                    lista.append(list2)
        k.write("Sample\tCRISPR-Cas Num\n")
        k.write("{}\t{}\n".format(sample, str(len(lista))))
        d.write("CRISPR ID\tSample Name\tLocation\tStart\tEnd\tCRISPR Length（bp)\tDR Num\tDR Average Len (bp)\tSPA Num\tSPA Average Len (bp)\n")
        for i in lista:
            d.write(i[5]+"\t"+sample+"\t"+i[0]+"\t"+str(i[7])+"\t"+str(i[4])+"\t"+str(i[1])+"\t"+str(i[2])+"\t"+str(i[8])+"\t"+str(i[6])+"\t"+str(i[3])+"\n")

def main():
    parser = OptionParser()
    parser.add_option('--d', dest='data', metavar='[data json file]')
    parser.add_option('--s', dest='sample', metavar='[sample name]')
    parser.add_option("--o", dest="prefix", metavar="[out file prefix]")
    (options,args) = parser.parse_args()
    if not options.data or not options.sample or not options.prefix:
        print "python crisprcasfinder_stat.py --d data --s sample_name --o prefix "
        return
    chang_result(options.data,  options.sample, options.prefix)

if __name__=='__main__':
    main()

# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
import json
from pandas.core.frame import DataFrame
import pandas as pd
import shutil


def chang_result(data, sample, prefix):
    with open(data, "r") as f:
        data = json.load(f, encoding='utf-8')
        st = ''
        list1 = []
        for key, value in data['mlst']['results'].items():
            if key == "sequence_type":
                if value != '':
                    st = "ST-" + str(value)
                else:
                    st = "un"
            if key == "allele_profile":
                for key2, value2 in value.items():
                    list2 = []
                    list2.append(key2)
                    for key3 in value2.keys():
                        list2.append(value2[key3])
                    list1.append(list2)
        data = DataFrame(list1)
        data.columns = ['locus', "align_len", "sbj_len", "gaps", "coverage", "allele", "allele_name", "identity"]
        data.replace("", "-", inplace=True)
        data['allele_name'].replace("No hit found", "-", inplace=True)
        data2 = data[["locus", "allele"]]
        data2.fillna("-")
        names = list(data2["locus"])
        data3 = data2.set_index('locus')
        data3 = data3.T
        if "-" in list(data2["allele"]):
            allele = "no"
        else:
            allele = "all"
        if st != "un":
            data3['ST'] = st
        else:
            if allele == "all":
                data3['ST'] = "ST-new"
            elif allele == "no":
                data3['ST'] = "unknown"
        data3['sample'] = sample
        ss = ["sample", "ST"] + names
        data3.to_csv(prefix + ".mlst.ST.xls", sep='\t', header=True, index=False, columns=ss)
        del data['allele']
        data.to_csv(prefix + ".mlst.detail.xls", sep='\t', header=True, index=False,
                    columns=["locus", "align_len", "sbj_len", "allele_name", "gaps", "coverage", "identity"])


def main():
    parser = OptionParser()
    parser.add_option('--d', dest='data', metavar='[data json file]')
    parser.add_option('--s', dest='sample', metavar='[sample name]')
    parser.add_option("--o", dest="prefix", metavar="[out file prefix]")
    (options,args) = parser.parse_args()
    if not options.data or not options.sample or not options.prefix:
        print "python mlst_stat.py --d data --s sample_name --o prefix"
        return
    chang_result(options.data,  options.sample, options.prefix)

if __name__=='__main__':
    main()
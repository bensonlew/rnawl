# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re, os
import shutil
import copy
import pandas as pd


def filter_snp(martix, sample_list, step, out, ref):
    snps = pd.read_table(martix, sep='\t', header=0)
    ###只保留含参考在内的两个碱基
    snps2 = snps[snps.apply(lambda x: cent_fun(x) == 2, axis=1)]
    snps2 = snps2.reset_index(drop=True)
    samples = sample_list.split(";")
    data2 = snps2.ix[:, ["Site", ref, "Seq"]]
    for sample in samples:
        data = snps2.ix[:, ["Site", ref, sample]]
        for i in range(0, data.shape[0] - 1):
            da = data.ix[i:i]
            da = da.reset_index(drop=True)
            print(da)
            if da.ix[0, ref] != da.ix[0, sample] or da.ix[0, sample] == 'N':
                da2 = data.ix[i + 1:i + 1]
                da2 = da2.reset_index(drop=True)
                print(da2)
                if (int(da2.ix[0, "Site"]) - int(da.ix[0, "Site"])) < int(step):
                    data.ix[i + 1, 2] = "N"
        data = data[sample]
        data2 = pd.concat([data2, data], axis=1)
    ss = copy.copy(samples)
    ss.append(ref)
    data3 = data2.ix[:, ss].T
    with open(out + ".snp.fasta", "w") as d:
        d.write(">{}\n{}\n".format(ref, "".join(data3.loc[ref])))
        for sample in samples:
            data4 = data2.ix[:, [ref, sample]]
            data4 = data4[data4.apply(lambda x:function(x[ref], x[sample]) == 1, axis=1)]
            print(data4.shape[0])
            num = data4.shape[0]
            d.write(">{}\n{}\n".format(sample + "[" + str(num) + "]", "".join(data3.loc[sample])))
    data2.to_csv(out + ".martix.xls", sep='\t', header=True, index=False)


def cent_fun(values):
    all = []
    for v in values[2:]:
        all.append(v)
    return int(len(set(all)))

def function(a,b):
    if a != b and b != "N":
        return 1
    else:
        return 0

def main():
    parser = OptionParser()
    parser.add_option('--m', dest='martix',metavar='[snp martix file]')
    parser.add_option("--s", dest="sample_list", metavar="[samples of list]")
    parser.add_option("--l", dest="length", metavar="[length of two snp]")
    parser.add_option("--o", dest="prefix", metavar="[outfile prefix]")
    parser.add_option("--r", dest="ref", metavar="[ref name]")
    (options, args) = parser.parse_args()
    if not options.martix or not options.sample_list or not options.length or not options.prefix or not options.ref:
        print "python filter_snp.py --m martix --s sample_list --l length --r ref --o prefix "
        return
    filter_snp(options.martix, options.sample_list, options.length, options.prefix, options.ref)


if __name__=='__main__':
    main()
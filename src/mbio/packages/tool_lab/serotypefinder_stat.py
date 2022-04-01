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
        list_o = []
        list_h = []
        for key, value in data['serotypefinder']['results'].items():
            if key == "O_type":
                if value == "No hit found":
                   continue
                else:
                   for key2, value2 in value.items():
                      list1 = []
                      list1.append(key2)
                      for key3 in sorted(value2.keys()):
                         if key3 == "positions_in_contig":
                            ss =value2[key3].split("..")
                            for i in ss:
                               list1.append(i)
                         else:
                            list1.append(value2[key3])
                      list_o.append(list1)
            if key == "H_type":
                if value == "No hit found":
                    continue
                else:
                    for key2, value2 in value.items():
                        list2 = []
                        list2.append(key2)
                        for key3 in sorted(value2.keys()):
                            if key3 == "positions_in_contig":
                                ss =value2[key3].split("..")
                                for	i in ss:
                                    list2.append(i)
                            else:
                                list2.append(value2[key3])
                    list_h.append(list2)
        if len(list_o) >0 and len(list_h) >0:
            data_o = DataFrame(list_o)
            data_h = DataFrame(list_h)
            data_o.columns = ['locus', "HSP_length", "accession", "cong_id", "coverage", "gene","hit_id", "identity", "position_in_ref","start","end","serotype","template_length"]
            data_h.columns = ['locus', "HSP_length", "accession", "cong_id", "coverage", "gene","hit_id", "identity", "position_in_ref","start","end","serotype","template_length"]
            o_type = ";".join(list(set(data_o["serotype"])))
            o_gene = ";".join(list(set(data_o["gene"])))
            h_type = ";".join(list(set(data_h["serotype"])))
            h_gene = ";".join(list(set(data_h["gene"])))
            serotypes =''
            if len(list(set(data_o["serotype"]))) >=2 or len(list(set(data_h["serotype"]))) >=2:
                serotypes = "Multiple serotypes"
            else:
                serotypes = o_type+":"+h_type
            with open(prefix+".SerotypeStat.xls","w") as f:
                f.write("Sample\to_type\to_gene\th_type\th_gene\tserotypes\n")
                f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sample,o_type,o_gene,h_type,h_gene,serotypes))
            data = pd.concat([data_o,data_h])
            data = data[['locus',"serotype","gene",'start','end','accession',"identity",'template_length',"HSP_length"]]
            data.to_csv(prefix + ".SerotypeDetail.xls", sep='\t', header=True, index=False)
        elif len(list_h) >0 and len(list_o) ==0:
            data_h = DataFrame(list_h)
            data_h.columns = ['locus', "HSP_length", "accession", "cong_id", "coverage", "gene", "hit_id", "identity",
                              "position_in_ref", "start", "end", "serotype", "template_length"]
            o_type = "-"
            o_gene = "-"
            h_type = ";".join(list(set(data_h["serotype"])))
            h_gene = ";".join(list(set(data_h["gene"])))
            serotypes = ''
            if len(list(set(data_h["serotype"]))) >= 2:
                serotypes = "Multiple serotypes"
            else:
                serotypes = h_type
            with open(prefix + ".SerotypeStat.xls", "w") as f:
                f.write("Sample\to_type\to_gene\th_type\th_gene\tserotypes\n")
                f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, o_type, o_gene, h_type, h_gene, serotypes))
            data = data_h[
                ['locus', "serotype", "gene", 'start', 'end', 'accession', "identity", 'template_length', "HSP_length"]]
            data.to_csv(prefix + ".SerotypeDetail.xls", sep='\t', header=True, index=False)
        elif len(list_h) ==0 and len(list_o) >0:
            data_o = DataFrame(list_o)
            data_o.columns = ['locus', "HSP_length", "accession", "cong_id", "coverage", "gene", "hit_id", "identity",
                              "position_in_ref", "start", "end", "serotype", "template_length"]
            o_type = ";".join(list(set(data_o["serotype"])))
            o_gene = ";".join(list(set(data_o["gene"])))
            h_type = "-"
            h_gene = "-"
            serotypes = ''
            if len(list(set(data_o["serotype"]))) >= 2:
                serotypes = "Multiple serotypes"
            else:
                serotypes = o_type
            with open(prefix + ".SerotypeStat.xls", "w") as f:
                f.write("Sample\to_type\to_gene\th_type\th_gene\tserotypes\n")
                f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, o_type, o_gene, h_type, h_gene, serotypes))
            data = data_o[
                ['locus', "serotype", "gene", 'start', 'end', 'accession', "identity", 'template_length', "HSP_length"]]
            data.to_csv(prefix + ".SerotypeDetail.xls", sep='\t', header=True, index=False)
        elif len(list_h) ==0 and len(list_o) ==0:
            o_type = "-"
            o_gene = "-"
            h_type = "-"
            h_gene = "-"
            serotypes = '-'
            with open(prefix + ".SerotypeStat.xls", "w") as f:
                f.write("Sample\to_type\to_gene\th_type\th_gene\tserotypes\n")
                f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, o_type, o_gene, h_type, h_gene, serotypes))

def main():
    parser = OptionParser()
    parser.add_option('--d', dest='data', metavar='[data json file]')
    parser.add_option('--s', dest='sample', metavar='[sample name]')
    parser.add_option("--o", dest="prefix", metavar="[out file prefix]")
    (options,args) = parser.parse_args()
    if not options.data or not options.sample or not options.prefix:
        print "python mlst_stat.py --d data --s sample_name --o prefix "
        return
    chang_result(options.data,  options.sample, options.prefix)

if __name__=='__main__':
    main()

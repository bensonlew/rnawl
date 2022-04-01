# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
import sys
import pandas as pd

def get_value(df, profile, head, des):
    filter_list = df["Gene_list"].split(",")
    result = profile.loc[profile.index.isin(filter_list)].sum()
    result.name = df.name
    result.columns = head
    if des:
        result["Description"] = df.tolist().pop()
    return result

if __name__ == '__main__':
    try:
        if len(sys.argv) != 4:
            sys.exit(1)
    except:
        raise Exception("python %s [category.genes.list] [gene.profile] [output.file.list]" % sys.argv[0])
    script, cate, gpro, out = sys.argv
    ca = cate.split(",")
    ou = out.split(",")
    data = pd.read_table(gpro, index_col=0)
    print "data read pass"
    heads = data.columns.tolist()
    for index,value in enumerate(ca):
        infile_path = value
        oufile_path = ou[index]
        cate_data = pd.read_table(infile_path, index_col=0)
        print "cate_data read pass"
        newhead = "\t".join([cate_data.head().index.name] + heads)
        if "class_stat" in infile_path or "family_stat" in infile_path:
            des = True
        else:
            des = False
        anno_result = cate_data.apply(get_value, axis=1, args=(data, newhead, des,))
        print "anno_result pass"
        anno_result.to_csv(oufile_path, sep="\t", index=True)


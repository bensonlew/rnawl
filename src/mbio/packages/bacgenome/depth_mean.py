# -*- coding: utf-8 -*-
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', metavar='[table]', required=True, help='depth table')
parser.add_argument("-n", metavar="[mean number]", required=False, help="depth number per sequence, default 100")
parser.add_argument('-o', metavar='[output]', required=True, help='outfilename')
args = parser.parse_args()

inputfile = args.i
out = args.o
number = int(args.n) if args.n else 100

def run():
    data = pd.read_table(inputfile, header=None)
    seq_map= {}
    for index,i in enumerate(data[0].drop_duplicates()):
        seq_map[i] = i.replace("_pilon", "")
    datac = data.groupby(0).count()
    datac['size'] = datac[1]//number
    # print datac
    output_data = pd.DataFrame()
    for i in datac.index:
        tmp_data = data.loc[data[0]==i,]
        for j in range(number):
            start = j * datac.loc[i, "size"]
            end = (j+1) * datac.loc[i, "size"]
            if j + 1 == number or datac.loc[i, "size"] == 0:
                mean_depth = tmp_data.loc[start:, 2].mean()
                rang_str = str(start+1) + "-"
            else:
                mean_depth = tmp_data.iloc[start:end, 2].mean()
                rang_str = "%s-%s" % (start+1, end)
            #print mean_depth
            # output_data.append(pd.DataFrame([i, rang_str, mean_depth]),ignore_index=True)
            output_data = output_data.append(pd.DataFrame({1:[seq_map[i]], 2:[rang_str], 3: [mean_depth]}))
    # print output_data
    output_data.to_csv(out, sep="\t", header=None, index=None)

if __name__ == '__main__':
    run()

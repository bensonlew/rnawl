import pandas as pd



def main(args):
    f = pd.read_table(args.input)
    bsjreads = []
    for i in range(0,len(f)):
        a = max(f['junction_reads_x'].fillna(0)[i], f['junction_reads_y'].fillna(0)[i])
        bsjreads.append(a)

    f['junction_reads'] = bsjreads
    f.to_csv(args.output,index = False,sep = '\t')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='BSJreads')
    parser.add_argument('-i', action='store', required=True,
                        help='input file', dest='input')
    parser.add_argument('-o', action='store', required=True,
                        help='output file ', dest='output')

    args = parser.parse_args()
    main(args)
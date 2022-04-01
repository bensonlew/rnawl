# -*- coding: UTF-8 -*-
import pandas as pd
import re

def make_dict(GTF):
    chrom_id_start_end = dict()
    with open(GTF) as f:  # 打开gtf文件
        for lines in f:
            if not lines.startswith('#'):
                (chr, source, type, start, end, score, strand, phase, attribute_string) = lines.rstrip().split('\t')
                gene_id = re.findall(r'gene_id\s\"(.*?)\";', lines)
                gene_id_str = str(gene_id[0])
                # f1 = attribute_string.split(';')  # 将attribute_string按分号分割
                # f2 = f1[0].split(' ')  # 将第一列按空格分割
                if chr not in chrom_id_start_end:
                    chrom_id_start_end[chr] = {gene_id_str: [int(start), int(end)]}
                else:
                    if gene_id_str not in chrom_id_start_end[chr]:
                        chrom_id_start_end[chr].update({gene_id_str: [int(start), int(end)]})
                    else:
                        chrom_id_start_end[chr][gene_id_str].append(int(start))
                        chrom_id_start_end[chr][gene_id_str].append(int(end))  # 生成染色体：基因id：[起始，终止]
                        chrom_id_start_end[chr][gene_id_str].sort()
        return chrom_id_start_end


def circ_type(chrom_id_start_end,merge,output):
    f = pd.read_table(merge)
    chrome_none = list()
    for i in range(0, len(f)):
        chrome = str(f['chr'][i])
        if chrome not in chrom_id_start_end:
            chrome_none.append(chrome)
    new_f = f[~f['chr'].isin(chrome_none)]
    circRNA_type = []
    host_id = []
    for i in list(new_f.index):
        chrome = str(new_f['chr'][i])
        circ_start = new_f['circRNA_start'][i]
        circ_end = new_f['circRNA_end'][i]
        print(circ_start)
        print(circ_end)

        id_start = filter(lambda x: circ_start in x[1], chrom_id_start_end[chrome].items())   #circ_start和circ_end均在外显子上
        id_end = filter(lambda x: circ_end in x[1], chrom_id_start_end[chrome].items())
        a = list(id_start)
        b = list(id_end)
        if a and b:  #如果文件不为空
            gene_id_start1 = list()
            gene_id_end1 = list()
            for (key1, value1) in a:
                gene_id_start1.append(key1)
            for (key2, value2) in b:
                gene_id_end1.append(key2)
            id_list1 = list(set(gene_id_start1).intersection(set(gene_id_end1)))  #将circ_start所在的gene_id与circ_end所在的gene_id取交集
            if id_list1:
                id_list1 = ",".join(str(x) for x in id_list1)
                host_id.append(id_list1)
                circRNA_type.append('exon')
            else:
                circRNA_type.append('intergenic_region')
                host_id.append('')
        if a and not b: #如果circ_start在外显子上，circ_end不在外显子上
            id_end2 = filter(lambda x: x[1][0] < circ_end < x[1][-1], chrom_id_start_end[chrome].items())  #找到circ_end所在的间区
            d = list(id_end2)
            if d:
                gene_id_end2 = list()
                gene_id_start2 = list()
                for (key3, value3) in d:
                    gene_id_end2.append(key3)
                for (key4, value4) in a:
                    gene_id_start2.append(key4)
                id_list2 = list(set(gene_id_start2).intersection(set(gene_id_end2)))  # 将circ_start所在的gene_id与circ_end所在的gene_id取交集
                if id_list2:
                    id_list2 = ",".join(str(x) for x in id_list2)
                    circRNA_type.append('intron')
                    host_id.append(id_list2)
                else:
                    circRNA_type.append('intergenic_region')
                    host_id.append('')
            else:
                circRNA_type.append('intergenic_region')
                host_id.append('')
        if not a and b:
            id_start3 = filter(lambda x: x[1][0] < circ_start < x[1][-1], chrom_id_start_end[chrome].items())
            e = list(id_start3)
            if e:
                gene_id_start3 = list()
                gene_id_end3 = list()
                for (key5, value5) in e:
                    gene_id_start3.append(key5)
                for (key6, value6) in b:
                    gene_id_end3.append(key6)
                id_list3 = list(set(gene_id_start3).intersection(set(gene_id_end3)))  # 将circ_start所在的gene_id与circ_end所在的gene_id取交集
                if id_list3:
                    id_list3 = ",".join(str(x) for x in id_list3)
                    circRNA_type.append('intron')
                    host_id.append(id_list3)
                else:
                    circRNA_type.append('intergenic_region')
                    host_id.append('')
            else:
                circRNA_type.append('intergenic_region')
                host_id.append('')

        if  not a and not b:
            id_end4 = filter(lambda x: x[1][0] < circ_end < x[1][-1], chrom_id_start_end[chrome].items())
            id_start4 = filter(lambda x: x[1][0] < circ_start < x[1][-1], chrom_id_start_end[chrome].items())
            k = list(id_end4)
            g = list(id_start4)
            if k and g:
                gene_id_end4 = list()
                gene_id_start4 = list()
                for (key7, value7) in k:
                    gene_id_end4.append(key7)
                for (key8, value8) in g:
                    gene_id_start4.append(key8)
                id_list4 = list(set(gene_id_start4).intersection(set(gene_id_end4)))
                if id_list4:
                    id_list4 = ",".join(str(x) for x in id_list4)
                    circRNA_type.append('intron')
                    host_id.append(id_list4)
                else:
                    circRNA_type.append('intergenic_region')
                    host_id.append('')
            else:
                circRNA_type.append('intergenic_region')
                host_id.append('')
    print(circRNA_type)
    print(host_id)
    new_f['circRNA_type'] = circRNA_type
    new_f['host_gene_id'] = host_id
    new_f_none = new_f.fillna("")
    new_f_none.to_csv(output,index=False,sep='\t')



def main(args):
    # GTF = '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/data/Homo_sapiens.GRCh38.96.gtf'
    make_dict(args.gtf)
    chrom_id_start_end = make_dict(args.gtf)
    # merge = '/mnt/ilustre/users/sanger-dev/workspace/20190924/Single_signal_1746_8723/Signal/output/circ_merge_signal.txt'
    circ_type(chrom_id_start_end,args.signal,args.type)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='circRNA_type')
    parser.add_argument('-g', action='store', required=True,
                        help='GTF', dest='gtf')
    parser.add_argument('-s', action='store', required=True,
                        help='merge-signal', dest='signal')
    parser.add_argument('-o', action='store', required=True,
                        help='signal-type ', dest='type')

    args = parser.parse_args()

    main(args)




# -*- coding: utf-8 -*-
# __author__ = "qiuping"
# last_modify:20160616


def get_diff_matrix(matrix, diff_list, output):
    with open(matrix, 'rb') as m, open(diff_list, 'rb') as d, open(output, 'wb') as w:
        mline = m.readline()
        w.write(mline)
        diff_gene = []
        for i in d.readlines():
            diff_gene.append(i.strip('\n'))
        while True:
            mline = m.readline()
            if not mline:
                break
            if mline.strip('\n').split('\t')[0] in diff_gene:
                w.write(mline)


def get_diff_list(edgerfile, output, fc=2, pvalue_padjust=None,diff_ci=0.05):
    """根据fc和diff_ci过滤差异基因/转录本"""

    import math
    from math import log
    from math import pow
    with open(edgerfile, 'rb') as r, open(output, 'wb') as w:
        line = r.readline()
        while True:
            line = r.readline()
            if not line:
                break
            line = line.strip('\n').split('\t')
            if pvalue_padjust == 'padjust':
                print 'edgerfile{}开始过滤差异基因'.format(edgerfile)
                if float(line[4]) <= float(diff_ci):
                    
                    if float(fc)>=1:
                        if pow(2,float(line[1])) > float(fc) or pow(2,float(line[1])) <= (float(1)/float(fc)):
                            w.write('%s\n' % line[0])
                    if float(fc)<1:
                        if pow(2,float(line[1])) < float(fc) or pow(2,float(line[1])) >= (float(1)/float(fc)):
                            w.write('%s\n' % line[0])
            elif pvalue_padjust == 'pvalue':
                
                if float(line[3]) <= float(diff_ci):
                    if float(fc)>=1:
                        if pow(2,float(line[1])) > float(fc) or pow(2,float(line[1])) <= (float(1)/float(fc)):
                            w.write('%s\n' % line[0])
                    if float(fc)<1:
                        if pow(2,float(line[1])) < float(fc) or pow(2,float(line[1])) >= (float(1)/float(fc)):
                            w.write('%s\n' % line[0])


def check_dispersion(genes, diff_num, diff_rate):
    dispersion = 0.1
    rate = diff_num / genes
    dis = rate / diff_rate
    if dis > 1:
        dispersion += (dis - 1) * 0.2
    else:
        dispersion *= dis
    if dispersion > 0.95:
        dispersion = 0.95
    elif dispersion == 0:
        dispersion = 0.001
    return dispersion


def get_gene_list(fpkm_file, output):
    with open(fpkm_file, 'rb') as r, open(output, 'wb') as w:
        w.write('gene_id\tgene_id\n')
        lines = r.readlines()
        for line in lines[1:]:
            w.write('{}\t{}\n'.format(line.split('\t')[0], line.split('\t')[0]))
if __name__ == "__main__":
    edgerfile = "/mnt/ilustre/users/sanger-dev/workspace/20170524/Single_diff_transcript_12/DiffExp/edger_result/transcripts.counts.matrix.CD_vs_HGD.edgeR.DE_results"
    output = "/mnt/ilustre/users/sanger-dev/workspace/20170524/Single_diff_transcript_12/DiffExp/diff_list_dir/CD_vs_HGD.new"
    fc=1
    pvalue_padjust='padjust'
    diff_ci = 0.01
    get_diff_list(edgerfile, output, fc, pvalue_padjust, diff_ci)

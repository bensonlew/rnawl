# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
import argparse
import pandas as pd
from Bio import SeqIO
from mummer import Mummer


def get_order(samples, mash):
    '''使用mash计算各样本和参考样本间的距离
    :param samples: str, 命令行参数传入的样本列表，逗号分隔，参考在第一个
    :param mash: mash 的软件路径
    '''
    sp = samples.split(',')
    cmd = mash + ' dist ' + ' '.join(sp) + '> mash_out.xls'
    ret_code = os.system(cmd)
    if ret_code != 0 or not os.path.exists('mash_out.xls'):
        exit('mash 距离计算出错')
    ms = pd.read_csv('mash_out.xls', sep='\t', header=None).sort_values(by=2)
    return [ms.at[0, 0], ] + list(ms[1])


def extract_seq(ref, region=[]):
    extracted_to = os.path.basename(ref)
    for seq in SeqIO.parse(ref, 'fasta'):
        if seq.id == region[0]:
            start = int(region[1]) - 1
            end = int(region[2])
            new_seq = seq[start:end]
            SeqIO.write(new_seq, extracted_to, 'fasta')
            return extracted_to


def _main():

    mm = Mummer()
    mm.mum_path(args.mum_path)
    # 选择的参考基因组区域 [chr, start, end]
    region = args.region.split(',')
    output = []

    sp_order = get_order(args.samples, args.mash_path)
    for i in range(0, len(sp_order)-1):
        extracted_to = extract_seq(sp_order[i], region)
        output.append([os.path.splitext(extracted_to)[0], ] + region)
        region = mm.run('nucmer', extracted_to, sp_order[i+1],
                        one=True, s=True, region_mode=True)

    last_sp = os.path.splitext(os.path.basename(sp_order[-1]))[0]
    output.append([last_sp, ] + region)
    df = pd.DataFrame(data=output)
    df.columns = ['Sample', 'SeqID', 'Start', 'End']
    df.to_csv('regions.xls', sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-s', '--samples', required=True,)
    parser.add_argument('-r', '--region', required=True,)
    parser.add_argument('-m', '--mum_path', required=True,)
    parser.add_argument('-n', '--mash_path', required=True,)

    args = parser.parse_args()

    _main()

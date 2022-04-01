# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
import argparse
import pandas as pd
from itertools import combinations
from multiprocessing import Pool
from functools import wraps
# from Bio import SeqIO


class Mummer(object):
    def __init__(self):
        self._mum_path = ''
        self._prefix = 'out'

    def mum_path(self, value):
        self._mum_path = value + '/'

    def prefix(self, value):
        self._prefix = value

    # mummer 运行流程 nucmer -> filter -> coords
    def run_mummer(self, mummer, ref, query):
        '''
        nucmer 或 promer 进行比对
        --maxmatch 同时允许在 ref 和 query 基因组存在多个匹配
        '''
        cmd = self._mum_path + mummer + ' --maxmatch -p {} {} {}'.format(self._prefix, ref, query)
        print(cmd)
        ret_code = os.system(cmd)
        if ret_code:
            print(mummer + '运行出错')

    def run_filter(self, one=False):
        '''
        按条件对结果过滤
        '''
        cmd = self._mum_path + 'delta-filter '
        if one:
            # 1-to-1 alignment
            cmd += '-1 {}.delta > {}.deltaf'.format(self._prefix, self._prefix)
        else:
            # many-to-many
            cmd += '-m {}.delta > {}.deltaf'.format(self._prefix, self._prefix)

        ret_code = os.system(cmd)
        if ret_code:
            print('filter 步骤运行出错')

    def run_coords(self):
        # 一共13列
        # [S1] [E1] [S2] [E2] [LEN 1] [LEN 2] [% IDY]
        # [LEN R] [LEN Q] [COV R] [COV Q] [TAGS]
        cmd = self._mum_path +\
            'show-coords -rclTH {}.deltaf > {}.coords'.format(self._prefix, self._prefix)

        ret_code = os.system(cmd)
        if ret_code:
            print('coords 步骤运行出错')

    # 结果格式化输出
    def lm(self, row, pre, super_loc):
        if row[0] != pre[0] or row[1] - pre[2] > 1:
            if row[0] != pre[0]:
                pre[0:-1] = [0, 0, 0]
            pre[-1] += 1
            super_loc[str(pre[-1])] = [row[1], row[2]]
        if pre[2] < row[2]:
            pre[0:-1] = row
        if row[2] > super_loc[str(pre[-1])][-1]:
            super_loc[str(pre[-1])][-1] = row[2]
        return str(pre[-1])

    def get_super_block(self, df, sp, col):
        print('开始获取' + sp + ' super信息')
        cols = [col + i for i in ['_scaf', '_start', '_end']]
        df2 = df.sort_values(by=cols)
        super_loc = {}
        pre = [0, 0, 0, 0]
        tmp_df = df2[cols].apply(self.lm, args=(pre, super_loc), axis=1)
        df[col + '_super_id'] = sp + tmp_df
        df[col + '_super_start'], df[col + '_super_end'] = tmp_df.apply(lambda x: super_loc[x]).str
        # df[col + '_super_loc'] = tmp_df.apply(lambda x: super_loc[x])
        print('获取' + sp + ' super信息完成')
        return df

    def set_out(self, ref, query, region_mode=False, s=False):
        # 格式化输出结果
        print('## 格式化输出结果')
        header = [
            'ref_start', 'ref_end', 'query_start', 'query_end',
            'ref_len', 'query_len', 'identity', 'ref_scaf', 'query_scaf'
        ]
        df = pd.read_csv(self._prefix + '.coords', header=None, sep='\t')
        select = list(range(7)) + [-2, -1]
        df = df.iloc[:, select]
        df.columns = header

        df['ref'] = ref
        df['query'] = query

        df['direction'] = '+'
        t = df['query_start'] > df['query_end']
        if any(t):
            df.loc[t, 'direction'] = '-'
            df.loc[t, ['query_start', 'query_end']] = df.loc[t, ['query_end', 'query_start']].values
        if s:
            p = ref + '-' + query
            df = self.get_super_block(df, p + '-r', 'ref')
            df = self.get_super_block(df, p + '-q', 'query')
            if region_mode:
                #m_index = df.loc[:, 'query_super_loc'].apply(
                #                                             lambda x: x[1] - x[0]
                #                                             ).sort_values().index[-1]
                m_index = (df['query_super_end'] - df['query_super_start']).sort_values().index[-1]
                region = list(df.loc[m_index, ['query_scaf', 'query_super_start', 'query_super_end']])
                return [region[0], region[1], region[2]]
        return df

    def save(self, df):
        df.to_csv(self._prefix + '.output.xls', sep='\t', index=None)

    def run(self, mummer, ref, query, region_mode, one=False, s=False, out='.'):
        ref_name = os.path.splitext(os.path.basename(ref))[0]
        query_name = os.path.splitext(os.path.basename(query))[0]
        self._prefix = out + '/' + ref_name + '-' + query_name

        self.run_mummer(mummer, ref, query)
        self.run_filter(one)
        self.run_coords()
        return self.set_out(ref_name, query_name, region_mode, s)


def wrap(params):
    return params[0].run(*params[1:])


def _main():
    mm = Mummer()
    mm.mum_path(args.mum_path)

    sp_list = args.query_list.split(',')
    if args.ref:
        sp_list = [args.ref, ] + sp_list
    pairs = []
    if args.circle_mode:
        pairs = [(sp_list[0], sp_list[i]) for i in range(1, len(sp_list))]
    else:
        pairs = [(s, s) for s in sp_list]
        pairs += combinations(sp_list, 2)
    run_params = [(mm, args.mummer, s[0], s[1], args.region_mode, args.one, args.super, args.out_path) for s in pairs]

    pools = Pool(len(sp_list))
    pret = pools.map_async(wrap, run_params)
    pools.close()
    pools.join()

    rets = pret.get()
    df = pd.concat(rets)
    df.to_csv(args.out_path + '/mummer.output.xls', sep='\t', index=None)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-m', '--mummer', default='nucmer',)
    parser.add_argument('-r', '--ref',)
    parser.add_argument('-l', '--query_list', required=True)
    parser.add_argument('-p', '--mum_path', required=True)
    parser.add_argument('-o', '--out_path', default='.')
    parser.add_argument('-s', '--super', default=False, action='store_true')
    parser.add_argument('-one', default=False, action='store_true')
    parser.add_argument('-region_mode', default=False, action='store_true')
    parser.add_argument('-circle_mode', default=False, action='store_true')

    args = parser.parse_args()

    _main()

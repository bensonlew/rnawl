# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
import re
import argparse
import pandas as pd
from mbio.packages.statistical.metastat import *
from mbio.packages.bac_comp_genome.groups_CI import *
from mbio.packages.statistical.twosample_CI import *


def get_n_l(df):
    upperCI = [x for x in df.columns if re.search(r'upperCI$', x)]
    lowerCI = [x for x in df.columns if re.search(r'lowerCI$', x)]
    mean = [x for x in df.columns if re.search(r'mean$', x)]
    upperCI.sort()
    lowerCI.sort()

    n = df.apply(lambda x: round((x[upperCI].values - x[lowerCI].values).max() / x[mean].max()), axis=1)
    n[n<1] = 1
    df['n'] = n

    l = df.apply(lambda x: round(abs(x[lowerCI].min() / x['n']) + x[mean].max() + 3), axis=1)
    df['l'] = l

    return df


def twosample_ci(intable, testout, samp1, samp2, coverage, out, meth):
    meth = meth or 'DiffBetweenPropAsymptoticCC'
    if meth == "DiffBetweenPropAsymptoticCC":
        DiffBetweenPropAsymptoticCC(intable, testout, samp1, samp2,
                                    coverage, out)
    if meth == "DiffBetweenPropAsymptotic":
        DiffBetweenPropAsymptotic(intable, testout, samp1, samp2,
                                  coverage, out)
    if meth == "NewcombeWilson":
        NewcombeWilson(intable, testout, samp1, samp2, coverage, out)


def twogroup_ci(test, ci_file, testout, gfile, coverage):
    print('开始计算{}检验置信区间'.format(test))
    if test == 'student':
        student(testout, gfile, coverage)
    elif test == 'welch':
        welch(testout, gfile, coverage)
    elif test == 'mann':
        bootstrap(args.intable, gfile, coverage)
    if os.path.exists(ci_file):
        print('{}检验置信区间计算完成！'.format(test))
        return ci_file
    else:
        exit('{}检验计算置信区间出错！'.format(test))


def posthoc(testout, gfile, coverage, out, meth):
    meth = meth or 'tukeykramer'
    if meth == 'tukeykramer':
        tukeykramer(testout, gfile, coverage, out)
    if meth == 'gameshowell':
        gameshowell(testout, gfile, coverage, out)
    if meth == 'welchuncorrected':
        welchuncorrected(testout, gfile, coverage, out)
    if meth == 'scheffe':
        scheffe(testout, gfile, coverage, out)


def run_ci(args):
    test = args.test
    ci_file = test + '_CI.xls'
    if test in ['chi', 'fisher']:
        twosample_ci(args.intable, args.out, args.samp1, args.samp2,
                     args.coverage, ci_file, meth=None)
    elif test in ['kru_H', 'anova']:
        ci_file = test + '_CI'
        posthoc(args.out, args.gfile, args.coverage, ci_file, args.meth)
    elif test in ['student', 'welch', 'mann']:
        ci_file = twogroup_ci(test, ci_file, args.out,
                              args.gfile, args.coverage)


def meth_check(mtype, meth):
    meths = {
        'posthoc': [
            'tukeykramer', 'gameshowell', 'welchuncorrected', 'scheffe'
        ],
        'two': [
            "DiffBetweenPropAsymptoticCC",
            "DiffBetweenPropAsymptotic",
            "NewcombeWilson"
        ]
    }
    if not meth or (meth in meths[mtype]):
        return True
    else:
        exit()


def run_stat(args):
    test = args.test
    args.out = test + '_result.xls'

    if test in ['chi', 'fisher']:
        two_sample_test(args.intable, args.out, test,
                        args.samp1, args.samp2, ci=0.95,
                        test_type=args.testtype, mul_test='fdr')
    elif test in ['kru_H', 'anova']:
        mul_group_test(args.intable, args.out, test + '_box.xls',
                       args.gfile, test, mul_test='fdr')
    elif test in ['student', 'welch', 'mann']:
        two_group_test(args.intable, args.gfile, args.out,
                       test + '_box.xls', test, ci=0.95,  norm='T',
                       test_type=args.testtype, mul_test='fdr')
    r_script = '~/app/program/R-3.3.1/bin/Rscript run_' + test + '_test.r'

    return os.system(r_script)


def _main(args):
    files = os.listdir('.')
    for f in files:
        if re.match(r'.*result.xls|.*CI.*xls|.*box.xls|.*test.r', f):
            os.remove(f)
    r_code = run_stat(args)
    if r_code == 0:
        run_ci(args)
    else:
        exit('运行{}检验出错'.format(args.test))

    files = os.listdir('.')
    dfs = []
    for f in files:
        if re.match(r'.*result.xls|.*CI.*xls', f):
            df = pd.read_csv(f, index_col=0, sep='\t')
            dfs.append(df)
    if args.norm == 'T':
        df_in = None
        if not args.samp1:
            df_in = pd.read_csv(args.intable, index_col=0, sep='\t')
            df_in = df_in / df_in.sum() * 100
            df_in.columns = map(lambda x: x + '-propotion', df_in.columns)
    else:
        df_in = pd.read_csv(args.intable, index_col=0, sep='\t')
    dfs.append(df_in)

    df = pd.concat(dfs, axis=1)
    df.index.name = 'function'
    
    if args.test in ['kru_H', 'anova']:
        df = get_n_l(df)

    df.to_csv('overall_result.xls', sep='\t')


if __name__ == '__main__':
    parse = argparse.ArgumentParser()

    parse.add_argument('-t', '--test', required=True,)
    parse.add_argument('-tp', '--testtype', default='two.side')
    parse.add_argument('-mt', '--multitest', default='fdr')
    parse.add_argument('-ci', '--ci', type=float, default=0.95,)
    parse.add_argument('-i', '--intable', required=True,)
    parse.add_argument('-s1', '--samp1',)
    parse.add_argument('-s2', '--samp2',)
    parse.add_argument('-g', '--gfile',)
    parse.add_argument('-c', '--coverage', type=float, default=0.95)
    parse.add_argument('-m', '--meth', default=None)
    parse.add_argument('-n', '--norm', default='T')

    args = parse.parse_args()

    _main(args)

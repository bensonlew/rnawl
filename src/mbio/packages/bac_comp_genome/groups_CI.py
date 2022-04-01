# -*- coding: utf-8 -*-
# by xieshichang
# modified from packages.statistical.twogroup_CI
# and packages.statistical.mul_posthoc
import math
import scipy
import itertools
import pandas as pd
from scipy.stats.distributions import t
from numpy import mean
from collections import Counter
from mbio.packages.statistical.QTable import QTable


def group_detail(groupfile, get_member=False, mul=False):
    group = pd.read_csv(groupfile, sep='\t', header=None, comment='#')

    group_num = Counter(group[1])
    
    if mul:
        N = len(group[1])  # 样本的个数
        g_num = len(group_num)  # 分组的个数
        dfN = g_num - 1
        dfD = N - g_num
        return N, dfN, dfD, group_num

    if not get_member:
        return group_num

    group_member = {}
    for gname in group_num:
        group_member[gname] = list(group[group[1] == gname][0])
    return group_num, group_member


def stat_info(statfile, gnames):
    stat = pd.read_csv(statfile, sep='\t', index_col=0)
    mean_dict = {}
    sd_dict = {}
    taxon_list = stat.index

    for gname in gnames:
        gmean = gname + '-mean'
        gsd = gname + '-sd'
        mean_dict[gname] = stat[gmean]
        sd_dict[gname] = stat[gsd]

    return mean_dict, sd_dict, taxon_list


def student(statfile, groupfile, coverage):
    group_num_dict = group_detail(groupfile)
    gnames = sorted(group_num_dict.keys())
    (mean_dict, sd_dict, taxon_list) = stat_info(statfile, gnames)

    with open('student_CI.xls', 'w') as w:
        w.write('\teffectsize\tlowerCI\tupperCI\n')
        for tx in taxon_list:
            meanG1 = mean_dict[gnames[0]][tx]
            meanG2 = mean_dict[gnames[1]][tx]
            dp = meanG1 - meanG2
            varG1 = (sd_dict[gnames[0]][tx]**2)
            varG2 = (sd_dict[gnames[1]][tx]**2)
            n1 = group_num_dict[gnames[0]]
            n2 = group_num_dict[gnames[1]]

            dof = n1 + n2 - 2
            pooledVar = ((n1 - 1)*varG1 + (n2 - 1)*varG2) / dof
            sqrtPooledVar = math.sqrt(pooledVar)
            denom = sqrtPooledVar * math.sqrt(1.0/n1 + 1.0/n2)
            tCritical = t.isf(0.5 * (1.0-coverage), dof)
            lowerCI = dp - tCritical*denom
            upperCI = dp + tCritical*denom

            w.write('{}\t{}\t{}\t{}\n'.format(tx, dp, lowerCI, upperCI))


def welch(statfile, groupfile, coverage):
    group_num_dict = group_detail(groupfile)
    gnames = sorted(group_num_dict.keys())
    (mean_dict, sd_dict, taxon_list) = stat_info(statfile, gnames)

    with open('welch_CI.xls', 'w') as w:
        w.write('\teffectsize\tlowerCI\tupperCI\n')
        for tx in taxon_list:
            meanG1 = mean_dict[gnames[0]][tx]
            meanG2 = mean_dict[gnames[1]][tx]
            dp = meanG1 - meanG2
            varG1 = (sd_dict[gnames[0]][tx]**2)
            varG2 = (sd_dict[gnames[1]][tx]**2)
            n1 = group_num_dict[gnames[0]]
            n2 = group_num_dict[gnames[1]]

            normVarG1 = varG1 / n1
            normVarG2 = varG2 / n2
            unpooledVar = normVarG1 + normVarG2
            sqrtUnpooledVar = math.sqrt(unpooledVar)
            dof = (unpooledVar**2) / ((normVarG1**2) /
                                      (n1-1) + (normVarG2**2)/(n2-1))
            tCritical = t.isf(0.5 * (1.0 - coverage), dof)
            lowerCI = dp - tCritical*sqrtUnpooledVar
            upperCI = dp + tCritical*sqrtUnpooledVar
            w.write('{}\t{}\t{}\t{}\n'.format(tx, dp, lowerCI, upperCI))


def bootstrap(intable, groupfile, coverage):
    group_num_dict, group_member_dict = group_detail(
        groupfile, get_member=True)
    gnames = sorted(group_num_dict.keys())
    intable = pd.read_csv(intable, sep='\t', index_col=0)

    profile = (intable + 0.0) / intable.sum()

    scipy.random.seed(1234)
    with open('mann_CI.xls', 'w') as w:
        w.write('\teffectsize\tlowerCI\tupperCI\n')
        for index in profile.index:
            distribution = []
            for _ in xrange(0, 999):
                samplesGroup = {}
                for gname in gnames:
                    sampleSize = group_num_dict[gname]
                    samples = group_member_dict[gname]
                    choices = scipy.random.randint(0, sampleSize, sampleSize)
                    samplesGroup[gname] = profile.loc[index, samples][choices]
                diffOfMeanProp = samplesGroup[gnames[0]].mean() -\
                    samplesGroup[gnames[1]].mean()
                distribution.append(diffOfMeanProp*100)
            dp = profile.loc[index, group_member_dict[gnames[0]]].mean() -\
                profile.loc[index, group_member_dict[gnames[1]]].mean()
            distribution.sort()
            dp *= 100
            lowerCI = distribution[max(
                0, int(math.floor(0.5*(1.0-coverage)*len(distribution))))]
            upperCI = distribution[min(
                len(distribution) - 1,
                int(math.ceil((coverage+0.5*(1.0-coverage))*len(distribution)))
                )]
            w.write('{}\t{}\t{}\t{}\n'.format(index, dp, lowerCI, upperCI))


# 多组posthoc test计算CI
def scheffe(statfile, groupfile, coverage, outfile):
    (N, dfN, dfD, group_num_dict) = group_detail(groupfile, mul=True)
    gnames = group_num_dict.keys()
    (mean_dict, sd_dict, taxon_list) = stat_info(statfile, gnames)
    two_hoc = list(itertools.combinations(gnames, 2))
    for one in two_hoc:
        g = list(one)
        g.sort()
        groups = '-'.join(g)
        with open(outfile + '_scheffe_%s.xls' % groups, 'w') as w:
            cv = dfN*distributions.f.ppf(coverage, dfN, dfD)
            w.write('\t%s_effectsize\t%s_lowerCI\t%s_upperCI\t%s_pvalue\n' %
                    (groups, groups, groups, groups))
            for tx in taxon_list:
                # calculate within group variance
                withinGroupVar = 0
                for name in group_num_dict.keys():
                    withinGroupVar += (group_num_dict[name] -
                                       1)*(sd_dict[name][tx]**2)
                withinGroupVar /= dfD
                withinGroupStdDev = math.sqrt(withinGroupVar)
                if withinGroupVar == 0:
                    # degenerate case: within group variance is zero; set to 1e-6.
                    withinGroupVar = 1e-6
                es = mean_dict[g[0]][i] - mean_dict[g[1]][tx]
                invSampleSize = 1.0 / \
                    (group_num_dict[g[0]]) + 1.0/(group_num_dict[g[1]])
                Fs = (es * es) / (withinGroupVar*invSampleSize)
                pValue = 1.0 - distributions.f.cdf(Fs / dfN, dfN, dfD)
                # confidence interval
                confInter = math.sqrt(cv*invSampleSize)*withinGroupStdDev
                lowerCI = es - confInter
                upperCI = es + confInter
                w.write('%s\t%s\t%s\t%s\t%s\n' %
                        (tx, es, lowerCI, upperCI, pValue))


def welchuncorrected(statfile, groupfile, coverage, outfile):
    (N, dfN, dfD, group_num_dict) = group_detail(groupfile, mul=True)
    gnames = group_num_dict.keys()
    (mean_dict, sd_dict, taxon_list) = stat_info(statfile, gnames)
    # the numbers of post-hoc test
    two_hoc = list(itertools.combinations(gnames, 2))
    for one in two_hoc:
        g = list(one)
        g.sort()
        groups = '-'.join(g)
        with open(outfile + '_welchuncorrected_%s.xls' % groups, 'w') as w:
            cv = dfN*distributions.f.ppf(coverage, dfN, dfD)
            w.write('\t%s_effectsize\t%s_lowerCI\t%s_upperCI\t%s_pvalue\n' %
                    (groups, groups, groups, groups))
            for tx in taxon_list:
                meanG1 = mean_dict[g[0]][tx]
                meanG2 = mean_dict[g[1]][tx]
                dp = meanG1 - meanG2
                varG1 = sd_dict[g[0]][tx]**2
                varG2 = sd_dict[g[1]][tx]**2
                n1 = group_num_dict[g[0]]
                n2 = group_num_dict[g[1]]
                normVarG1 = varG1 / n1
                normVarG2 = varG2 / n2
                unpooledVar = normVarG1 + normVarG2
                sqrtUnpooledVar = math.sqrt(unpooledVar)
                if unpooledVar != 0:
                    # p-value
                    T_statistic = -1 * abs(meanG1 - meanG2) / sqrtUnpooledVar
                    dof = unpooledVar**2 / \
                        ((normVarG1**2)/(n1-1) + (normVarG2**2)/(n2-1))
                    pValue = t.cdf(T_statistic, dof) * 2
                    # CI
                    tCritical = t.isf(0.5 * (1.0-coverage), dof)
                    # 0.5 factor accounts from symmetric nature of distribution
                    lowerCI = dp - tCritical*sqrtUnpooledVar
                    upperCI = dp + tCritical*sqrtUnpooledVar
                else:
                    if meanG1 != meanG2:
                        pValue = 0.0
                        # the difference (at least according to these samples) must be true as there is no variance
                    else:
                        pValue = 0.5
                    lowerCI = dp
                    upperCI = dp
                w.write('%s\t%s\t%s\t%s\t%s\n' %
                        (tx, dp, lowerCI, upperCI, pValue))


def tukeykramer(statfile, groupfile, coverage, outfile, preferences=None):
    qtable = QTable(preferences)
    (N, dfN, dfD, group_num_dict) = group_detail(groupfile, mul=True)
    gnames = group_num_dict.keys()
    (mean_dict, sd_dict, taxon_list) = stat_info(statfile, gnames)
    k = len(group_num_dict)
    q_cv = qtable.cv(1.0-coverage, k, dfD)
    cv001 = qtable.cv(0.001, k, dfD)
    cv01 = qtable.cv(0.01, k, dfD)
    cv02 = qtable.cv(0.02, k, dfD)
    cv05 = qtable.cv(0.05, k, dfD)
    cv1 = qtable.cv(0.1, k, dfD)
    two_hoc = list(itertools.combinations(gnames, 2))
    for one in two_hoc:
        g = list(one)
        g.sort()
        groups = '-'.join(g)
        with open(outfile + '_tukeykramer_%s.xls' % groups, 'w') as w:
            w.write('\t%s_effectsize\t%s_lowerCI\t%s_upperCI\t%s_pvalue\n' %
                    (groups, groups, groups, groups))
            for tx in taxon_list:
                # calculate within group variance
                withinGroupVar = 0
                for name in group_num_dict.keys():
                    withinGroupVar += (group_num_dict[name] -
                                       1)*(sd_dict[name][tx]**2)
                withinGroupVar /= dfD
                withinGroupStdDev = math.sqrt(withinGroupVar)
                if withinGroupStdDev == 0:
                    # degenerate case: within group variance is zero; set to 1e-6.
                    withinGroupStdDev = 1e-6
                sqrtInvSampleSize = math.sqrt(
                    (1.0/group_num_dict[g[0]] + 1.0/group_num_dict[g[1]]) / 2.0
                )
                meanG1 = mean_dict[g[0]][tx]
                meanG2 = mean_dict[g[1]][tx]
                es = meanG1 - meanG2
                qs = abs(es) / (withinGroupStdDev*sqrtInvSampleSize)
                if qs > cv001:
                    pValue = '< 0.001'
                elif qs > cv01:
                    pValue = '< 0.01'
                elif qs > cv02:
                    pValue = '< 0.05'   # < 0.02
                elif qs > cv05:
                    pValue = '< 0.05'
                elif qs > cv1:
                    pValue = '< 0.1'
                else:
                    pValue = '>= 0.1'
                confInter = q_cv * withinGroupStdDev * sqrtInvSampleSize
                lowerCI = es - confInter
                upperCI = es + confInter
                w.write('%s\t%s\t%s\t%s\t%s\n' %
                        (tx, es, lowerCI, upperCI, pValue))


def gameshowell(statfile, groupfile, coverage, outfile, preferences=None):
    qtable = QTable(preferences)
    (N, dfN, dfD, group_num_dict) = group_detail(groupfile, mul=True)
    gnames = group_num_dict.keys()
    (mean_dict, sd_dict, taxon_list) = stat_info(statfile, gnames)
    k = len(group_num_dict)
    two_hoc = list(itertools.combinations(gnames, 2))
    for one in two_hoc:
        g = list(one)
        g.sort()
        groups = '-'.join(g)
        with open(outfile + '_gameshowell_%s.xls' % groups, 'w') as w:
            w.write('\t%s_effectsize\t%s_lowerCI\t%s_upperCI\t%s_pvalue\n' %
                    (groups, groups, groups, groups))
            for tx in taxon_list:
                meanG1 = mean_dict[g[0]][tx]
                meanG2 = mean_dict[g[1]][tx]
                # effect size
                es = meanG1 - meanG2
                varG1 = sd_dict[g[0]][tx]**2
                varG2 = sd_dict[g[1]][tx]**2
                n1 = group_num_dict[g[0]]
                n2 = group_num_dict[g[1]]
                vn1 = varG1 / n1
                vn2 = varG2 / n2
                if vn1 == 0:
                    vn1 = 1e-6
                if vn2 == 0:
                    vn2 = 1e-6
                df = (vn1 + vn2) * (vn1 + vn2)
                df /= (vn1*vn1)/(n1-1) + (vn2*vn2)/(n2-1)
                q_cv = qtable.cvInterpolate(1.0-coverage, k, df)
                cv001 = qtable.cvInterpolate(0.001, k, df)
                cv01 = qtable.cvInterpolate(0.01, k, df)
                cv02 = qtable.cvInterpolate(0.02, k, df)
                cv05 = qtable.cvInterpolate(0.05, k, df)
                cv1 = qtable.cvInterpolate(0.1, k, df)
                # calculate Games-Howell unequal variance adjustment
                varAdj = math.sqrt((vn1 + vn2) / 2.0)
                # p-value
                qs = abs(es) / varAdj
                if qs > cv001:
                    pValue = '< 0.001'
                elif qs > cv01:
                    pValue = '< 0.01'
                elif qs > cv02:
                    pValue = '< 0.05'  # < 0.02
                elif qs > cv05:
                    pValue = '< 0.05'
                elif qs > cv1:
                    pValue = '< 0.1'
                else:
                    pValue = '>= 0.1'
                # confidence interval
                confInter = q_cv * varAdj
                lowerCI = es - confInter
                upperCI = es + confInter
                w.write('%s\t%s\t%s\t%s\t%s\n' %
                        (tx, es, lowerCI, upperCI, pValue))

# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

import math
import itertools
from scipy.stats import distributions
from scipy.stats.distributions import t
from QTable import QTable


def group_detail(groupfile):
    with open(groupfile, 'r') as g:
        ginfo = g.readlines()
        N = len(ginfo) - 1  #sample numbers
        group_dict = {}
        group_num = {}
        for line in ginfo[1:]:
            line = line.split()
            group_dict[line[0]] = line[1]
        gnames = group_dict.values()
        for gname in gnames:
            group_num[gname] = gnames.count(gname)
        g_num = len(group_num)  #group type numbers
        dfN = g_num - 1
        dfD = N - g_num
        return N, dfN, dfD, group_num


def stat_info(statfile, groupfile):
    (N, dfN, dfD, group_num_dict) = group_detail(groupfile)
    with open(statfile, 'r') as s:
        shead = s.readline().strip('\n').split()
        # print shead
        mean_dict = {}
        sd_dict = {}
        taxon_list = []
        for gname in group_num_dict.keys():
            mean_dict[gname] = []
            sd_dict[gname] = []
        while True:
            sline = s.readline()
            sline_list = sline.strip('\n').split('\t')
            if not sline:
                break
            taxon_list.append(sline_list[0])
            for gname in group_num_dict.keys():
                gmean = gname + "-mean"
                for foo in shead:
                    foo_list = [foo]
                    if gmean in foo_list:
                        index_site = shead.index(foo)
                        # print index_site
                        mean_dict[gname].append(float(sline_list[index_site + 1].strip('\"')))
                        sd_dict[gname].append(float(sline_list[index_site + 2].strip('\"')))
        # print mean_dict, sd_dict, taxon_list
        return mean_dict, sd_dict, taxon_list
# stat_info("C:\\Users\\ping.qiu.MAJORBIO\\Desktop\\anova_result.xls", "C:\\Users\\ping.qiu.MAJORBIO\\Desktop\\anova_group")



def scheffe(statfile, groupfile, coverage, outfile):
    (mean_dict, sd_dict, taxon_list) = stat_info(statfile, groupfile)
    (N, dfN, dfD, group_num_dict) = group_detail(groupfile)
    two_hoc = list(itertools.combinations(group_num_dict.keys(), 2)) # the numbers of post-hoc test
    for one in two_hoc:
        g = list(one)
        g.sort()
        groups = '-'.join(g)
        with open(outfile + '_scheffe_%s.xls' % groups, 'w') as w:
            cv = dfN*distributions.f.ppf(coverage, dfN, dfD)
            w.write('\t%s_effectsize\t%s_lowerCI\t%s_upperCI\t%s_pvalue\n' % (groups, groups, groups, groups))
            for i in range(len(taxon_list)):
                # calculate within group variance
                withinGroupVar = 0
                for name in group_num_dict.keys():
                    withinGroupVar += (group_num_dict[name] -1)*(sd_dict[name][i]**2)
                withinGroupVar /= dfD
                withinGroupStdDev = math.sqrt(withinGroupVar)
                if withinGroupVar == 0:
                    withinGroupVar = 1e-6  #degenerate case: within group variance is zero; set to 1e-6.
                es = mean_dict[g[0]][i] - mean_dict[g[1]][i]
                invSampleSize = 1.0/(group_num_dict[g[0]]) + 1.0/(group_num_dict[g[1]])
                Fs = (es * es) / (withinGroupVar*invSampleSize)
                pValue = 1.0 - distributions.f.cdf(Fs / dfN, dfN, dfD)
                # confidence interval
                confInter = math.sqrt(cv*invSampleSize)*withinGroupStdDev
                lowerCI = es - confInter
                upperCI = es + confInter
                w.write('%s\t%s\t%s\t%s\t%s\n' % (taxon_list[i], '%0.4g' % es, '%0.4g' % lowerCI, '%0.4g' % upperCI, '%0.4g' % pValue))


def welchuncorrected(statfile, groupfile, coverage, outfile):
    (mean_dict, sd_dict, taxon_list) = stat_info(statfile, groupfile)
    (N, dfN, dfD, group_num_dict) = group_detail(groupfile)
    two_hoc = list(itertools.combinations(group_num_dict.keys(), 2)) # the numbers of post-hoc test
    for one in two_hoc:
        g = list(one)
        g.sort()
        groups = '-'.join(g)
        with open(outfile + '_welchuncorrected_%s.xls' % groups, 'w') as w:
            cv = dfN*distributions.f.ppf(coverage, dfN, dfD)
            w.write('\t%s_effectsize\t%s_lowerCI\t%s_upperCI\t%s_pvalue\n' % (groups, groups, groups, groups))
            for i in range(len(taxon_list)):
                meanG1 = mean_dict[g[0]][i]
                meanG2 = mean_dict[g[1]][i]
                dp = meanG1 - meanG2
                varG1 = sd_dict[g[0]][i]**2
                varG2 = sd_dict[g[1]][i]**2
                n1 = group_num_dict[g[0]]
                n2 = group_num_dict[g[1]]
                normVarG1 = varG1 / n1
                normVarG2 = varG2 / n2
                unpooledVar = normVarG1 + normVarG2
                sqrtUnpooledVar = math.sqrt(unpooledVar)
                if unpooledVar != 0:
                    # p-value
                    T_statistic = -1 * abs(meanG1 - meanG2) / sqrtUnpooledVar
                    dof = (unpooledVar*unpooledVar) / ( (normVarG1*normVarG1)/(n1-1) + (normVarG2*normVarG2)/(n2-1) )
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
                w.write('%s\t%s\t%s\t%s\t%s\n' % (taxon_list[i], '%0.4g' % dp, '%0.4g' % lowerCI, '%0.4g' % upperCI, '%0.4g' % pValue))


def tukeykramer(statfile, groupfile, coverage, outfile, preferences=None):
    qtable = QTable(preferences)
    (mean_dict, sd_dict, taxon_list) = stat_info(statfile, groupfile)
    (N, dfN, dfD, group_num_dict) = group_detail(groupfile)
    k = len(group_num_dict)
    q_cv = qtable.cv(1.0-coverage, k, dfD)
    cv001 = qtable.cv(0.001, k, dfD)
    cv01 = qtable.cv(0.01, k, dfD)
    cv02 = qtable.cv(0.02, k, dfD)
    cv05 = qtable.cv(0.05, k, dfD)
    cv1 = qtable.cv(0.1, k, dfD)
    two_hoc = list(itertools.combinations(group_num_dict.keys(), 2))
    for one in two_hoc:
        g = list(one)
        g.sort()
        groups = '-'.join(g)
        with open(outfile + '_tukeykramer_%s.xls' % groups, 'w') as w:
            w.write('\t%s_effectsize\t%s_lowerCI\t%s_upperCI\t%s_pvalue\n' % (groups, groups, groups, groups))
            for i in range(len(taxon_list)):
                # calculate within group variance
                withinGroupVar = 0
                for name in group_num_dict.keys():
                    withinGroupVar += (group_num_dict[name] -1)*(sd_dict[name][i]**2)
                withinGroupVar /= dfD
                withinGroupStdDev = math.sqrt(withinGroupVar)
                if withinGroupStdDev == 0:
                    withinGroupStdDev = 1e-6  #degenerate case: within group variance is zero; set to 1e-6.
                sqrtInvSampleSize = math.sqrt((1.0/group_num_dict[g[0]] + 1.0/group_num_dict[g[1]]) / 2.0)
                meanG1 = mean_dict[g[0]][i]
                meanG2 = mean_dict[g[1]][i]
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
                w.write('%s\t%s\t%s\t%s\t%s\n' % (taxon_list[i], '%0.4g' % es, '%0.4g' % lowerCI, '%0.4g' % upperCI, pValue))


def gameshowell(statfile, groupfile, coverage, outfile, preferences=None):
    qtable = QTable(preferences)
    (mean_dict, sd_dict, taxon_list) = stat_info(statfile, groupfile)
    (N, dfN, dfD, group_num_dict) = group_detail(groupfile)
    k = len(group_num_dict)
    two_hoc = list(itertools.combinations(group_num_dict.keys(), 2))
    for one in two_hoc:
        g = list(one)
        g.sort()
        groups = '-'.join(g)
        with open(outfile + '_gameshowell_%s.xls' % groups, 'w') as w:
            w.write('\t%s_effectsize\t%s_lowerCI\t%s_upperCI\t%s_pvalue\n' % (groups, groups, groups, groups))
            for i in range(len(taxon_list)):
                meanG1 = mean_dict[g[0]][i]
                meanG2 = mean_dict[g[1]][i]
                # effect size
                es = meanG1 - meanG2
                varG1 = sd_dict[g[0]][i]**2
                varG2 = sd_dict[g[1]][i]**2
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
                varAdj = math.sqrt( (vn1 + vn2) / 2.0)
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
                confInter = q_cv*varAdj
                lowerCI = es - confInter
                upperCI = es + confInter
                w.write('%s\t%s\t%s\t%s\t%s\n' % (taxon_list[i], '%0.4g' % es, '%0.4g' % lowerCI, '%0.4g' % upperCI, pValue))

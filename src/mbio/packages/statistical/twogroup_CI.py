# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
import math
import scipy
from scipy.stats.distributions import t
from numpy import mean


def group_detail(groupfile, get_member = False):
    with open(groupfile, 'r') as g:
        ginfo = g.readlines()
        group_dict = {}
        group_num = {}
        group_member = {}
        for line in ginfo[1:]:
            line = line.split()
            group_dict[line[0]] = line[1]
        gnames = group_dict.values()
        for gname in gnames:
            group_num[gname] = gnames.count(gname)
            group_member[gname] = []
            for member in group_dict.keys():
                if group_dict[member] == gname:
                    group_member[gname].append(member)
        if get_member:
            return group_num, group_member
        return group_num


def stat_info(statfile, groupfile):
    group_num_dict = group_detail(groupfile)
    with open(statfile, 'r') as s:
        shead = s.readline().strip('\n').split('\t')
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
                #for foo in shead:
                    #if gmean in foo:
                index_site = shead.index(gmean)
                mean_dict[gname].append(float(sline_list[index_site].strip('\"')))
                sd_dict[gname].append(float(sline_list[index_site + 1].strip('\"')))
        return mean_dict, sd_dict, taxon_list


def student(statfile, groupfile, coverage):
    """
    计算影响大小及置信区间
    """
    (mean_dict, sd_dict, taxon_list) = stat_info(statfile, groupfile)
    group_num_dict = group_detail(groupfile)
    group_names = sorted(group_num_dict.keys())
    length = len(taxon_list)
    with open('student_CI.xls', 'w') as w:
        w.write('\teffectsize\tlowerCI\tupperCI\n')
        for i in range(length):
            # calculate effect size, and CI
            meanG1 = mean_dict[group_names[0]][i]
            meanG2 = mean_dict[group_names[1]][i]
            dp = meanG1 - meanG2
            varG1 = (sd_dict[group_names[0]][i])**2
            varG2 = (sd_dict[group_names[1]][i])**2
            n1 = group_num_dict[group_names[0]]
            n2 = group_num_dict[group_names[1]]

            dof = n1 + n2 -2
            dof = n1 + n2 - 2
            pooledVar = ((n1 - 1)*varG1 + (n2 - 1)*varG2) / (n1 + n2 - 2)
            sqrtPooledVar = math.sqrt(pooledVar)
            denom = sqrtPooledVar * math.sqrt(1.0/n1 + 1.0/n2)
            tCritical = t.isf(0.5 * (1.0-coverage), dof)  # 0.5 factor accounts from symmetric nature of distribution
            lowerCI = dp - tCritical*denom
            upperCI = dp + tCritical*denom
            w.write('%s\t%s\t%s\t%s\n' % (taxon_list[i], '%0.4g' % (dp), '%0.4g' % (lowerCI), '%0.4g' % (upperCI)))


def welch(statfile, groupfile, coverage):
    (mean_dict, sd_dict, taxon_list) = stat_info(statfile, groupfile)
    group_num_dict = group_detail(groupfile)
    group_names = sorted(group_num_dict.keys())
    length = len(taxon_list)
    with open('welch_CI.xls', 'w') as w:
        w.write('\teffectsize\tlowerCI\tupperCI\n')
        for i in range(length):
            meanG1 = mean_dict[group_names[0]][i]
            meanG2 = mean_dict[group_names[1]][i]
            dp = meanG1 - meanG2

            varG1 = (sd_dict[group_names[0]][i])**2
            varG2 = (sd_dict[group_names[1]][i])**2
            n1 = group_num_dict[group_names[0]]
            n2 = group_num_dict[group_names[1]]

            normVarG1 = varG1 / n1
            normVarG2 = varG2 / n2
            unpooledVar = normVarG1 + normVarG2
            sqrtUnpooledVar = math.sqrt(unpooledVar)
            dof = (unpooledVar*unpooledVar) / ( (normVarG1*normVarG1)/(n1-1) + (normVarG2*normVarG2)/(n2-1) )
            # CI
            tCritical = t.isf(0.5 * (1.0-coverage), dof) # 0.5 factor accounts from symmetric nature of distribution
            lowerCI = dp - tCritical*sqrtUnpooledVar
            upperCI = dp + tCritical*sqrtUnpooledVar
            w.write('%s\t%s\t%s\t%s\n' % (taxon_list[i], '%0.4g' % (dp), '%0.4g' % (lowerCI), '%0.4g' % (upperCI)))

def bootstrap(profile, groupfile, coverage,method="T"):
    # (mean_dict, sd_dict, taxon_list) = stat_info(statfile, groupfile)
    group_num_dict, group_member_dict = group_detail(groupfile, get_member = True)
    group_names = sorted(group_num_dict.keys())
    profile_dict = {}
    # w = open('bootstrap_CI.xls', 'w')
    # w.write('\teffectsize\tlowerCI\tupperCI\n')
    with open(profile, 'r') as s:
        shead = s.readline().strip('\n').split('\t')
        group_index = {}
        sum = [0]
        # profile_dict['samples'] = {}
        for gname in group_names:
            # profile_dict[gname] = []
            # profile_dict['samples'][gname] = []
            group_index[gname] = []
            for name in group_member_dict[gname]:
                for i in xrange(0, len(shead)):
                #for i in xrange(0, len(shead) + 1):
                    if name == shead[i].split('-')[-1]:
                        #group_index[gname].append(i)
                        group_index[gname].append(i + 1)
                        sum.append(0)
        while True:
            sline = s.readline()
            if not sline: break
            sline_list = sline.strip('\n').split('\t')
            profile_dict[sline_list[0]] = {}
            for index in xrange(1,len(sline_list)):
                sum[index] += float(sline_list[index])
            for gname in group_names:
                profile_dict[sline_list[0]][gname] = []
                for index_num in group_index[gname]:
                    # profile_dict[gname].append(float(sline_list[index_num]))
                    profile_dict[sline_list[0]][gname].append(float(sline_list[index_num]))

    w = open('bootstrap_CI.xls', 'w')
    w.write('\teffectsize\tlowerCI\tupperCI\n')
    for annotation in profile_dict.keys():
        for gname in group_names:
            for i, index in enumerate(group_index[gname]):
                if method == "T":
                    if sum[index] != 0:
                        profile_dict[annotation][gname][i] /= sum[index]
        distribution = []
        for _ in xrange(0, 999):
            samplesGroup = {}
            for gname in group_names:
                sampleSize = group_num_dict[gname]
                choices = scipy.random.random_integers(0, sampleSize - 1, sampleSize)
                g_scipy = scipy.array(profile_dict[annotation][gname], copy = 0)
                samplesGroup[gname] = g_scipy[choices]
            diffOfMeanProp = mean(samplesGroup[group_names[0]]) - mean(samplesGroup[group_names[1]])
            dp = mean(profile_dict[annotation][group_names[0]]) - mean(profile_dict[annotation][group_names[1]])
            if method == "T":
                distribution.append(diffOfMeanProp * 100)
            else:
                distribution.append(diffOfMeanProp)
        distribution.sort()
        if method == "T":
            dp = dp * 100
        lowerCI = distribution[max(0, int(math.floor(0.5*(1.0-coverage)*len(distribution))))]
        upperCI = distribution[min(len(distribution)-1, int(math.ceil((coverage + 0.5*(1.0-coverage))*len(distribution))))]
        w.write('%s\t%s\t%s\t%s\n' % (annotation, '%0.4g' % (dp), '%0.4g' % (lowerCI), '%0.4g' % (upperCI)))

        '''
            distribution = []
            for _ in xrange(0, 999):
                samplesGroup = {}
                for gname in group_names:
                    sampleSize = group_num_dict[gname]
                    choices = scipy.random.random_integers(0, sampleSize - 1, sampleSize)
                    g_scipy = scipy.array(profile_dict[gname], copy = 0)
                    samplesGroup[gname] = g_scipy[choices]
                diffOfMeanProp = mean(samplesGroup[group_names[0]]) - mean(samplesGroup[group_names[1]])
                dp = mean(profile_dict[group_names[0]]) - mean(profile_dict[group_names[1]])
                distribution.append(diffOfMeanProp)
            distribution.sort()
            lowerCI = distribution[max(0, int(math.floor(0.5*(1.0-coverage)*len(distribution))))]
            upperCI = distribution[min(len(distribution)-1, int(math.ceil((coverage + 0.5*(1.0-coverage))*len(distribution))))]
            w.write('%s\t%s\t%s\t%s\n' % (sline_list[0], '%0.4g' % (dp), '%0.4g' % (lowerCI), '%0.4g' % (upperCI)))
        '''
    w.close()
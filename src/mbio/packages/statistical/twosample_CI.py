# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

import math
from scipy.stats import norm
from scipy.special import erfinv


def standardNormalCDF(z):
    return norm.cdf(z)


def zScore(area):
    if area == 0.9:
        return 1.6448536269514722
    if area == 0.95:
        return 1.959963984540054
    if area == 0.98:
        return 2.3263478740408408
    if area == 0.99:
        return 2.5758293035489004
    return erfinv(area) * math.sqrt(2.0)


def otu_info(otufile, sample1, sample2):
    with open(otufile, 'r') as o:
        head = o.readline().strip().split('\t')
        totalSeq1 = 0
        totalSeq2 = 0
        s1_list = []
        s2_list = []
        taxon_list = []
        while True:
            line = o.readline()
            if not line:
                break
            taxon_list.append(line.split('\t')[0])
            for sam in head:
                if sam == sample1:
                    index_1 = head.index(sam)
                    totalSeq1 += float((line.strip('\n').split('\t')[index_1]).strip('\"'))
                    s1_list.append(float((line.strip('\n').split('\t')[index_1]).strip('\"')))
                if sam == sample2:
                    index_2 = head.index(sam)
                    totalSeq2 += float((line.strip('\n').split('\t')[index_2]).strip('\"'))
                    s2_list.append(float((line.strip('\n').split('\t')[index_2]).strip('\"')))
        return taxon_list, s1_list, s2_list, totalSeq1, totalSeq2

# (taxon_list, s1_list, s2_list, totalSeq1, totalSeq2) = otu_info('C:\Users\ping.qiu.MAJORBIO\Desktop\otu_file.xls', 'sample2', 'sample3')
# print taxon_list,s1_list,s2_list,totalSeq1,totalSeq2
def DiffBetweenPropAsymptoticCC(otufile, statfile, sample1, sample2, coverage, outfile):
    (taxon_list, s1_list, s2_list, totalSeq1, totalSeq2) = otu_info(otufile, sample1, sample2)
    with open(outfile, 'w') as w, open(statfile, 'r') as s:
        head = s.readline().strip().split('\t')
        stat_taxon = []
        while True:
            line = s.readline()
            if not line:
                break
            stat_taxon.append(line.split('\t')[0])
        w.write('\teffectsize\tlowerCI\tupperCI\n')
        for name in stat_taxon:
            i = taxon_list.index(name)
            taxon = taxon_list[i]
            R1 = s1_list[i] / totalSeq1
            R2 = s2_list[i] / totalSeq2
            diff = R1 - R2
            stdErr = math.sqrt((R1*(1-R1)) / totalSeq1 + (R2*(1-R2)) / totalSeq2) + (1.0/totalSeq1 + 1.0/totalSeq2)/2
            offset = zScore(coverage) * stdErr
            lowerCI = diff - offset
            upperCI = diff + offset
            diff *= 100
            lowerCI *= 100
            upperCI *= 100
            w.write('%s\t%s\t%s\t%s\n' % (taxon, '%0.4g' % diff, '%0.4g' % lowerCI, '%0.4g' % upperCI))

def DiffBetweenPropAsymptotic(otufile, statfile, sample1, sample2, coverage, outfile):
    (taxon_list, s1_list, s2_list, totalSeq1, totalSeq2) = otu_info(otufile, sample1, sample2)
    with open(outfile, 'w') as w, open(statfile, 'r') as s:
        head = s.readline().strip().split('\t')
        stat_taxon = []
        while True:
            line = s.readline()
            if not line:
                break
            stat_taxon.append(line.split('\t')[0])
        w.write('\teffectsize\tlowerCI\tupperCI\n')
        for name in stat_taxon:
            i = taxon_list.index(name)
            taxon = taxon_list[i]
            R1 = s1_list[i] / totalSeq1
            R2 = s2_list[i] / totalSeq2
            diff = R1 - R2
            stdErr = math.sqrt((R1*(1-R1)) / totalSeq1 + (R2*(1-R2)) / totalSeq2) 
            offset = zScore(coverage) * stdErr
            lowerCI = diff - offset
            upperCI = diff + offset
            diff *= 100
            lowerCI *= 100
            upperCI *= 100
            w.write('%s\t%s\t%s\t%s\n' % (taxon, '%0.4g' % diff, '%0.4g' % lowerCI, '%0.4g' % upperCI))


def NewcombeWilsonFindRoots(seq, totalSeq, z):
    '''
    Find roots required by Newcombe-Wilson CI method
    '''
    value = 0.0
    stepSize = 1.0 / max(totalSeq,1000)
    steps = int(1.0 / stepSize)
    prevP = z*math.sqrt(value*(1.0-value) / totalSeq) - abs(value - float(seq) / totalSeq)
    prevValue = value
    roots = []
    for dummy in xrange(0,steps):
        p = z*math.sqrt(value*(1.0-value) / totalSeq) - abs(value - float(seq) / totalSeq)
        if p*prevP < 0 or (p == 0 and value == 0) or (p == 0 and value == 1.0):
            # we have found a root since there is a sign change
            if abs(p)+abs(prevP) != 0:
                root = prevValue + stepSize*(1.0 - (abs(p)/(abs(p)+abs(prevP))))
            else:
                root = prevValue
            roots.append(root)

            if len(roots) == 2:
                break

        prevP = p
        prevValue = value
        value += stepSize
    # check if we have a double root
    if len(roots) == 1:
        roots.append(roots[0])
    return roots


def NewcombeWilson(otufile, statfile, sample1, sample2, coverage, outfile):
    (taxon_list, s1_list, s2_list, totalSeq1, totalSeq2) = otu_info(otufile, sample1, sample2)
    with open(outfile, 'w') as w, open(statfile, 'r') as s:
        head = s.readline().strip().split('\t')
        stat_taxon = []
        while True:
            line = s.readline()
            if not line:
                break
            stat_taxon.append(line.split('\t')[0])
        w.write('\teffectsize\tlowerCI\tupperCI\n')
        z = zScore(coverage)
        for name in stat_taxon:
            i = taxon_list.index(name)
            taxon = taxon_list[i]
            roots1 = NewcombeWilsonFindRoots(s1_list[i], totalSeq1, z)
            roots2 = NewcombeWilsonFindRoots(s2_list[i], totalSeq2, z)
            diff = s1_list[i]/totalSeq1 - s2_list[i]/totalSeq2
            lowerCI = z*math.sqrt(roots1[0]*(1-roots1[0])/totalSeq1 + roots2[1]*(1-roots2[1])/totalSeq2)
            upperCI = z*math.sqrt(roots1[1]*(1-roots1[1])/totalSeq1 + roots2[0]*(1-roots2[0])/totalSeq2)
            diff *= 100
            lowerCI *= 100
            upperCI *= 100
            w.write('%s\t%s\t%s\t%s\n' % (taxon, '%0.4g' % diff, '%0.4g' % (diff-lowerCI), '%0.4g' % (diff+upperCI)))



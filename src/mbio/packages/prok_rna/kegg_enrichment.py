# -*- coding:utf-8 -*-
from __future__ import print_function
from collections import defaultdict
from scipy.stats import hypergeom
import numpy as np
import argparse
import glob
import re
from textwrap import wrap
import matplotlib as mpl
import matplotlib.pyplot as plt
#plt.style.use('ggplot')
__author__ = "gdq"

def parse_gene_annot(map_file, header=True):
    """
    parse gene_classification file such as:
    ----------------------
    GOLocus Tag Accession
    PA0001  GO:0006260
    PA0001  GO:0003677
    ...
    PA0009  GO:0003688
    **Note, also support line likes: "PA0001  GO:0006260;GO:0003677"
    ----------------------
    """
    gene_class = defaultdict(set)
    class_gene = defaultdict(set)
    d_r = defaultdict(set)
    f = open(map_file)
    if header:
        head = f.readline()
    for line in f:
        if line.strip():
            if len(line.strip().split())<=1:
                continue
            a, b = line.strip().split()
            if ";" in b:
                cls = b.split(";")
                gene_class[a].update(cls)
                for each_class in cls:
                    class_gene[each_class].add(a)
            else:
                gene_class[a].add(b)
                class_gene[b].add(a)
    f.close()
    return class_gene, gene_class


def read_diff_genes(deg_file):
    """ get diff gene list"""
    with open(deg_file) as f:
        deg_list = [x.strip().split()[0:2] for x in f if x.strip()]
    return dict(deg_list)

def parse_br08901(brite):
    """ The brite file can be download from http://www.kegg.jp/kegg-bin/get_htext?br08901 """
    f = open(brite)
    cls_dict = dict()
    for line in f:
        if not line:
            continue
        if line.startswith('A<'):
            A = line[4:-5]
        if line.startswith('B'):
            B = line[1:].strip()
        if line.startswith('C'):
            C = line.split()[1]
            C_detail = ' '.join(line.split()[2:])
            cls_dict['path:ko'+C]=(C_detail, B, A)
            cls_dict['ko'+C]=(C_detail, B, A)
            cls_dict['map'+C]=(C_detail, B, A)
            cls_dict[C]=(C_detail, B, A)
            cls_dict['path:map'+C]=(C_detail, B, A)
    f.close()
    return cls_dict

def enrichment_bar(xdata, ydata, category, stat_value, stat_cutoff=[0.001, 0.01, 0.05], cmap='PiYG', fig_name='enrichment_bar', fig_size=(10, 6),
    wrap_length=75, dpi=300, label_size=5, rotation=45, ylabel='Enrichment Ratio', stat_type='pvalue',
    xlabel=None, title=None, category_detail_dict=None):
    """
    Plot bar graph
    :param xdata: x axis tick label information
    :param ydata: bar height information
    :param category: information for classify bar
    :param stat_value: pvalue or FDR
    :param stat_cutoff:
    :param cmap:
    :param fig_name:
    :param fig_size:
    :param wrap_length:
    :param dpi:
    :param label_size:
    :param rotation:
    :param ylabel:
    :param stat_type: pvalue or FDR
    :param xlabel:
    :param title:
    :return:
    """
    # -log10 stat_value  and normalization
    #stat_value = -np.log10(stat_value)
    stat_value = np.array(stat_value)
    bar_number = len(stat_value)
    # print(bar_number)
    max_value = np.sort(stat_value)[int(bar_number*0.98)]
    #mcn = mpl.colors.Normalize(vmax=round(max_value),clip=False)
    mcn = mpl.colors.Normalize(vmax=max_value,clip=False)
    # combine data into np.array to facilitate sorting
    data = np.array([xdata, ydata, category, stat_value])
    data = data.transpose()
    data = sorted(data, key=lambda X: (X[2], X[3]))
    data = np.array(data)
    fig, ax = plt.subplots(figsize=fig_size)
    # plot each bar
    cmap = plt.get_cmap(cmap)
    stat_value = np.array(data[:,-1], dtype=float)
    color_values = mcn(stat_value)
    bar_colors = cmap(color_values)
    bar_heights = np.array(data[:,1], dtype=float)
    bar_locs = range(1, bar_number+1)
    # print(bar_locs)
    # print(bar_heights)
    ax.bar(bar_locs, bar_heights, color=bar_colors, width=0.8, align='center', linewidth=0.2)
    # turn ticks off
    plt.tick_params(top='off', right='off', direction='out', labelsize='small')
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # change background color
    ax.grid(axis='y', linestyle='dotted', alpha=0.3)
    # draw '*' on significant bar
    min_y, max_y = ax.get_ylim()
    ytick_num = len(ax.get_yticks())
    unit = float(max_y-min_y)/(ytick_num-1)
    for i in range(len(bar_locs)):
        cutoff_list = stat_cutoff
        for n, cutoff in enumerate(cutoff_list):
            if stat_value[i] <= cutoff:
                sig_pos = bar_heights[i]
                for p in range(0, len(cutoff_list)-n):
                    sig_pos += unit*0.13
                    plt.text(bar_locs[i]-0.22, sig_pos, '*', fontsize=7, color='gray')
                    break
          
    # add text of category detail
    if category_detail_dict:
        pos = max_y*0.96 + float(max_y-min_y)*0.3/len(category_detail_dict.keys())
        for c in category_detail_dict:
            pos = pos - float(max_y-min_y)*0.3/len(category_detail_dict.keys())
            plt.text(0.5, pos, c + ': '+category_detail_dict[c], color='k',fontsize=8,alpha=0.5)
    # plot color bar
    ax2,kw = mpl.colorbar.make_axes_gridspec(ax, pad=0.02)
    cbar = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, spacing='proportional', norm=mcn)
    #cbar.ax.set_xlabel('$-log10('+stat_type+')$', fontsize=8)
    cbar.ax.set_xlabel(stat_type, fontsize=8)
    cbar.ax.tick_params(labelsize='small')
    #cbar.set_ticks([mcn(x)/max_value for x in range(1,int(max(color_values))+1)]+[1])
    #cbar.set_ticklabels(list(range(1,int(max(color_values))+1))+[round(max_value,2)])
    # set xtick label, description
    xlabels = data[:, 0]
    label_lens = [len(x) for x in xlabels]
    label_lens.sort()
    # tricks for xtick labels when they are too long
    if label_lens[-1] > wrap_length :
        max_lens = label_lens[int(len(xlabels) * 0.85)]
        xlabels = ['\n'.join(wrap(x, max_lens)) for x in xlabels]
    ax.set_xticks(bar_locs)
    ax.set_xticklabels(xlabels, rotation= rotation, fontsize=label_size, ha='right')
    ax.set_xlim(left=0, right=bar_number+1, auto=True)
    # set spines of ax

    # set ylabel and fig name
    ax.set_ylabel(ylabel, fontsize=9, weight='bold')
    if xlabel: ax.set_xlabel(xlabel)
    if title: fig.suptitle(title)
    # save fig
    fig.tight_layout()
    fig.savefig(fig_name, dpi=dpi, bbox_inches='tight')
    plt.close('all')


def prepare_hypergeom_data(class_gene_dict, gene_class_dict, deg_dict, total_gene_number, gene2ec_dict, gene2k_dict, path_annot_dict):
    pop_number = total_gene_number
    study_number = len(deg_dict)
    # get all DE gene associated classification, named considered_classes 
    considered_classes = set()
    for gene in deg_dict.keys():
        if gene not in gene_class_dict:
            # print(gene, 'not annotated')
            pass
        else:
            considered_classes.update(gene_class_dict[gene])
    # print(considered_classes)

    for each_class in considered_classes:
        associated_genes = class_gene_dict[each_class]
        pop_hitnumber = len(associated_genes)
        associated_diff_genes = associated_genes.intersection(set(deg_dict.keys()))
        associated_diff_info = list()
        ks, colors = list(), list()
        # annote associated genes
        for each_gene in associated_diff_genes:
            regulate = deg_dict[each_gene]
            # gene -> k id
            if each_gene in gene2k_dict:
                k_id = gene2k_dict[each_gene]
                # color = 'red' if deg_dict[each_gene].lower() == 'up' else 'green'
                if regulate.lower() == 'up':
                    color = 'red'
                elif regulate.lower() == 'up':
                    color = 'green'
                else:
                    color = 'red' # edited by shijin
                for each_kid in k_id:
                    if each_kid in ks:  # becasuse many genes -> one Kxxxxx
                        if colors[ks.index(each_kid)] == color:
                            pass
                        else:
                            if colors[ks.index(each_kid)] == 'purple':
                                colors[ks.index(each_kid)] = color
                            elif color == 'purple':
                                pass
                            else:
                                colors[ks.index(each_kid)] = 'yellow'
                    else:
                        ks.append(each_kid)
                        colors.append(color)
            else:
                k_id = ''
            # gene -> EC number
            if each_gene in gene2ec_dict:
                enzymes = gene2ec_dict[each_gene]
            else:
                enzymes = ''

            if k_id or enzymes:
                associated_diff_info.append(each_gene + '|'+ regulate + '|' + '|'.join(k_id) + '|' + '|'.join(enzymes))  
            else:
                associated_diff_info.append(each_gene + '|'+ regulate)
            
        if ks:
            mark = ''
            for k, c in zip(ks, colors):
                mark = mark + '/' + k + '%09' + c
            if "path:" in each_class:
                link = "http://www.genome.jp/kegg-bin/show_pathway?{p}{h}".format(p=each_class[5:], h=mark)
            else:
                link = "http://www.genome.jp/kegg-bin/show_pathway?{p}{h}".format(p=each_class, h=mark)
        else:
            link ='None'

        each_class2 = re.findall(r'\d+',each_class)[0]

        path_name, typeII, typeI = path_annot_dict[each_class2]
        associated_diff_info = '|'.join([x.split("|")[0] for x in associated_diff_info])
        study_hitnumber = len(associated_diff_genes)
        yield study_hitnumber, pop_number, pop_hitnumber, study_number, each_class, associated_diff_info, link, path_name, typeII, typeI


def multtest_correct(p_values, methods=3):
    """
    1. Bonferroni
    2. Bonferroni Step-down(Holm)
    3. Benjamini and Hochberg False Discovery Rate
    4. FDR Benjamini-Yekutieli
    :param pvalue_list:
    :param methods:
    :return: np.array
    """
    def fdrcorrection(pvals, alpha=0.05, method='indep', is_sorted=False):
        '''pvalue correction for false discovery rate
        This covers Benjamini/Hochberg for independent or positively correlated and
        Benjamini/Yekutieli for general or negatively correlated tests. Both are
        available in the function multipletests, as method=`fdr_bh`, resp. `fdr_by`.
        Parameters
        ----------
        pvals : array_like
            set of p-values of the individual tests.
        alpha : float
            error rate
        method : {'indep', 'negcorr')
        Returns
        -------
        rejected : array, bool
            True if a hypothesis is rejected, False if not
        pvalue-corrected : array
            pvalues adjusted for multiple hypothesis testing to limit FDR
        Notes
        -----
        If there is prior information on the fraction of true hypothesis, then alpha
        should be set to alpha * m/m_0 where m is the number of tests,
        given by the p-values, and m_0 is an estimate of the true hypothesis.
        (see Benjamini, Krieger and Yekuteli)
        The two-step method of Benjamini, Krieger and Yekutiel that estimates the number
        of false hypotheses will be available (soon).
        Method names can be abbreviated to first letter, 'i' or 'p' for fdr_bh and 'n' for
        fdr_by.
        '''
        def _ecdf(x):
            '''
            no frills empirical cdf used in fdrcorrection
            '''
            nobs = len(x)
            return np.arange(1, nobs+1)/float(nobs)
            
        pvals = np.asarray(pvals)
        if not is_sorted:
            pvals_sortind = np.argsort(pvals)
            pvals_sorted = np.take(pvals, pvals_sortind)
        else:
            pvals_sorted = pvals  # alias
    
        if method in ['i', 'indep', 'p', 'poscorr']:
            ecdffactor = _ecdf(pvals_sorted)
        elif method in ['n', 'negcorr']:
            cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))   #corrected this
            ecdffactor = _ecdf(pvals_sorted) / cm
        else:
            raise ValueError('only indep and negcorr implemented')
        reject = pvals_sorted <= ecdffactor*alpha
        if reject.any():
            rejectmax = max(np.nonzero(reject)[0])
            reject[:rejectmax] = True
    
        pvals_corrected_raw = pvals_sorted / ecdffactor
        pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
        del pvals_corrected_raw
        pvals_corrected[pvals_corrected>1] = 1
        if not is_sorted:
            pvals_corrected_ = np.empty_like(pvals_corrected)
            pvals_corrected_[pvals_sortind] = pvals_corrected
            del pvals_corrected
            reject_ = np.empty_like(reject)
            reject_[pvals_sortind] = reject
            return reject_, pvals_corrected_
        else:
            return reject, pvals_corrected
    
    pvalue_list = list(p_values)
    n = len(pvalue_list)
    if methods == 1:
        fdr = [eachP*n for eachP in pvalue_list]
    elif methods == 2:
        sorted_pvalues = sorted(pvalue_list)
        fdr = [eachP*(n - sorted_pvalues.index(eachP)) for eachP in pvalue_list]
    elif methods == 3:
        sorted_pvalues = sorted(pvalue_list)
        fdr = [eachP*n/(sorted_pvalues.index(eachP)+1) for eachP in pvalue_list]
    elif methods == 4:
        _, fdr = fdrcorrection(pvalue_list, alpha=0.05, method='negcorr', is_sorted=False)
    fdr = np.array(fdr)
    fdr[fdr > 1] = 1.
    return fdr


def hypergeom_test(data, sort_fdr=True, correct_method=3):
    p_value_list = []
    ratio_in_study_list = []
    ratio_in_pop_list = []
    classes = []
    hit_genes = []
    hit_links = []
    study_number_list = []
    path_names, typeIIs, typeIs = list(), list(), list()
    for study_hitnumber, pop_number, pop_hitnumber, study_number, each_class, associated_diff_info, link, path_name, typeII, typeI in data:
        study_number_list.append(study_hitnumber)
        p_value = 1 - hypergeom.cdf(study_hitnumber-1,pop_number, pop_hitnumber, study_number)
        ratio_in_study = str(study_hitnumber)+'/'+str(study_number)
        ratio_in_pop = str(pop_hitnumber)+'/'+str(pop_number)
        p_value_list.append(p_value)
        ratio_in_study_list.append(ratio_in_study)
        ratio_in_pop_list.append(ratio_in_pop)
        classes.append(each_class)
        hit_genes.append(associated_diff_info)
        hit_links.append(link)
        path_names.append(path_name)
        typeIIs.append(typeII)
        typeIs.append(typeI)

    q_value_list = multtest_correct(p_value_list, methods=correct_method)
    number = len(q_value_list)
    #print(number)
    databases = ['KEGG PATHWAY']*number
    result = zip(study_number_list, path_names, databases, classes,ratio_in_study_list,ratio_in_pop_list, p_value_list, q_value_list, hit_genes, hit_links, typeIIs, typeIs)
    if sort_fdr:
        sorted_result = sorted(result, key=lambda x: (x[7],x[6]))
    else:
        sorted_result = sorted(result, key=lambda x: (x[6],x[7]))
    #num = sum([1 for x in p_value_list if x < 1])
    #return sorted_result[0:num]
    return sorted_result

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-deg', required=True, help='file with two columns: gene\tup/down')
    parser.add_argument('-g2p', required=True, help='file with two columns: gene\tpath:konumber')
    parser.add_argument('-g2k', required=True, help='file with two columns: gene\tgene_Knumber')
    parser.add_argument('-bgn', required=True, type=int, help='int, total background gene number')
    parser.add_argument('-k2e', 
        default='/mnt/ilustre/users/deqing.gu/Databases/Enrichment/K2enzyme.tab',
        help='K and enzyme mapping file; Default file is in ~deqing.gu/Databases/Enrichment/')
    parser.add_argument('-brite', 
        default='/mnt/ilustre/users/deqing.gu/Databases/Enrichment/br08901.txt',
        help='The brite file can be download from http://www.kegg.jp/kegg-bin/get_htext?br08901; Default file is in ~deqing.gu/Databases/Enrichment/')
    parser.add_argument('--FDR', default=False, action='store_true', help='if used, FDR will be used for plotting')
    parser.add_argument('-dn', metavar='draw_number', default=0, type=int, help='plot with the top "dn" terms. Default: draw all terms whose pvalue/fdr<1')
    parser.add_argument('--rm_HD_DD', default=False, action='store_true', help='Do not draw term classified as HumanDisease and DrugDevelopement')
    parser.add_argument('-correct', metavar='correct_method', default=3, type=int, help='multi-test method')
    args = parser.parse_args()

    g2p_file = args.g2p
    deg_file = args.deg
    gene_number = args.bgn
    g2k_file = args.g2k
    k2e_file = args.k2e
    brite_file = args.brite
    FDR = args.FDR
    draw_number = args.dn
    rm_HD_DD = args.rm_HD_DD
    correct_method=args.correct

    # calculating result
    p_gene, gene_p = parse_gene_annot(g2p_file)
    k_gene, gene_k = parse_gene_annot(g2k_file)
    enzyme_k, k_enzyme = parse_gene_annot(k2e_file) 
    path_annot_dict = parse_br08901(brite_file)

    g2e_dict = defaultdict(set)
    for gene, ks in gene_k.items():
        for k in ks:
            for e in k_enzyme[k]:
                g2e_dict[gene].add(e)
    
    deg_files = glob.glob(deg_file)
    # print(deg_files)
    for deg_file in deg_files:
        deg_list = read_diff_genes(deg_file)
        data = prepare_hypergeom_data(p_gene, gene_p, deg_list, gene_number, g2e_dict, gene_k, path_annot_dict)
        result = hypergeom_test(data, sort_fdr=FDR, correct_method=correct_method)
        f = open(deg_file + '.kegg_enrichment.xls', 'w')
        # print(deg_file + '.kegg_enrichment.xls')
        f.write('#Study_num\tTerm\tDatabase\tID\tRatio_in_study\tRatio_in_pop\tP-Value\tCorrected P-Value\tGenes\tHyperlink\ttypeII\ttypeI\n')
        for each in result:
            f.write('\t'.join([str(x) for x in each])+'\n')
        f.close()
    """
        ## plot bar graph
        # path_names(0), databases(1), classes(2),ratio_in_study_list(3),ratio_in_pop_list(4), p_value_list(5), q_value_list(6), 
        # hit_genes(7), hit_links(8), typeIIs(9), typeIs(10)
    tmp_result = zip(*result)
    tmp_result = np.array(tmp_result)
    tmp_result = tmp_result.transpose()
    if rm_HD_DD:
        tmp_result = tmp_result[tmp_result[:,10] != 'Human Diseases']
        tmp_result = tmp_result[tmp_result[:,10] != 'Drug Development']
    if draw_number == 0:
            if not FDR:
            num = sum([1 for x in tmp_result[:,5] if float(x) < 1])
        else:
        num = sum([1 for x in tmp_result[:,6] if float(x) < 1])
        filtered_result = tmp_result[0:num, :]
    else:
        filtered_result = tmp_result[0:draw_number, :]

        result = filtered_result
        # print(result)
        # FDR or pvaue to be used
        if not FDR:
            stat_value = np.array(result[:, 5], dtype=float)
            stat_type = 'pvalue'
        else:
            stat_value = np.array(result[:, 6], dtype=float)
            stat_type = 'FDR'
        #print(stat_value)
        # get x coor
        xdata = list(result[:, 0])
        # get y coor
        study_hit = np.array([float(x.split('/')[0]) for x in result[:, 3]])
    study_number = np.array([float(x.split('/')[1]) for x in result[:, 3]])
        pop_hit = np.array([float(x.split('/')[0]) for x in result[:, 4]])
    pop_number = np.array([float(x.split('/')[1]) for x in result[:, 4]])
    #for odd ratio: 
    #ydata = (study_hit/(study_number-study_hit))/((pop_hit-study_hit)/(pop_number-study_number-pop_hit+study_hit))
        #ydata = (study_hit/study_number)/(pop_hit/pop_number)
    
    # if you prefer the bar height is study_hitnumber/pop_hitnumber, just un comment the following 3 lines
        study_hitnumber = [float(x.split('/')[0]) for x in result[:, 3]]
        pop_hitnumber = [float(x.split('/')[0]) for x in result[:, 4]]
        ydata = np.array(study_hitnumber)/np.array(pop_hitnumber)
        category = result[:, 10]
        category[category=='Metabolism'] = 'M'
        category[category=='Genetic Information Processing'] = 'GIP'
        category[category=='Environmental Information Processing'] = 'EIP'
        category[category=='Cellular Processes'] = 'CP'
        category[category=='Organismal Systems'] = 'OS'
        category[category=='Human Diseases'] = 'HD'
        category[category=='Drug Development'] = 'DD'
    
        category_dict = dict()
        category_dict['M'] = 'Metabolism'
        category_dict['GIP'] = 'Genetic Information Processing'
        category_dict['EIP'] = 'Environmental Information Processing'
        category_dict['CP'] = 'Cellular Processes'
        category_dict['OS'] = 'Organismal Systems'
        category_dict['DD'] = 'Drug Development'
        category_dict['HD'] = 'Human Diseases'
    for each in category_dict.keys():
        if not each in list(category):
        category_dict.pop(each)

        new_xdata = []
        for i, c in enumerate(list(category)):
            #print(i, c)
            if c == 'M':
                new_xdata.append(xdata[i]+'<M>')
            elif c == 'GIP':
                new_xdata.append(xdata[i]+'<GIP>')
            elif c == 'EIP':
                new_xdata.append(xdata[i]+'<EIP>')
            elif c == 'CP':
                new_xdata.append(xdata[i] + '<CP>')
            elif c == 'OS':
                new_xdata.append(xdata[i] + '<OS>')
            elif c == 'HD':
                new_xdata.append(xdata[i] + '<HD>')
            elif c == 'DD':
                new_xdata.append(xdata[i] + '<DD>')
            else:
                new_xdata.append(xdata[i])
                print(xdata[i], " is not classfied?")
        
    # draw 
        enrichment_bar(
            new_xdata, 
            ydata, 
            category, 
            stat_value, 
            category_detail_dict=category_dict, # classification detail
            cmap='PiYG', # color of bar
            fig_name=deg_file+'.kegg_enrichment.pdf', # out put PDF
            ylabel='Enrichment Ratio', # bar height
            stat_type=stat_type, # pvalue or FDR
            stat_cutoff=[0.001, 0.01, 0.05], 
            fig_size=(10, 6),
            wrap_length=75, 
            dpi=300, # resolution for figure type such as png or tiff
            label_size=5, # tick label size
            rotation=45, # tick label direction
            )
    """

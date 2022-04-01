# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20180524

from mbio.packages.metabolome.scripts.metab_cluster_toolbox import hcluster, r_hcluster, kmeans
import fastcluster as hclust
import sys
import argparse
import os, glob
import pandas as pd
import shutil
import scipy
from scipy import stats
from sklearn import decomposition, preprocessing


class CorrClust(object):
    '''
    计算相关性和pvalue，并聚类
    '''
    def __init__(self):
        self.sel_length_table = ""

    def corr_and_pvalue(self, profile, cor_method, outDir, trans=False):
        """
        计算相关性和pvalue
        :param profile: DataFrame of expression matrix
        :param cor_method: Choices: [‘pearson’, ‘kendall’, ‘spearman’]
        """
        profile = pd.read_table(profile, sep='\t', header=0, index_col=0)
        if trans:
            profile = profile.T
        # table_matrix = table.T.corr(coefficient)
        length = len(profile.columns)
        corr_table = []
        p_table = []
        samples = []
        if length < 2:
            raise Exception('行数小于2，无法计算相关性!')
        for i in range(0, length):
            correlations = []
            pvalues = []
            sample = profile.columns[i]
            samples.append(sample)
            for j in range(0, length):
                (correlation, pvalue) = self.coeff_cal(cor_method, profile.iloc[:, i], profile.iloc[:, j])
                correlations.append(correlation)
                pvalues.append(pvalue)
            corr_frame = pd.DataFrame(correlations, columns=[sample])
            p_frame = pd.DataFrame(pvalues, columns=[sample])
            corr_table.append(corr_frame)
            p_table.append(p_frame)
        final1 = pd.concat(corr_table, axis=1)
        final2 = pd.concat(p_table, axis=1)
        final1.index = samples
        final2.index = samples
        final1.to_csv(outDir + "/ori_corr.xls", sep="\t")
        final2.to_csv(outDir + "/ori_pvalue.xls", sep="\t")
        return final1, final2

    def coeff_cal(self, coefficient, x, y):
        if coefficient == "spearman":
            result = scipy.stats.spearmanr(x, y)
        elif coefficient == "pearson":
            result = scipy.stats.pearsonr(x, y)
        elif coefficient == "kendall":
            result = scipy.stats.kendalltau(x, y)
        correlation = round(result[0], 5)
        pvalue = result[1]
        return (correlation, pvalue)

    def cor_cluster(self, type, exp_pd, use_c, use_r, c_clu_type, r_clu_type, c_clu_method, r_clu_method,
                    c_dist_method, r_dist_method, outDir, n_clusters=0, sn_clusters=0):
        """
        聚类
        :param exp_pd: DataFrame of expression matrix
        :param cluster_method: Choices: ['single', 'average', 'weighted', 'centroid', 'complete', 'median', 'ward']
        :param dist_method: ['braycurtis', 'canberra', 'chebyshev', 'cityblock',
        'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski',
        'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao',
        'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule', ]
        """
        if len(exp_pd.columns) < 2:
            use_c = None
        if len(exp_pd) < 2:
            use_r = None
        use_rscr = False
        use_sub = False
        result = []
        if use_rscr:
            cluster = r_hcluster
        else:
            cluster = hcluster
        if use_c:
            # for sample
            if type == "corr":
                prefix = "corr."
            elif type == "clu":
                prefix = "col."
            if c_clu_type == 'hierarchy':
                s_cluster_result = cluster(exp_pd, transpose=True, n_clusters=sn_clusters, method=c_clu_method,
                                           metric=c_dist_method, output=outDir, prefix=prefix, use_sub=use_sub)
            else:
                s_cluster_result = kmeans(exp_pd, n_clusters=sn_clusters, transpose=True, output=outDir,
                                          prefix=prefix)
            result.append(s_cluster_result)
        if use_r:
            # for metabolome or gene
            if type == "corr":
                prefix = "corr."
            elif type == "clu":
                prefix = "row."
            if r_clu_type == 'hierarchy':
                use_sub = True  #20190618
                r_cluster_result = cluster(exp_pd, transpose=False, n_clusters=n_clusters, method=r_clu_method,
                                           metric=r_dist_method, output=outDir, prefix=prefix, use_sub=use_sub)
            else:
                r_cluster_result = kmeans(exp_pd, n_clusters=n_clusters, transpose=False, output=outDir, prefix=prefix)
            result.append(r_cluster_result)
        return result

    def order_sample(self, exp_table, outfile=None, coltree=None, rowtree=None):
        """
        根据tree文件重排序表达矩阵顺序
        """
        table = exp_table
        if coltree:
            with open(coltree, "r") as f:
                tree = f.next()
                col_tree_order = f.next().strip().split(";")
            table = table[col_tree_order]
        if rowtree and rowtree != coltree:
            with open(rowtree, "r") as f:
                tree = f.next()
                row_tree_order = f.next().strip().split(";")
        elif rowtree and rowtree == coltree:
            row_tree_order = col_tree_order
        if rowtree:
            table = table.loc[row_tree_order]
        if outfile:
            table.to_csv(outfile, sep="\t", index=True)
        return table

    def order_kmeans(self, input_table, outfile=None, outDir=None, corr=False, cluster=False, cluster_col=False):
        """
        根据kmeans的subcluster结果重排序
        """
        all_subs = []
        table = input_table
        if corr:
            subclusters = glob.glob(outDir + '/corr.subcluster*')
            sub_number = len(subclusters)
            for i in range(0, sub_number):
                name = "/*" + "corr.subcluster_" + str(i) + "*"
                subfile = glob.glob(outDir + name)
                subtable = pd.read_table(subfile[0], sep='\t', header=0, index_col=0)
                metabs = subtable.columns.tolist()
                all_subs = all_subs + metabs
            table = table[all_subs]
            table = table.loc[all_subs]
        cluster_subs = []
        if cluster:
            subclusters = glob.glob(outDir + '/*row.subcluster*')
            sub_number = len(subclusters)
            for i in range(0, sub_number):
                name = "/*" + "row.subcluster_" + str(i) + "*"
                subfile = glob.glob(outDir + name)
                subtable = pd.read_table(subfile[0], sep='\t', header=0, index_col=0)
                metabs = subtable.index.tolist()
                cluster_subs = cluster_subs + metabs
            table = table.loc[cluster_subs]
            #print table.head()
        if cluster_col:
            col_subs = []
            subclusters = glob.glob(outDir + '/col.subcluster*')
            sub_number = len(subclusters)
            for i in range(0, sub_number):
                name = "/*" + "col.subcluster_" + str(i) + "*"
                subfile = glob.glob(outDir + name)
                subtable = pd.read_table(subfile[0], sep='\t', header=0, index_col=0)
                metabs = subtable.columns.tolist()
                col_subs = col_subs + metabs
            table = table[col_subs]
        if outfile:
            table.to_csv(outfile, sep="\t")
        return table

    def remove_file(self, outDir):
        subclusters = glob.glob(outDir + '/*subcluster*')
        for each in subclusters:
            os.remove(each)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-exp', type=str, metavar="exp_matrix_file", required=True,
                        help="expression matrix file")
    #parser.add_argument('-group', type=str, metavar="group_file", default=None, help="file with two col: sample\tgroup")
    parser.add_argument('-out', type=str, default=None, help="default is local dir. Output directory name")
    parser.add_argument('--T', action='store_true', help="need transform exp_file")
    parser.add_argument('-n_clusters', type=int, default=10, metavar="cluster_num", help="expected_cluster_number")
    parser.add_argument('-sn_clusters', type=int, default=0, metavar="cluster_num",
                        help="sample expected_cluster_number")
    parser.add_argument('--nsc', action='store_true', help="'no-sample-cluster', do not perform sample cluster")
    parser.add_argument('--ngc', action='store_true', help="'no-gene-cluster', do not perform gene cluster")
    parser.add_argument('-sct', metavar='sample-cluster-type', default='hierarchy', help="hierarchy or kmeans cluster")
    parser.add_argument('-gct', metavar='gene-cluster-type', default='hierarchy', help="hierarchy or kmeans cluster")
    parser.add_argument('-scm', metavar='sample-cluster-method', default='complete', help="hclust method")
    parser.add_argument('-gcm', metavar='gene-cluster-method', default='complete', help="hclust method")
    parser.add_argument('-scd', metavar='sample-cluster-distance', default='correlation', help="hclust distance metric")
    parser.add_argument('-gcd', metavar='gene-cluster-distance', default='euclidean', help="hclust distance metric")
    parser.add_argument('--corr', action='store_true', help="'sample_correlation', do sample correlation")
    parser.add_argument('-corr_method', default='pearson', help="pearson, kendall, spearman are supported")
    args = parser.parse_args()
    exp_pd = args.exp
    outDir = args.out
    s_cluster, g_cluster = True, True
    trans = False
    if args.nsc:
        s_cluster = False
    if args.ngc:
        g_cluster = False
    if args.sn_clusters == 0:
        sn_clusters = args.n_clusters
    else:
        sn_clusters = args.sn_clusters
    CorrClust = CorrClust()
    CorrClust.remove_file(outDir)
    if args.corr:
        if args.T:
            trans = True
        corr_table, pvalue_table = CorrClust.corr_and_pvalue(exp_pd, args.corr_method, outDir, trans=trans)
        run = CorrClust.cor_cluster("corr", corr_table, s_cluster, g_cluster, c_clu_type=args.sct,
                                    r_clu_type=args.gct, c_clu_method=args.scm, r_clu_method=args.gcm,
                                    c_dist_method=args.scd, r_dist_method=args.gcd, outDir=outDir,
                                    n_clusters=args.n_clusters, sn_clusters=sn_clusters)
        if not args.nsc:
            if args.sct == "hierarchy" and len(corr_table.columns) > 1:
                tree_file = outDir + "/corr.cluster_tree.txt"
                CorrClust.order_sample(corr_table, outDir + "/corr.xls", coltree=tree_file, rowtree=tree_file)
                CorrClust.order_sample(pvalue_table, outDir + "/pvalue.xls", coltree=tree_file, rowtree=tree_file)
            elif args.sct == "kmeans" and len(corr_table) > 1:
                CorrClust.order_kmeans(corr_table, outDir + "/corr.xls", outDir, corr=True)
                CorrClust.order_kmeans(pvalue_table, outDir + "/pvalue.xls", outDir, corr=True)
            else:
                CorrClust.order_sample(corr_table, outDir + "/corr.xls")
                CorrClust.order_sample(pvalue_table, outDir + "/pvalue.xls")
        else:
            CorrClust.order_sample(corr_table, outDir + "/corr.xls")
            CorrClust.order_sample(pvalue_table, outDir + "/pvalue.xls")
    else:
        exp_table = pd.read_table(exp_pd, sep='\t', header=0, index_col=0)
        run = CorrClust.cor_cluster("clu", exp_table, s_cluster, g_cluster, c_clu_type=args.sct,
                                    r_clu_type=args.gct, c_clu_method=args.scm, r_clu_method=args.gcm,
                                    c_dist_method=args.scd, r_dist_method=args.gcd, outDir=outDir,
                                    n_clusters=args.n_clusters, sn_clusters=sn_clusters)
        col_tree_file = None
        row_tree_file = None
        result_table = exp_table
        outfile = outDir + "/cluster_exp.xls"
        if not args.nsc and args.sct == "hierarchy" and len(exp_table.columns) > 1:
            col_tree_file = outDir + "/col.cluster_tree.txt"
        if not args.ngc and args.gct == "hierarchy" and len(exp_table) > 1:
            row_tree_file = outDir + "/row.cluster_tree.txt"
        if args.sct == "hierarchy" or args.gct == "hierarchy":
            result_table = CorrClust.order_sample(result_table, outfile=None, coltree=col_tree_file,
                                                  rowtree=row_tree_file)
        if args.gct == "kmeans" and len(result_table) > 1:
            result_table = CorrClust.order_kmeans(result_table, outfile=None, outDir=outDir, cluster=True)
        result_table.to_csv(outfile, sep="\t", quoting=3)

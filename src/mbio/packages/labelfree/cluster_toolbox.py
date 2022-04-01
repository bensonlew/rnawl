# coding=utf-8
import os
import re
import math
import pandas as pd
import numpy as np
from collections import OrderedDict
import scipy.cluster.hierarchy as sch
from sklearn.cluster import KMeans
from sklearn import decomposition, preprocessing
import fastcluster as hclust
import matplotlib
matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# from sklearn.preprocessing import StandardScaler


def hcluster(exp_pd, transpose=False, n_clusters=10, method='average',
             metric='correlation', output=None, prefix=''):
    """
    'fastcluster' was used, http://www.danifold.net/fastcluster.html?section=3.
    scipy.cluster.hierarchy.linkage could also do the same thing but slower.
    However, the documentation of scipy.cluster.hierarchy.linkage is pretty good.
    >> https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html

    :param exp_pd: pandas DataFrame, expression matrix
    :param transpose: if to transpose the expression matrix
    :param n_clusters: int, optional, default: 8. The number of clusters to generate.
    :param method: methods for calculating the distance between the newly formed clusters.
        Choices: ['single', 'average', 'weighted', 'centroid', 'complete', 'median', 'ward']
    :param metric: The distance metric to use.
        Choices: ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation',
        'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski',
        'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao',
        'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule', ]
    :param output: output directory of subcluster information
    :param prefix: outfile prefix letters
    :return: tree -> tuple(tree_str, tree_list),
             subcluster -> {0:[s1,s2], 1:[s4,s5], ...}
             z -> cluster result of hclust.linkage
    """
    if transpose:
        exp_pd = exp_pd.transpose()
    if n_clusters > exp_pd.shape[0]:
        print("n_clusters is bigger than sample number!!")
        n_clusters = exp_pd.shape[0]
    try:
        z = hclust.linkage(exp_pd, method=method, metric=metric)
    except FloatingPointError as e:
        print("fastcluster failed as : {}".format(e))
        print('it seems that (at least) one of the vectors you want to cluster is all zeros, '
              'so when it tries to compute the cosine distances to it there is a division by zero,'
              ' hence nan is stored in your distance array, and that leads to your error.')
        print('Anyway, we will remove the special rows for you now.')
        check = exp_pd[exp_pd.sum(axis=1) == 0]
        if check.shape[0] >= 1:
            print('Actually, we detected that some genes have zero expression across all samples.')
            print('such as {} has zero expression across all sample'.format(check.index[0]))
            exp_pd = exp_pd[exp_pd.sum(axis=1) > 0]
        try:
            z = hclust.linkage(exp_pd, method=method, metric=metric)
        except:
            print("enhen? fastcluster failed again, we will try scipy.cluster.hierarchy")
            z = sch.linkage(exp_pd, method=method, metric=metric)
    except:
        print("fastcluster failed, we will try scipy.cluster.hierarchy")
        z = sch.linkage(exp_pd, method=method, metric=metric)
    labels = exp_pd.index
    tree = to_tree(z, labels)
    # tree = R_to_tree(exp_pd, metric, method)
    subcluster = get_subcluster(z, labels, num=n_clusters)
    # write out subcluster
    if output is not None:
        if not os.path.exists(output):
            os.mkdir(output)
    else:
        output = os.getcwd()
    for k, v in subcluster.items():
        out_dir = os.path.join(output, prefix + 'subcluster_{}_{}.xls'.format(k, len(v)))
        sub = exp_pd.loc[v, :]
        if transpose:
            sub = sub.transpose()
        sub.to_csv(out_dir, sep='\t', header=True, index=True)
    # write out tree
    out_dir = os.path.join(output, prefix + "cluster_tree.txt")
    with open(out_dir, 'w') as f:
        f.write(tree[0]+'\n')
        f.write(";".join(tree[1])+'\n')
    # write out cluster result z
    out_dir = os.path.join(output, prefix + "linkage_result")
    pd.DataFrame(z).to_csv(out_dir, sep='\t')
    return tree, subcluster, z


def to_tree(zclust, labels):
    """
    tree, subcluster, z
    :param zclust： hclust.linkage result
    :param labels: leaf label from DataFrame.columns
    :return: a string denotes cluster tree with distance information contained.
    """
    clust_step_detail = dict()
    sample_num = zclust.shape[0] + 1
    new_node = zclust.shape[0]
    for i, j, h, _ in zclust:
        i, j = int(i), int(j)
        new_node += 1
        if i < sample_num:
            i = labels[i]+':{:7f}'.format(h)
        if j < sample_num:
            j = labels[j]+':{:7f}'.format(h)
        if i in clust_step_detail:
            i = clust_step_detail[i]
        if j in clust_step_detail:
            j = clust_step_detail[j]
        clust_step_detail[new_node] = '({i},{j}):{h:7f}'.format(i=i, j=j, h=h)
    else:
        tree = clust_step_detail[new_node]
    # print(clust_step_detail)
    tree_list = labels[sch.leaves_list(zclust)]
    # return '(' + tree + ')', tree_list
    end_ind = tree.rfind('):')
    return tree[:end_ind+1], tree_list


def r_hcluster(exp_pd, transpose=False, n_clusters=10, method='average',
               metric='euclidean', output=None, prefix='', do_scale=False):
    """
        'fastcluster' was used, http://www.danifold.net/fastcluster.html?section=3.
        scipy.cluster.hierarchy.linkage could also do the same thing but slower.
        However, the documentation of scipy.cluster.hierarchy.linkage is pretty good.
        >> https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html

        :param exp_pd: pandas DataFrame, expression matrix
        :param transpose: if to transpose the expression matrix
        :param n_clusters: int, optional, default: 8. The number of clusters to generate.
        :param method: the agglomeration method to be used. This must be (an unambiguous abbreviation of)
                       one of ‘"single"’, ‘"complete"’, ‘"average"’, ‘"mcquitty"’, ‘"ward.D"’, ‘"ward.D2"’,
                       ‘"centroid"’ or ‘"median"’.
        :param metric: the distance measure to be used. This must be one of
                       ‘"euclidean"’, ‘"maximum"’, ‘"manhattan"’, ‘"canberra"’,
                       ‘"binary"’ or ‘"minkowski"’.  Any unambiguous substring can be given.
                       In addition,  correlation method = ("pearson", "kendall", "spearman") is also supported
        :param output: output directory of subcluster information
        :param prefix: outfile prefix letters
        :param do_scale: if do normalization with sklearn.preprocessing.scale
        :return: tree_str , tree_list
                 subcluster -> {0:[s1,s2], 1:[s4,s5], ...}
                 hclust_result -> cluster result of hclust.linkage
        """

    def r_hclust_func(exp_pd, method='complete', metric='euclidean', n_clusters=None):
        from rpy2.robjects.packages import importr
        from rpy2.robjects import r, pandas2ri
        importr("fastcluster")
        importr("ape")
        pandas2ri.activate()
        r_dist = r['dist']
        r_hclust = r['hclust']
        if metric == "cityblock":
            metric = "manhattan"
        if metric == "correlation":
            metric = "pearson"
        if metric == "ward":
            metric = "ward.D"
        if metric in ("pearson", "kendall", "spearman"):
            r_dist = r['as.dist']
            corr_result = 1 - exp_pd.transpose().corr(method=metric)
            corr_ds = pandas2ri.py2ri(corr_result)
            ds = r_dist(corr_ds)
        else:
            exp_df = pandas2ri.py2ri(exp_pd)
            ds = r_dist(exp_df, method=metric)
        hc = r_hclust(ds, method=method)
        # get tree
        r_phylo = r['as.phylo']
        r_write_tree = r['write.tree']
        tree = r_phylo(hc)
        tree_str = str(r_write_tree(tree)).split('"')[1][:-1]
        tree_list = re.findall('[(,]([^(]*?):', tree_str)
        # get subcluster
        if n_clusters:
            r_cut_tree = r['cutree']
            r_subcluster = r_cut_tree(hc, n_clusters)
            subcluster_dict = dict()
            for k, v in r_subcluster.items():
                subcluster_dict.setdefault(v, list())
                subcluster_dict[v].append(k)
        else:
            subcluster_dict = None
        return tree_str, tree_list, subcluster_dict, hc
    if transpose:
        raw_exp = exp_pd.transpose()
    else:
        raw_exp = exp_pd
    if do_scale:
        exp_pd = exp_pd.apply(preprocessing.scale, axis=0)
    if transpose:
        exp_pd = exp_pd.transpose()
    if n_clusters > exp_pd.shape[0]:
        print("n_clusters is bigger than sample number!!")
        n_clusters = exp_pd.shape[0]
    try:
        z = r_hclust_func(exp_pd, method=method, metric=metric, n_clusters=n_clusters)
    except Exception as e:
        print("fastcluster failed as : {}".format(e))
        print('it seems that (at least) one of the vectors you want to cluster is all zeros, '
              'so when it tries to compute the cosine distances to it there is a division by zero,'
              ' hence nan is stored in your distance array, and that leads to your error.')
        print('Anyway, we will remove the special rows for you now.')
        check = exp_pd[exp_pd.sum(axis=1) == 0]
        if check.shape[0] >= 1:
            print('Actually, we detected that some genes have zero expression across all samples.')
            print('such as {} has zero expression across all sample'.format(check.index[0]))
            exp_pd = exp_pd[exp_pd.sum(axis=1) > 0]
        z = r_hclust_func(exp_pd, method=method, metric=metric, n_clusters=n_clusters)

    tree_str, tree_list, subcluster, hclust_result = z
    # write out subcluster
    if output is not None:
        if not os.path.exists(output):
            os.mkdir(output)
    else:
        output = os.getcwd()
    for k, v in subcluster.items():
        out_dir = os.path.join(output, prefix + 'subcluster_{}_{}.xls'.format(k, len(v)))
        # sub = exp_pd.loc[v, :]
        sub = raw_exp.loc[v, :]
        if transpose:
            sub = sub.transpose()
        sub.to_csv(out_dir, sep='\t', header=True, index=True)
    # write out tree
    out_dir = os.path.join(output, prefix + "cluster_tree.txt")
    with open(out_dir, 'w') as f:
        f.write(tree_str + '\n')
        f.write(";".join(tree_list) + '\n')
    return z


def get_subcluster(zclust, labels, num=2):
    """
    get leave ids for each sub cluster
    :param zclust hclust.linkage result
    :param num: the number of sub-clusters specified.
    :param labels: leaf label from DataFrame.columns
    :return: dict with list of samples as element.
    """
    cluster = sch.cut_tree(zclust, num)
    tmp_pd = pd.DataFrame(cluster)
    tmp_pd['label'] = labels
    result = tmp_pd.groupby(0).groups
    subcluster = dict()
    for k in result:
        subcluster[k] = list(labels[result[k]])
    return subcluster


def kmeans(exp_pd, n_clusters=10, transpose=False, output=None, prefix=''):
    """
    http://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html
    :param exp_pd:
    :param n_clusters: int, optional, default: 8
           The number of clusters to form as well as the number of centroids to generate.
    :param transpose: if to transpose the input data
    :param output: output directory of subcluster information
    :param prefix: outfile prefix letters
    :return: dict with list of samples as element.
    """
    if transpose:
        exp_pd = exp_pd.transpose()
    if n_clusters > exp_pd.shape[0]:
        print("n_clusters is bigger than sample number!!")
        n_clusters = exp_pd.shape[0]
    cluster = KMeans(n_clusters=n_clusters, ).fit_predict(exp_pd)
    # k-means clustering and give subcluster
    tmp_pd = pd.DataFrame(cluster)
    labels = exp_pd.index
    tmp_pd['label'] = labels
    result = tmp_pd.groupby(0).groups
    subcluster = dict()
    for k in result:
        subcluster[str(int(k) + 1)] = list(labels[result[k]])
    if output is not None:
        if not os.path.exists(output):
            os.mkdir(output)
    else:
        output = os.getcwd()
    for k, v in subcluster.items():
        out_dir = os.path.join(output, prefix+'subcluster_{}_{}.xls'.format(k, len(v)))
        sub = exp_pd.loc[v, :]
        if transpose:
            sub = sub.transpose()
        sub.to_csv(out_dir, sep='\t', header=True, index=True)
    # write out cluster dict
    cluster_out = os.path.join(output, prefix+'kmeans_cluster.txt')
    with open(cluster_out, 'w') as f:
        for k, v in subcluster.items():
            f.write('{}\t{}\n'.format(k, ";".join(v)))
    return subcluster


def corr(exp_pd, method='pearson', output=None):
    """
    Correlation calculation with pandas.
    :param exp_pd: DataFrame of expression matrix
    :param method: Choices: [‘pearson’, ‘kendall’, ‘spearman’]
                    pearson : standard correlation coefficient
                    kendall : Kendall Tau correlation coefficient
                    spearman : Spearman rank correlation
    :param output: output directory
    :return: correlation matrix
    """
    if output is not None:
        if not os.path.exists(output):
            os.mkdir(output)
    else:
        output = os.getcwd()
    out_dir = os.path.join(output, 'sample_correlation.xls')
    result = exp_pd.corr(method=method)
    result.index.name = 'sample'
    result.to_csv(out_dir, sep='\t', header=True, index=True)
    return result


def pca(all_exp_pd, output=None):
    """
    PCA analysis using sk-learn package of python
    :param all_exp_pd: DataFrame of expression matrix
    :return: DataFrame with samples as rows and component values as columns.
    """
    pca = decomposition.PCA()
    data = all_exp_pd.transpose()
    pca.fit(data)
    _ratio = list(enumerate(pca.explained_variance_ratio_, start=1))
    total_ratio, n_components = 0, 0
    for ind, each in _ratio:
        total_ratio += each
        if total_ratio >= 0.95:
            n_components = ind
            break
    if n_components <= 1:
        n_components = 2
    _ratio = _ratio[:n_components]
    pc_ratio = {'PC'+str(n): r for n, r in _ratio}
    result = pd.DataFrame(pca.transform(data), index=data.index)
    result = result.iloc[:, :n_components]
    result.index.name = 'sample'
    # result.columns = ['PC'+str(n)+'('+'{:.2f}%'.format(r*100)+')' for n, r in _ratio]
    result.columns = ['PC'+str(n) for n in range(1, result.shape[1] + 1)]
    # write out
    if output is not None:
        if not os.path.exists(output):
            os.mkdir(output)
    else:
        output = os.getcwd()
    out_dir = os.path.join(output, 'PCA.xls')
    result.to_csv(out_dir, sep='\t', header=True, index=True)
    out_dir2 = os.path.join(output, 'Explained_variance_ratio.xls')
    with open(out_dir2, 'w') as f:
        for each in sorted(pc_ratio.keys(), key=lambda a:int(a[2:])):
            f.write(str(each) + '\t' + str(pc_ratio[each]) + '\n')
    return result, pc_ratio


def exp_corr_cluster(exp_pd, method='pearson', n_clusters=10, s_method="complete",
                     s_metric="correlation", output=None):
    # for sample cluster
    s_cluster = hcluster(exp_pd, transpose=True, n_clusters=n_clusters, method=s_method,
                         metric=s_metric, output=output, prefix="sample.")
    # for sample corr
    correlation = corr(exp_pd, method=method, output=output)
    return correlation, s_cluster


def process_exp_matrix(exp_matrix, log_base=None, group_dict=None, group_method="mean"):
    if type(exp_matrix) == str or type(exp_matrix) == bytes:
        # all_exp_pd = pd.read_table(exp_matrix, index_col=0, header=0, dtype={0:str})
        # df.index = df.index.map(str)
        all_exp_pd = pd.read_table(exp_matrix, header=0, dtype={0:str})
        all_exp_pd = all_exp_pd.set_index(all_exp_pd.columns[0])
    else:
        print('Input is Not a str ? then exp_matrix is assumed to be a pandas DataFrame Object')
        all_exp_pd = exp_matrix
    all_exp_pd.index.name = 'accession_id'
    #
    if group_dict is not None:
        group_exp = list()
        for g in group_dict:
            if group_method == 'median':
                g_exp = all_exp_pd.loc[:, group_dict[g]].median(axis=1)
            else:
                g_exp = all_exp_pd.loc[:, group_dict[g]].mean(axis=1)
            g_exp.name = g
            group_exp.append(g_exp)
        all_exp_pd = pd.concat(group_exp, axis=1)
    #
    if log_base:
        if log_base == math.e:
            all_exp_pd = np.log(all_exp_pd + 0.0001)
        elif log_base == 2:
            all_exp_pd = np.log2(all_exp_pd + 0.0001)
        elif log_base == 10:
            all_exp_pd = np.log10(all_exp_pd + 0.0001)
        else:
            raise Exception('log base of {} is not supported'.format(log_base))

    #由封一统添加的可以对表达表进行scale的方法——2018-07-12
    if all_exp_pd.shape[1] == 2:   #当以分组均值计算或者一共两个样本时，不除以方差
        scaler = preprocessing.StandardScaler(with_std=False)
    else:
        scaler = preprocessing.StandardScaler()
    scaler.fit(all_exp_pd.T)
    all_exp_pd1 = pd.DataFrame(scaler.transform(all_exp_pd.T)).T
    all_exp_pd1.columns = all_exp_pd.columns
    all_exp_pd1.insert(0, 'accession_id', all_exp_pd.index.tolist())
    all_exp_pd1 = all_exp_pd1.set_index('accession_id')
    # print(all_exp_pd1)

    # return
    return all_exp_pd, all_exp_pd1


def parse_group(group_file):
    # group_info -> dict, group_name as key, list of sample names as values. {group:[s1,s2,]}
    sample_list = set()
    # group_dict = dict()
    group_dict = OrderedDict()
    with open(group_file) as f:
        header_line = f.readline()
        if not header_line.startswith('#'):
            print('Note: we assume header line exist')
        for line in f:
            if not line.strip():
                continue
            tmp_list = line.strip().split()
            sample_list.add(tmp_list[0])
            for g in tmp_list[1:]:
                group_dict.setdefault(g, list())
                group_dict[g].append(tmp_list[0])
        for g in group_dict.keys():
            group_dict[g] = sorted(list(set(group_dict[g])))#
            # 分组样本还会产生重复和空，以及排序的需求是什么
    if not group_dict:
        raise Exception("分组信息内容为空")
    return group_dict
def exp_cluster(exp_scale, n_clusters=10, s_cluster=True, g_cluster=True,
                s_cluster_type='hierarchy', g_cluster_type='hierarchy',
                s_method="complete", g_method="average",
                s_metric="correlation", g_metric="euclidean",
                use_r=True,output=None):
    result = []
    if use_r:
        hcluster = r_hcluster
    if s_cluster:
        # for sample
        if s_cluster_type == 'hierarchy':
            s_cluster_result = hcluster(exp_scale, transpose=True, n_clusters=n_clusters, method=s_method,
                                        metric=s_metric, output=output, prefix="sample.")
        else:
            s_cluster_result = kmeans(exp_scale, n_clusters=n_clusters, transpose=True, output=output, prefix='sample.')
        result.append(s_cluster_result)

    if g_cluster:
        # for gene
        if g_cluster_type == 'hierarchy':
            g_cluster_result = hcluster(exp_scale, transpose=False, n_clusters=n_clusters, method=g_method,
                                        metric=g_metric, output=output, prefix="seq.")
        else:
            g_cluster_result = kmeans(exp_scale, n_clusters=n_clusters, transpose=False, output=output, prefix='seq.')
        result.append(g_cluster_result)
    # write out exp matrix
    if output is not None:
        if not os.path.exists(output):
            os.mkdir(output)
    else:
        output = os.getcwd()
    out_dir = os.path.join(output, 'expression_matrix.xls')
    exp_scale.to_csv(out_dir, sep='\t', header=True, index=True)
    return result


if __name__ == '__main__':
    import argparse
    a = "Refer https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html for hierarchy cluster arguments"
    b = "Refer http://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html for kmeans arguments"
    parser = argparse.ArgumentParser(description=a+'; '+b)
    parser.add_argument('-exp', type=str, metavar="exp_matrix_file", required=True,
                        help="expression matrix file")
    parser.add_argument('-log_base', type=int, metavar="log_base", default=10, help="log base for transformation exp ")
    parser.add_argument('-group', type=str, metavar="group_file", default=None, help="file with two col: sample\tgroup")
    parser.add_argument('-out', type=str, default=None, help="default is local dir. Output directory name")
    parser.add_argument('-n_clusters', type=int, default=10, metavar="cluster_num", help="expected_cluster_number")
    parser.add_argument('--nsc', action='store_true', help="'no-sample-cluster', do not perform sample cluster")
    parser.add_argument('--ngc', action='store_true', help="'no-gene-cluster', do not perform gene cluster")
    parser.add_argument('-sct', metavar='sample-cluster-type', default='hierarchy', help="hierarchy or kmeans cluster")
    parser.add_argument('-gct', metavar='gene-cluster-type', default='hierarchy', help="hierarchy or kmeans cluster")
    parser.add_argument('-scm', metavar='sample-cluster-method', default='complete', help="hclust method")
    parser.add_argument('-gcm', metavar='gene-cluster-method', default='average', help="hclust method")
    parser.add_argument('-scd', metavar='sample-cluster-distance', default='correlation', help="hclust distance metric")
    parser.add_argument('-gcd', metavar='gene-cluster-distance', default='euclidean', help="hclust distance metric")
    parser.add_argument('--corr', action='store_true', help="'sample_correlation', do sample correlation")
    parser.add_argument('--pca', action='store_true', help="Perform primary component analysis")
    parser.add_argument('-group_method', metavar='group-exp-method', default='mean', help="mean median")
    parser.add_argument('-corr_method', default='pearson', help="pearson, kendall, spearman are supported")

    #
    args = parser.parse_args()
    if args.group:
        group_dict = parse_group(args.group)
    else:
        group_dict = None
    #
    exp_pd, exp_scale = process_exp_matrix(args.exp, log_base=args.log_base, group_dict=group_dict, group_method=args.group_method)
    #
    s_cluster, g_cluster = True, True
    if args.nsc:
        s_cluster = False
    if args.ngc:
        g_cluster = False
    exp_cluster(exp_scale, n_clusters=args.n_clusters,
                s_cluster=s_cluster, g_cluster=g_cluster,
                s_cluster_type=args.sct, g_cluster_type=args.gct,
                s_method=args.scm, g_method=args.gcm,
                s_metric=args.scd, g_metric=args.gcd,
                output=args.out)
    if args.corr:
        exp_cluster(exp_pd, n_clusters=args.n_clusters,
                    s_cluster=s_cluster, g_cluster=g_cluster,
                    s_cluster_type=args.sct, g_cluster_type=args.gct,
                    s_method=args.scm, g_method=args.gcm,
                    s_metric=args.scd, g_metric=args.gcd,
                    output=args.out)
        corr(exp_pd, args.corr_method, args.out)
    if args.pca:
        exp_cluster(exp_pd, n_clusters=args.n_clusters,
                    s_cluster=s_cluster, g_cluster=g_cluster,
                    s_cluster_type=args.sct, g_cluster_type=args.gct,
                    s_method=args.scm, g_method=args.gcm,
                    s_metric=args.scd, g_metric=args.gcd,
                    output=args.out)
        pca(exp_pd, output=args.out)

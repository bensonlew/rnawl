#!/mnt/ilustre/users/sanger/app/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "qindanhua"
from mako.template import Template
import os


this_file_dir = os.path.dirname(os.path.realpath(__file__))


def correlation(inputfile, corr_matrix, pvalue_out, tvalue_out, heatmap_out, col_tree, row_tree, method="pearson",
                clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", clustering_method = "complete"):
    """
    生成并运行R脚本，进行相关系数分析
    :param inputfile: 输入文件
    :param corr_matrix: 输出的相关系数矩阵
    :param pvalue_out: 输出的相关系数相应Pvalue值
    :param tvalue_out: 输出的相关系数响应tvalue值
    :param heatmap_out: 输出的相关系数热图
    :param col_tree: 输出的相关系数列的树文件
    :param row_tree: 输出的相关系数行的树文件
    :param method: 距离算法 pearson spearman
    :param clustering_distance_cols: 列计算距离方式
    :param clustering_distance_rows: 行计算距离方式
    :param clustering_method: 聚类方式
    """
    f = Template(filename=this_file_dir + '/correlation_new.r', strict_undefined=True)
    mul_test = f.render(inputfile=inputfile, corr_matrix=corr_matrix, pvalue_out=pvalue_out, tvalue_out=tvalue_out,
                        heatmap_out=heatmap_out, col_tree=col_tree, row_tree=row_tree, method=method,
                        clustering_distance_cols=clustering_distance_cols, clustering_distance_rows=clustering_distance_rows, clustering_method=clustering_method
                        )
    with open("run_correlation_new.r", 'w') as rfile:
        rfile.write("%s" % mul_test)


# def corr_heatmap(inputfile, col_tree, row_tree, col_cluster_method, row_cluster_method,
#                  col_distance_method, row_distance_method):
#     """
#     输入相关系数矩阵，输出矩阵树文件
#     :param correlation:
#     :param col_tree:
#     :param row_tree:
#     :return:
#     """
#     f = Template(filename=this_file_dir + '/corr_heatmap.r')
#     heatmap = f.render(inputfile=inputfile, col_tree=col_tree, row_tree=row_tree, col_cluster_method=col_cluster_method,
#                        row_cluster_method=row_cluster_method, col_distance_method=col_distance_method,
#                        row_distance_method=row_distance_method)
#     with open("run_corr_heatmap.r", 'w') as rfile:
#         rfile.write("%s" % heatmap)

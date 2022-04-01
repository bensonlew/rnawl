#!/mnt/ilustre/users/sanger/app/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "qiuping"


from mako.template import Template
from biocluster.config import Config
import os


def clust(input_matrix, sub_num, method="both", lognorm=10, genes_distance_method="complete", samples_distance_method = 'complete', cltype="both"):
    """
    生成并运行R脚本，进行两组样品的差异性分析，包括student T检验，welch T检验，wilcox秩和检验
    :param input_matrix: 差异表达量矩阵
    :param method: 聚类方法选择，hclust，kmeans
    :param lognorm: 对数底数选择，10或2
    :param cltype：热图排序方式
    :param sub_num: 子聚类数目
    :param distance_method: 距离矩阵算法  默认是complete(对于层次结构聚类算法来说) kmeans默认是euclidean
    """
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    f = Template(filename=this_file_dir + '/cluster.r')
    clust = f.render(input_matrix=input_matrix, sub_num=sub_num, lognorm=lognorm, genes_distance_method=genes_distance_method, samples_distance_method =  samples_distance_method, method=method, cltype=cltype)
    with open("clust.r", 'w') as rfile:
        rfile.write("%s" % clust)

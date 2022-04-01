#!/mnt/ilustre/users/sanger/app/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "qiuping"


from mako.template import Template
from biocluster.config import Config
import os


def clust(input_matrix, sub_num, method="both", lognorm=10, genes_distance_method="complete", samples_distance_method = "average", cltype="both",genes_distance_algorithm="euclidean",samples_distance_algorithm="pearson"):
    """
    生成并运行R脚本，进行两组样品的差异性分析，包括student T检验，welch T检验，wilcox秩和检验
    :param input_matrix: 差异表达量矩阵
    :param method: 聚类方法选择，hclust，kmeans
    :param lognorm: 对数底数选择，10或2
    :param cltype：热图排序方式
    :param sub_num: 子聚类数目
    :param distance_method: 距离方式  默认是complete(对于层次结构聚类算法来说) kmeans默认是euclidean
    :param distance_algorith, 距离算法，此参数只用于hclust聚类
    """
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    f = Template(filename=this_file_dir + '/cluster.r')
    clust = f.render(input_matrix=input_matrix, sub_num=sub_num, lognorm=lognorm, genes_distance_method=genes_distance_method, samples_distance_method =  samples_distance_method, method=method, cltype=cltype,
                     genes_distance_algorithm=genes_distance_algorithm,samples_distance_algorithm=samples_distance_algorithm)
    with open("clust.r", 'w') as rfile:
        rfile.write("%s" % clust)

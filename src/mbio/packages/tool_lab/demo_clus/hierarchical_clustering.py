# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from sklearn import datasets
from sklearn import cluster

iris = datasets.load_iris()
print 'iris: {}'.format(iris)
print 'iris type: {}'.format(type(iris))

sk = cluster.AgglomerativeClustering(4)
sk.fit(iris.data)
print(sk.labels_)

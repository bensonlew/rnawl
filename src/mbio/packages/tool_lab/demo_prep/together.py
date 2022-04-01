# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from sklearn import preprocessing
import numpy as np

X_train = np.array([[3.0, 9.0, 8.0, -1.0],
                    [6.0, 0.0, 0.0,  5.0],
                    [2.0, 2.0, 0.0, -3.0],
                    [7.0, 6.0, 1.0, -0.0],
                    [3.0, 7.0, 2.0,  2.0],
                    [9.0, 7.0, 4.0, -1.0]])
print 'X train:\n{}'.format(X_train)
X_scaled = preprocessing.scale(X_train)
X_minmax_scaled = preprocessing.minmax_scale(X_train)
X_maxabs_scaled = preprocessing.maxabs_scale(X_train)
X_normalized = preprocessing.normalize(X_train)
print 'X scaled:\n{}'.format(X_scaled)
print 'X minmax scaled:\n{}'.format(X_minmax_scaled)
print 'X maxabs scaled:\n{}'.format(X_maxabs_scaled)
print 'X normalized:\n{}'.format(X_normalized)

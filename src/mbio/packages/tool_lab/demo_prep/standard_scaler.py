# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from sklearn import preprocessing
import numpy as np

X_train = np.array([[3.0, 9.0, 8.0],
                    [6.0, 0.0, 0.0],
                    [2.0, 2.0, 0.0],
                    [7.0, 6.0, 1.0],
                    [3.0, 7.0, 2.0]])
print 'X train:\n{}'.format(X_train)
scaler = preprocessing.StandardScaler().fit(X_train)
print 'scaler:\n{}'.format(scaler)
print 'scaler mean:\n{}'.format(scaler.mean_)
print 'scaler std:\n{}'.format(scaler.scale_)
X_scaled = scaler.transform(X_train)
print 'X scaled:\n{}'.format(X_scaled)
X_test = np.array([[-1.0, 1.0, 0.0]])
print 'X test:\n{}'.format(X_test)
X_test_scaled = scaler.transform(X_test)
print 'X test scaled:\n{}'.format(X_test_scaled)

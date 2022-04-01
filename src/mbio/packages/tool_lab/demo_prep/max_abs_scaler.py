# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from sklearn import preprocessing
import numpy as np

X_train = np.array([[3.0, 9.0, 8.0, -1.0],
                    [6.0, 0.0, 0.0,  5.0],
                    [2.0, 2.0, 0.0, -3.0],
                    [7.0, 6.0, 1.0, -4.0],
                    [3.0, 7.0, 2.0,  2.0]])
print 'X train:\n{}'.format(X_train)
max_abs_scaler = preprocessing.MaxAbsScaler()
print 'scaler:\n{}'.format(max_abs_scaler)
X_train_maxabs = max_abs_scaler.fit_transform(X_train)
print 'X scaled:\n{}'.format(X_train_maxabs)
print 'scaler max abs: {}'.format(max_abs_scaler.max_abs_)
print 'scaler scale: {}'.format(max_abs_scaler.scale_)
X_test = np.array([[-3.0, -1.0, 4.0, 0.0]])
print 'X test:\n{}'.format(X_test)
X_test_maxabs = max_abs_scaler.transform(X_test)
print 'X test minmax:\n{}'.format(X_test_maxabs)

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
min_max_scaler = preprocessing.MinMaxScaler()
print 'scaler:\n{}'.format(min_max_scaler)
X_train_minmax = min_max_scaler.fit_transform(X_train)
print 'X scaled:\n{}'.format(X_train_minmax)
print 'scaler max: {}'.format(min_max_scaler.data_max_)
print 'scaler min: {}'.format(min_max_scaler.data_min_)
X_test = np.array([[-3.0, -1.0, 4.0]])
print 'X test:\n{}'.format(X_test)
X_test_minmax = min_max_scaler.transform(X_test)
print 'X test minmax:\n{}'.format(X_test_minmax)

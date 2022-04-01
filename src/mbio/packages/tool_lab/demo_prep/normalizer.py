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
normalizer = preprocessing.Normalizer()
print 'normalizer:\n{}'.format(normalizer)
X_train_normalized = normalizer.fit_transform(X_train)
print 'X normalized:\n{}'.format(X_train_normalized)
print 'X square root of row quadratic sum:\n{}'.format(
    (X_train ** 2).sum(axis=1) ** 0.5
)

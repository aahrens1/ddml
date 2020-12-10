# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 23:08:33 2020

@author: ecomes
"""

from sklearn.linear_model import Ridge
import numpy as np
n_samples, n_features = 100, 5
rng = np.random.RandomState(0)
y = rng.randn(n_samples)
X = rng.randn(n_samples, n_features)

y100 = y*100
clf100 = Ridge(alpha=100.0)
clf100.fit(X, y100)
# rescale coefficients and print
print(clf100.coef_/100)

y10 = y*10
clf10 = Ridge(alpha=10.0)
clf10.fit(X, y10)
# rescale coefficients and print
print(clf10.coef_/10)

clf = Ridge(alpha=1.0)
clf.fit(X, y)
print(clf.coef_)

y_10 = y/10
clf_10 = Ridge(alpha=0.1)
clf_10.fit(X, y_10)
# rescale coefficients and print
print(clf_10.coef_*10)

y_100 = y/100
clf_100 = Ridge(alpha=0.01)
clf_100.fit(X, y_100)
# rescale coefficients and print
print(clf_100.coef_*100)

y_1000 = y/1000
clf_1000 = Ridge(alpha=0.001)
clf_1000.fit(X, y_1000)
# rescale coefficients and print
print(clf_1000.coef_*1000)

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
clf = Ridge(alpha=1.0)
clf.fit(X, y)
print(clf.coef_)
y10 = y*10
clf10 = Ridge(alpha=10.0)
clf10.fit(X, y10)
# rescale coefficients and print
print(clf10.coef_/10)

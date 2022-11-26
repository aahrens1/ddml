
// https://www.stata.com/support/faqs/statistics/linear-regression-with-interval-constraints/#ex6
// https://www.stata.com/support/faqs/statistics/regression-with-interval-constraints/

*** ddml cross-fitting for the interactive model & LATE
program nnreg, eclass sortpreserve

	version 13
	syntax varlist(numeric min=2) [if] [in] [, /// 
							 	gen(name) ///
							 	double /// datatype for fitted value
								///ml /// use ML instead of NL
								VERBose ///
								SE ///
								* ///
							]

	marksample touse

	if "`verbose'"=="" local qui qui

	local yvar : word 1 of `varlist'
	local xvars : list varlist - yvar
	local k : word count `xvars'

	python: run_nnreg("`yvar'","`xvars'","`touse'")

    ereturn clear
    *ereturn post `bhat' `vhat', depname(`yvar') obs(`N') esample(`touse')
    *ereturn display
    *ereturn matrix nlbhat = `bnl'
    *ereturn local cmd _ddml_nnls
    *ereturn local predict _ddml_nnls_p

end  


python:

import sfi
import numpy as np
import __main__
from sklearn.base import TransformerMixin,BaseEstimator
from scipy.optimize import minimize 
from scipy.optimize import nnls 
from sklearn.utils import check_X_y,check_array

class ConstrLS(BaseEstimator):
    _estimator_type="regressor"
    def fit(self, X, y):

        X,y = check_X_y(X,y, accept_sparse=True)
        xdim = X.shape[1]

        #Use nnls to get initial guess
        coef0, rnorm = nnls(X,y)

        #Define minimisation function
        def fn(coef, X, y):
            return np.linalg.norm(X.dot(coef) - y)

        #Constraints and bounds
        cons = {'type': 'eq', 'fun': lambda coef: np.sum(coef)-1}
        bounds = [[0.0,1.0] for i in range(xdim)] 

        #Do minimisation
        fit = minimize(fn,coef0,args=(X, y),method='SLSQP',bounds=bounds,constraints=cons)
        self.coef_ = fit.x
        self.is_fitted_ = True
        return self
        
    def predict(self, X):
        X = check_array(X, accept_sparse=True)
        check_is_fitted(self, 'is_fitted_')
        return np.matmul(X,self.coef_)

def run_nnreg(
    yvar, # outcome
    xvars, # predictors (temp names)
    touse # sample
    ):
    
    y = np.array(sfi.Data.get(yvar,selectvar=touse))
    x = np.array(sfi.Data.get(xvars,selectvar=touse))
    if x.ndim == 1:
        x=np.reshape(x,(-1,1))

    print(ConstrLS().fit(x,y).coef_)

end
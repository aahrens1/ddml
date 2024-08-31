*! ddml v1.4.4
*! last edited: 30aug2024
*! authors: aa/ms

program _ddml_nnls
	version 16
    syntax [anything] [if] [in] [aw fw iw] , [ stata SE VERBose * ]

    local dostata = "`stata'"~=""
    if "`verbose'"=="" local qui qui
    
    if `dostata'==0 {
        // first check that python is installed
        cap python query
        if _rc > 0 {
            local dostata = 1
            `qui' di as res "Python not available; using Stata -nl-..."
        }
        else {
            // now check that sklearn is available
            cap python which sklearn
            if _rc > 0 {
                local dostata = 1
                `qui' di as res "Python module sklearn not available; using Stata -nl-..."
            }
        }
    }
    if `dostata'==0 & "`se'"~="" {
        di as res "se option available only with stata option; using Stata -nl-..."
        local dostata = 1
    }
  
    if `dostata' {
        _ddml_nnls_stata `anything' `if' `in' [`weight' `exp'], `options' `se' `verbose'
    }
    else {
        _ddml_nnls_python `anything' `if' `in' [`weight' `exp'], `options' `verbose'
    }
end

******************************************************************************

program _ddml_nnls_stata, eclass sortpreserve

    version 16
    syntax varlist(numeric min=2)            ///
        [if] [in] [aw fw iw] [,              /// 
                                    VERBose  ///
                                    SE       ///
                                    *        ///
                                    ]

    marksample touse

    if "`verbose'"=="" local qui qui

    local yvar : word 1 of `varlist'
    local xvars : list varlist - yvar
    local k : word count `xvars'
    
    ** rescale
    tempvar yvar_t
    qui sum `yvar' if `touse'
    qui gen double `yvar_t' = (`yvar' - r(mean))/r(sd) if `touse'
    foreach v of varlist `xvars' {
        tempvar v_t
        qui gen double `v_t' = (`v' - r(mean))/r(sd) if `touse'
        local xvars_t `xvars_t' `v_t'
    }

    ** create equation part 1
    local j = 1
    foreach var of varlist `xvars_t' {
        local coef b`var'
        if `j'==1 {
            local denom 1
            local nom`j' 1
            local nlcom_denom 1
            local nlcom_nom`j' 1 
        }
        else {
            local denom `denom' + exp({`coef'}) 
            local nom`j' exp({`coef'})
            local nlcom_denom `nlcom_denom' + exp(_b[`coef':_cons]) 
            local nlcom_nom`j' exp(_b[`coef':_cons]) 
        }
        local j = `j'+1
    }

    ** create equation part 2
    local j = 1
    foreach var of varlist `xvars_t' {
        if `j'==1 {
            local eq (`nom`j''/(`denom'))*`var'
            local nlcom_coef`j' (`nlcom_nom`j''/(`nlcom_denom'))
        }
        else {
            local eq `eq' + (`nom`j''/(`denom'))*`var'
            local nlcom_coef`j' (`nlcom_nom`j''/(`nlcom_denom'))
        }
        local j = `j'+1
    }

    ** estimation
    `qui' di "nl (`yvar' = `eq') if `touse' [`weight' `exp'], `options'"
    `qui' nl (`yvar_t' = `eq') if `touse' [`weight' `exp'], `options'
    tempname bnl
    mat `bnl' = e(b)
    ** fix column names to use original rather than temp varnames
    local cnames : colfullnames `bnl'
    local j = 1
    foreach v_t of varlist `xvars_t' {
        local v : word `j' of `xvars'
        local cnames : subinstr local cnames "`v_t'" "`v'"
    }
    mat coleq `bnl' = `cnames'
    mat rownames `bnl' = `yvar'
    local N = e(N)

    ** create nlcom eq
    if "`se'"!="" {
        local nlcom_eq 
        local j = 1
        foreach var of varlist `xvars_t' {
            local nlcom_eq `nlcom_eq' (`var': `nlcom_coef`j'')
            local j = `j'+1
        }

        ** get coefficients
        `qui' di "`nlcom_eq'"
         `qui' nlcom `nlcom_eq'
       
        tempname bhat vhat
        matrix `bhat' = r(b)
        matrix `vhat' = r(V)
        matrix colnames `bhat' = `xvars'
        matrix rownames `bhat' = `yvar'
        matrix colnames `vhat' = `xvars'
        matrix rownames `vhat' = `xvars'
    }
    else {
        tempname bhat  
        mat `bhat' = J(1,`k',.)
        forvalues j = 1(1)`k' {
            mat `bhat'[1,`j'] = `nlcom_coef`j''
        }
        matrix colnames `bhat' = `xvars'
        matrix rownames `bhat' = `yvar'
    }

    ereturn clear
    ereturn post `bhat' `vhat', depname(`yvar') obs(`N') esample(`touse')
    ereturn display
    ereturn matrix nlbhat = `bnl'
    ereturn local cmd _ddml_nnls
    ereturn local predict _ddml_nnls_p

end  

******************************************************************************

program _ddml_nnls_python, eclass sortpreserve

	version 16
	syntax varlist(numeric min=2)					///
		[if] [in] [aw fw iw] [,						/// 
									NOIsily			///
									finalest(name)	///
									stype(name)		///
									*				///
									]

	* defaults
	if "`noisily'"==""		local qui qui
	if "`finalest'"==""		local finalest nnls1
	if "`stype'"==""		local stype reg

	marksample touse
	qui count if `touse'
	local N = r(N)
	
	local yvar : word 1 of `varlist'
	local xvars : list varlist - yvar
	
    tempvar wt
    if "`weight'"~="" {
        qui gen double `wt' `exp'
    }
    else {
    	qui gen byte `wt' = 1
    }

	`qui' python: py_get_stack_weights("`yvar'","`xvars'","`touse'","`wt'","`finalest'","`stype'")
	tempname bhat
	// python returns column vectors
	mat `bhat' = r(b)'
	matrix colnames `bhat' = `xvars'
	matrix rownames `bhat' = `yvar'
	ereturn clear
	ereturn post `bhat', depname(`yvar') obs(`N') esample(`touse')
	ereturn scalar N = `N'
	ereturn local finalest `finalest'
	ereturn display

end

version 16
python:

import sfi
from sfi import Data,Matrix,Scalar,SFIToolkit
from sklearn.linear_model import LinearRegression,RidgeCV
from sklearn.base import TransformerMixin,BaseEstimator
from sklearn.utils import check_X_y,check_array
from sklearn.utils.validation import _check_sample_weight
import numpy as np
from scipy.optimize import minimize 
from scipy.optimize import nnls 

def py_get_stack_weights(yvar,xvars,touse,wvar,finalest,stype):

    X = Data.get(xvars,selectvar=touse)
    y = Data.get(yvar,selectvar=touse)
    w = Data.get(wvar,selectvar=touse)

    if finalest == "nnls0" and stype == "class": 
        fin_est = LinearRegressionClassifier(fit_intercept=False,positive=True)
    elif finalest == "nnls_sk" and stype == "class": 
        fin_est = LinearRegressionClassifier(fit_intercept=False,positive=True)
    elif finalest == "nnls1" and stype == "class": 
        fin_est = ConstrLSClassifier()
    elif finalest == "ridge" and stype == "class": 
        fin_est = LogisticRegression()
    elif finalest == "nnls0" and stype == "reg": 
        fin_est = LinearRegression2(fit_intercept=False,positive=True)
    elif finalest == "nnls_sk" and stype == "reg": 
        fin_est = LinearRegression2(fit_intercept=False,positive=True)
    elif finalest == "nnls1" and stype == "reg": 
        fin_est = ConstrLS()
    elif finalest == "ridge" and stype == "reg": 
        fin_est = RidgeCV2()
    elif finalest == "avg" and stype == "reg": 
        fin_est = AvgEstimator()
    elif finalest == "avg" and stype == "class": 
        fin_est = AvgClassifier()
    elif finalest == "singlebest" and stype == "reg": 
        fin_est = SingleBest()
    elif finalest == "singlebest" and stype == "class": 
        fin_est = SingleBestClassifier()
    elif finalest == "ols" and stype == "class": 
        fin_est = LinearRegressionClassifier()    
    elif finalest == "ols" and stype == "reg": 
        fin_est = LinearRegression2()    
    elif finalest == "ls1" and stype == "reg":
        fin_est = ConstrLS(unit_interval=False)    
    elif finalest == "ls1" and stype == "class":
        fin_est = ConstrLSClassifier(unit_interval=False)    
    else:
        sfi.SFIToolkit.stata('di as err "specified final estimator not supported"')
        #"
        sfi.SFIToolkit.error(198)
    fin_est.fit(X, y, w)
    b = fin_est.coef_
    Matrix.store("r(b)", b)

class ConstrLS(BaseEstimator):
    _estimator_type="regressor"
    def fit(self, X, y, w):

        X,y = check_X_y(X,y, accept_sparse=True)
        xdim = X.shape[1]

        #Use nnls to get initial guess
        #coef0, rnorm = nnls(X,y)
        #Use LinearRegression to get initial guess
        initial_est = LinearRegression(positive=True,fit_intercept=False)
        initial_est.fit(X, y, w)
        coef0 = initial_est.coef_

        #Define minimisation function
        def fn(coef, X, y):
            return np.linalg.norm(X.dot(coef) - y)
        
        #Constraints and bounds
        cons = {'type': 'eq', 'fun': lambda coef: np.sum(coef)-1}
        if self.unit_interval==True:
            bounds = [[0.0,1.0] for i in range(xdim)] 
        else:
            bounds = None

        w = _check_sample_weight(w, X)
        #If weights vector=1, no weighting needed
        if np.all(w==1):
            #Do minimisation
            fit = minimize(fn,coef0,args=(X, y),method='SLSQP',bounds=bounds,constraints=cons)
        else:
            #Use additional precision
            Xw = np.multiply(X,w,dtype=np.longdouble)
            yw = np.multiply(y,w,dtype=np.longdouble)
            #Do minimisation
            fit = minimize(fn,coef0,args=(Xw, yw),method='SLSQP',bounds=bounds,constraints=cons)
        
        self.coef_ = fit.x
        if abs(sum(self.coef_)-1)>1e-12:
            #Renormalize if coefs still don't sum to 1
            self.coef_ = self.coef_ * 1/sum(self.coef_)
        self.is_fitted_ = True
        return self
        
    def predict(self, X):
        X = check_array(X, accept_sparse=True)
        check_is_fitted(self, 'is_fitted_')
        return np.matmul(X,self.coef_)

    def __init__(self, unit_interval=True):
        self.unit_interval = unit_interval

class SingleBest(BaseEstimator):
    _estimator_type="regressor"
    def fit(self, X, y, w):
        X, y = check_X_y(X, y, accept_sparse=True)
        self.is_fitted_ = True
        ncols = X.shape[1]
        lowest_mse = np.Inf
        for i in range(ncols):
            this_mse=np.mean((y-X[:, i]) ** 2)
            if this_mse < lowest_mse:
                lowest_mse = this_mse
                best = i
        self.best = best
        coef = np.zeros(ncols)
        coef[best] = 1
        self.coef_ = coef
        self.cvalid=X
        return self
    def predict(self, X):
        X = check_array(X, accept_sparse=True)
        check_is_fitted(self, 'is_fitted_')
        return X[:,self.best]

class AvgEstimator(BaseEstimator):
    """
    Avg of learners
    """
    _estimator_type="regressor"
    def fit(self, X, y, w):
        X, y = check_X_y(X, y, accept_sparse=True)
        self.is_fitted_ = True
        ncols = X.shape[1]
        self.coef_ = np.repeat(1/ncols,ncols)
        self.cvalid=X
        return self
    def predict(self, X):
        X = check_array(X, accept_sparse=True)
        check_is_fitted(self, 'is_fitted_')
        return X.mean(axis=1)

class ConstrLSClassifier(ConstrLS):
    _estimator_type="classifier"
    def predict_proba(self, X):
        return self.predict(X)

class SingleBestClassifier(SingleBest):
    _estimator_type="classifier"
    def predict_proba(self, X):
        return self.predict(X)

class AvgClassifier(AvgEstimator):
    _estimator_type="classifier"
    def predict_proba(self, X):
        return self.predict(X)

class LinearRegressionClassifier(LinearRegression):
    _estimator_type="classifier"
    def predict_proba(self, X):
        self.cvalid=X
        return self.predict(X)

class LinearRegression2(LinearRegression):
    def fit(self,X,y,w):
        self.cvalid=X
        return LinearRegression.fit(self,X,y,w)

class RidgeCV2(RidgeCV):
    def fit(self,X,y,w):
        self.cvalid=X
        return RidgeCV.fit(self,X,y,w)

end

**** end python section

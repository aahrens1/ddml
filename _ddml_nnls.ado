*! ddml v1.2
*! last edited: 20 feb 2023
*! authors: aa/ms

program _ddml_nnls
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

program _ddml_nnls_python, eclass sortpreserve

    version 16
    syntax varlist(numeric min=2)                ///
        [if] [in] [aw fw iw] [,                  /// 
                                    VERBose      ///
                                    *            ///
                                    ]

    if "`verbose'"==""   local qui qui

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

    `qui' python: py_nnls("`yvar'","`xvars'","`touse'","`wt'")
    tempname bhat
    // python returns column vectors
    mat `bhat' = r(b)'
    matrix colnames `bhat' = `xvars'
    matrix rownames `bhat' = `yvar'
    ereturn clear
    ereturn post `bhat', depname(`yvar') obs(`N') esample(`touse')
    ereturn display    

end


version 16.0
python:

from sfi import Data,Matrix,Scalar,SFIToolkit
from sklearn.linear_model import LinearRegression
from sklearn.utils import check_X_y
from sklearn.base import BaseEstimator
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import nnls

def py_nnls(yvar,xvars,touse,wvar):

    X = Data.get(xvars,selectvar=touse)
    y = Data.get(yvar,selectvar=touse)
    w = Data.get(wvar,selectvar=touse)
    
    # standard nnls doesn't impose constraint that coefs sum to 1
    # reg_nnls = LinearRegression(positive=True,fit_intercept=False)
    # so use scipy minimize instead
    reg_nnls = ConstrLS()
    reg_nnls.fit(X, y, w)       
    b = reg_nnls.coef_
    Matrix.store("r(b)", b)

class ConstrLS(BaseEstimator):
    _estimator_type="regressor"
    def fit(self, X, y, w):

        X,y = check_X_y(X,y, accept_sparse=True)
        xdim = X.shape[1]

        #Use nnls to get initial guess
        #coef0, rnorm = nnls(X,y)
        #Use LinearRegression below to get initial guess, incorporating weights
        initial_est = LinearRegression(positive=True,fit_intercept=False)
        initial_est.fit(X, y, w)
        coef0 = initial_est.coef_

        #Define minimisation function
        def fn(coef, X, y):
            return np.linalg.norm(X.dot(coef) - y)
        
        #Constraints and bounds
        cons = {'type': 'eq', 'fun': lambda coef: np.sum(coef)-1}
        bounds = [[0.0,1.0] for i in range(xdim)]

        w = np.reshape(np.sqrt(w),(-1,1))
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

end

**** end python section

program _ddml_nnls_stata, eclass sortpreserve

    version 13
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

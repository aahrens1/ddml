*! ddml v1.2
*! last edited: 4 feb 2023
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
    syntax varlist(numeric min=2)            ///
        [if] [in] [aw fw iw] [,              /// 
                                    VERBose  ///
                                    *        ///
                                    ]

    if "`verbose'"=="" local qui qui

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
import numpy as np

def py_nnls(yvar,xvars,touse,wvar):

    X = Data.get(xvars,selectvar=touse)
    y = Data.get(yvar,selectvar=touse)
    w = Data.get(wvar,selectvar=touse)
    
    reg_nnls = LinearRegression(positive=True,fit_intercept=False)
    reg_nnls.fit(X, y, w)
    # renormalize:
    b = reg_nnls.coef_ * 1/sum(reg_nnls.coef_)
    Matrix.store("r(b)", b)

end


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


program _ddml_optiv, eclass

	version 13

	syntax [anything] , yvar(varname) ///
				dvar(varlist) ///
				dtilde(varlist) ///
                dhtilde(varlist) ///
                ytilde(varname) ///
				[touse(varname)]

    if "`touse'"=="" {
        tempvar touse
        mark `touse'
        markout `touse' `yvar' `dvar' `dtilde' `dhtilde' `ytilde'
    }

    tempname b
    tempname V 

    tempvar zvar 
    tempvar dx
    gen double `zvar' = `dhtilde'-`dtilde'
    gen double `dx' = `dvar'-`dtilde'

    qui ivreg2 `yvar' (`dx'=`zvar') if `touse', nocons `robust' noheader nofooter

    mat `b' = e(b)
    mat `V' = e(V)

    // display
    local N = `e(N)'
    matrix colnames `b' = "`dvar'"
    matrix rownames `b' = "`yvar'"
    matrix colnames `V' = "`dvar'"
    matrix rownames `V' = "`dvar'"
    ereturn clear
    tempvar esample
    // ereturn post esample(vname) moves - doesn't copy! - vname into e(sample), so create a copy
    qui gen byte `esample'=`touse'
    ereturn post `b' `V', depname(`yvar') obs(`N') esample(`esample')
    ereturn display

end


program _ddml_optiv, eclass

	version 13

	syntax [anything] , yvar(varname) ///
				dvar(varname) ///
				dtilde(varname) ///
                dhtilde(varname) ///
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
    ereturn post `b' `V', depname(`yvar') obs(`N') esample(`touse')
    ereturn display

end

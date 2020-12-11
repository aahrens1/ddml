
program _ddml_optiv, eclass

	version 13

	syntax [anything] , yvar(varname) ///
				dvar(varlist) ///
				dtilde(varlist) ///
                dhtilde(varlist) ///
                ytilde(varname) ///
				[touse(varname) ///
                debug]

    if "`touse'"=="" {
        tempvar touse
        mark `touse'
        markout `touse' `yvar' `dvar' `dtilde' `dhtilde' `ytilde'
    }

    tempname b
    tempname V 

    ** generate IVs
    local NumEndog : word count `dtilde'
    local NumInstr : word count `dhtilde'

    if ("`debug'"~="") {
        di "Names: `dtilde'"
        di "Names: `dhtilde'"
        di "Names: `dvar'"
    }

    local dlist
    local zlist
    forvalues i = 1/`NumEndog' {
        tempvar zvar`i' 
        tempvar dx`i'
        local dh : word `i' of `dhtilde'
        local dt : word `i' of `dtilde'
        local dd : word `i' of `dvar'
        if ("`debug'"~="") {
            di "endog regressor `i' = `dd'-`dt'"
            di "instr `i' = `dh'-`dt'"
        }
        gen double `zvar`i'' = `dh'-`dt'
        gen double `dx`i'' = `dd'-`dt' 
        local dlist `dlist' `zvar`i''
        local zlist `zlist' `dx`i''
    }

    qui ivreg2 `ytilde' (`dlist'=`zlist') if `touse', nocons `robust' noheader nofooter

    mat `b' = e(b)
    mat `V' = e(V)

    // display
    local N = `e(N)'
    matrix colnames `b' = `dvar'
    matrix rownames `b' = `yvar'
    matrix colnames `V' = `dvar'
    matrix rownames `V' = `dvar'
    ereturn clear
    tempvar esample
    // ereturn post esample(vname) moves - doesn't copy! - vname into e(sample), so create a copy
    qui gen byte `esample'=`touse'
    ereturn post `b' `V', depname(`yvar') obs(`N') esample(`esample')
    ereturn display

end

* notes:
* is detailed message on creation of orthogonalized D and optimal IV needed?

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
	* model is exactly identified
	local NumEndog : word count `dtilde'
	local NumInstr : word count `dhtilde'
	if (`NumEndog'!=`NumInstr') {
		di as err "something is wrong. Number of IVs != number of endog. regressors"
	}

	local dlist
	local zlist
	forvalues i = 1/`NumEndog' {
		tempvar zvar`i' 
		tempvar dx`i'
		local dh : word `i' of `dhtilde'
		local dt : word `i' of `dtilde'
		local dd : word `i' of `dvar'
//		if `i'==1 {
			local ntext1 ="(`i') "
			local ntext2 ="    "
//		}
//		else {
//			local ntext1
//			local ntext2
//		}
		di as text "`ntext1'Orthogonalised `dd' created as `dd'-`dt'."
		di as text "`ntext2'Optimal instrument for `dd' created as `dh'-`dt'."
		gen double `zvar`i'' = `dh'-`dt' // E[D|ZX]-E[D|X] = instrument
		gen double `dx`i'' = `dd'-`dt' // D-E[D|X] = endogenous regressor
		local dlist `dlist' `dx`i''
		local zlist `zlist' `zvar`i''
	}

	qui regress `ytilde' `dlist' (`zlist') if `touse', nocons `robust'

	mat `b' = e(b)
	mat `V' = e(V)

	// display
	local N = e(N)
	matrix colnames `b' = `dvar'
	matrix rownames `b' = `yvar'
	matrix colnames `V' = `dvar'
	matrix rownames `V' = `dvar'
	ereturn clear
	tempvar esample
	// ereturn post esample(vname) moves - doesn't copy! - vname into e(sample), so create a copy
	qui gen byte `esample'=`touse'
	ereturn post `b' `V', depname(`yvar') obs(`N') esample(`esample')

end

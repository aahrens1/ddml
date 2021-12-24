// does OLS/IV and reports with substitute yname and dnames
program define _ddml_reg, eclass
	syntax [anything] [if] [in] , [ y(varname) yname(name) d(varlist) dnames(namelist) robust * ]

	if ~replay() {
		// estimate
		
		marksample touse
	
		qui reg `y' `d' if `touse', `robust' `options'
	
		tempname b
		tempname V
		mat `b' = e(b)
		mat `V' = e(V)
		matrix colnames `b' = `dnames'
		matrix rownames `b' = `yname'
	 	matrix colnames `V' = `dnames'
		matrix rownames `V' = `dnames'
		local N = e(N)
		ereturn clear
		ereturn post `b' `V', depname(`yname') obs(`N') esample(`touse')
		if "`robust'"=="robust" {
			ereturn local vce		robust
			ereturn local vctype	Robust
		}
		ereturn local cmd _ddml_reg
		ereturn local y `y'
		ereturn local d `d'
	}
	
	di as text "E[y|X] = " as res "`e(y)'" _c
	di as text _col(52) "Number of obs   =" _col(70) as res %9.0f `e(N)'
	di as text "E[D|X] = " as res "`e(d)'"
	ereturn display				

end

*! ddml v1.1
*! last edited: 28 dec 2022
*! authors: aa/ms

*** ddml cross-fitting for the interactive model & LATE
program _ddml_nnls, eclass sortpreserve

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

	** create equation part 1
	local j = 1
	foreach var of varlist `xvars' {
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
	foreach var of varlist `xvars' {
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
	`qui' di "nl (`yvar' = `eq') if `touse', `options'"
	`qui' nl (`yvar' = `eq') if `touse', `options'
	tempname bnl
	mat `bnl' = e(b)
	local N = e(N)

	** create nlcom eq
	if "`se'"!="" {
		local nlcom_eq 
		local j = 1
		foreach var of varlist `xvars' {
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
	    matrix rownames `bhat' = "`yvar'"
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
	    matrix rownames `bhat' = "`yvar'"
	}

    ereturn clear
    ereturn post `bhat' `vhat', depname(`yvar') obs(`N') esample(`touse')
    ereturn display
    ereturn matrix nlbhat = `bnl'
    ereturn local cmd _ddml_nnls
    ereturn local predict _ddml_nnls_p

end  

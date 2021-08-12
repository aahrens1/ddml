
*** ddml cross-fitting for the interactive model & LATE
program _ddml_nnls, eclass sortpreserve

	syntax [anything] [if] [in] , /// 
							[ * ///
							]

	marksample touse

	local yvar : word 1 of `anything'
	local xvars : list anything - yvar
	local k : word count `xvars'

	** create equation part 1
	local j = 1
	foreach var of varlist `xvars' {
		local coef b`var'
		if `j'==1 {
			local denom 1
			local nom`j' = 1 
		}
		else {
			local denom `denom' + exp({`coef'}) 
			local nom`j' exp({`coef'}) 
		}
		local j = `j'+1
	}

	** create equation part 2
	local j = 1
	foreach var of varlist `xvars' {
		if `j'==1 {
			local eq (`nom`j''/(`denom'))*`var'
		}
		else {
			local eq `eq' + (`nom`j''/(`denom'))*`var'
		}
		local j = `j'+1
	}

	** estimation
	qui nl (`yvar' = `eq'), `options'
	tempname bnl
	mat `bnl' = e(b)

	** calculate denominator
	local j = 1
	local denomhat = 1
	foreach var of varlist `xvars' {
		if `j'>1 {
			local thisbeta = `bnl'[1,`j'-1]
			local denomhat = `denomhat' + exp(`thisbeta')
		}
		local j = `j'+1
	} 

	** obtain coefficients
	local j = 1
	tempname bhat
	mat `bhat' = J(1,`k',.)
	foreach var of varlist `xvars' {
		if `j'==1 {
			mat `bhat'[1,1] = 1 / `denomhat'
		}
		else {
			local thisbeta = `bnl'[1,`j'-1]
			mat `bhat'[1,`j'] = exp(`thisbeta') / `denomhat'
		}
		local j = `j'+1
	} 
 
    local N = `e(N)'
    matrix colnames `bhat' = `xvars'
    matrix rownames `bhat' = "`yvar'"
    ereturn clear
    ereturn post `bhat', depname(`yvar') obs(`N')
    ereturn display
end  

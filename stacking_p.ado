
program define stacking_p
	version 14
	
	syntax namelist(min=1 max=2) [if] [in] , [ xb Residuals NOIsily TRANSForm]
 
	if "`noisily'"=="" {
		local qui quietly
	}
	
	// prediction sample
	marksample touse, novarlist
	
	// estimation sample
	tempvar esample
	qui gen byte `esample' = e(sample)
	
	local nopts : word count `xb' `residuals'
	if `nopts' >1 {
		display "{err}only one statistic may be specified"
		exit 498
	}
	if `nopts' == 0 {
		local xb xb
	}
	
	// get var type & name
	tokenize `namelist'
	if "`2'"=="" {						//  only new varname provided
		local varname `1'
	}
	else {								//  datatype also provided
		local vtype `1'
		local varname `2'
	}

	local command=e(cmd)
	if ("`command'"~="stacking") {
		di as err "error: -stacking_p- supports only the -stacking- command"
		exit 198
	}
	local depvar `e(depvar)'
	local eqn `e(eqn)'
	
	local nlearners		`e(mcount)'
	local vtlist		`e(base_est)'
	
	// predicted values
	tempvar xbstacked
	qui gen `vtype' `xbstacked' = 0 if `touse'
	
	// stacking weights
	tempname s_weights
	mat `s_weights'		=e(weights)
	
	// learners will reset e(.) macros, so store them
	tempname ehold
	_estimates hold `ehold', copy
	
	forvalues i=1/`nlearners' {
		local vtilde : word `i' of `vtlist'
		_estimates unhold `ehold'
		local 0 `e(estring_`i')'
		_estimates hold `ehold', copy
		syntax [anything] [if] [in] , [*]
		local est_main `anything'
		local est_options `options'
		// estimate using estimation sample
		`qui' `est_main' if `esample', `est_options'
		tempvar xbv
		// use all observations
		qui predict `vtype' `xbv'
		if "`transform'"=="" {
			qui replace `xbstacked' = `xbstacked' + `xbv'*`s_weights'[`i',1]
		}
		else if "`xb'" != "" {
			qui replace `xbv' = . if ~`touse'
			qui rename `xbv' `varname'`i'
			label var `varname'`i' "Prediction `vtilde'"
			qui count if `varname'`i'==.
			if r(N)>0 {
				di as text "(`r(N)' missing values generated)"
			}
		}
		else {
			qui generate `vtype' `varname'`i' = `depvar' - `xbv' if `touse'
			label var `varname'`i' "Residuals `vtilde'"
			qui count if `varname'`i'==.
			if r(N)>0 {
				di as text "(`r(N)' missing values generated)"
			}
		}
	}
	
	// restore e(.) macros
	_estimates unhold `ehold'
	
	if "`transform'"=="" {
		if "`xb'" != "" {
			qui rename `xbstacked' `varname'
			label var `varname' "Prediction stacking"
		}
		else {
			qui generate `vtype' `varname' = `depvar' - `xbstacked' if `touse'
			label var `varname' "Residuals stacking"
		}
		qui count if `varname'==.
		if r(N)>0 {
			di as text "(`r(N)' missing values generated)"
		}
	}
	

end

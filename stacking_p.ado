
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
		local varlist `1'
	}
	else {								//  datatype also provided
		local vtype `1'
		local varlist `2'
	}

	local command=e(cmd)
	if ("`command'"~="stacking") {
		di as err "error: -stacking_p- supports only the -stacking- command"
		exit 198
	}
	local depvar `e(depvar)'
	local eqn `e(eqn)'
	
	mata: st_local("nlearners", strofreal(`eqn'.nlearners))
	mata: st_local("vtlist", invtokens(`eqn'.vtlist))
	
	// predicted values
	tempvar xbstacked
	`qui' gen double `xbstacked' = 0 if `touse'
	
	// stacking weights
	tempname ss_weights
	mata: st_matrix("`ss_weights'",return_result_item(`eqn',"`e(stack)'","ss_weights", "1"))
	
	// learners will reset e(.) macros, so store them
	tempname ehold
	_estimates hold `ehold', copy
	
	forvalues i=1/`nlearners' {
		local vtilde : word `i' of `vtlist'
		mata: st_local("est_main", return_learner_item(`eqn',"`vtilde'","est_main"))
		mata: st_local("est_options", return_learner_item(`eqn',"`vtilde'","est_options"))
		// estimate using estimation sample
		`qui' `est_main' if `esample', `est_options'
		tempvar xbv
		// use all observations
		`qui' predict double `xbv'
		if "`transform'"=="" {
			`qui' replace `xbstacked' = `xbstacked' + `xbv'*`ss_weights'[1,`i']
		}
		else if "`xb'" != "" {
			`qui' replace `xbv' = . if ~`touse'
			`qui' rename `xbv' `varlist'`i'
			label var `varlist'`i' "Prediction `vtilde'"
			`qui' count if `varlist'`i'==.
			di as text "(`r(N)' missing values generated)"
		}
		else {
			`qui' generate `typlist' `varlist'`i' = `depvar' - `xbv' if `touse'
			label var `varlist'`i' "Residuals `vtilde'"
			`qui' count if `varlist'`i'==.
			di as text "(`r(N)' missing values generated)"
		}
	}
	
	// restore e(.) macros
	_estimates unhold `ehold'
	
	if "`transform'"=="" {
		if "`xb'" != "" {
			qui rename `xbstacked' `varlist'
		}
		else {
			qui generate `typlist' `varlist' = `depvar' - `xbstacked' if `touse'
		}
		`qui' count if `varlist'==.
		di as text "(`r(N)' missing values generated)"
	}
	

end

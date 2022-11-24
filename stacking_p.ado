
program define stacking_p
	version 14
	
	syntax namelist(min=1 max=2) [if] [in] , [ xb Residuals NOIsily CValid TRANSForm]

	if "`noisily'"=="" {
		local qui quietly
	}
	
	local cvalidflag	=("`cvalid'"~="")
	// check
	if `cvalidflag' & e(cvalid)==0 {
		di as err "error - cross-validated OOS predictions not available; re-estimate using cvalid option"
		exit 198
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
	local base_est		`e(base_est)'
	local yhat_list		`e(base_yhat)'
	local cvalid_list	`e(base_cv)'
	
	// stacking weights
	tempname s_weights
	mat `s_weights'		=e(weights)
	
	// predicted values
	if "`transform'"=="" {
		// initialize stacked predictions; will accumulate by looping through learners
		qui gen `vtype' `varname' = 0
	}
	forvalues i=1/`nlearners' {
		local yhat			: word `i' of `yhat_list'
		local yhat_cv	: word `i' of `cvalid_list'
		local vtilde		: word `i' of `base_est'
		if "`transform'"=="" & `cvalidflag'==0 {
			// stacked predictions
			qui replace `varname' = `varname' + `yhat'*`s_weights'[`i',1]
		}
		else if "`transform'"=="" {
			// stacked CV OOS predictions
			qui replace `varname' = `varname' + `yhat_cv'*`s_weights'[`i',1]
		}
		else {
			// base learner predictions
			if `cvalidflag' {
				qui gen `vtype' `varname'`i' = `yhat_cv'
			}
			else {
				qui gen `vtype' `varname'`i' = `yhat'
			}
			if "`residuals'"~="" {
				qui replace `varname'`i' = `depvar' - `varname'`i'
				label var `varname'`i' "Residuals `vtilde'"
			}
			else {
				label var `varname'`i' "Prediction `vtilde'"
			}
			qui replace `varname'`i' = . if ~`touse'
			qui count if `varname'`i'==.
			if r(N)>0 {
				di as text "(`r(N)' missing values generated)"
			}
		}
	}
	
	// stacked predictions - finish up
	if "`transform'"=="" {
		if "`xb'" != "" {
			label var `varname' "Prediction stacking"
		}
		else {
			qui replace `varname' = `depvar' - `varname'
			label var `varname' "Residuals stacking"
		}
		qui count if `varname'==.
		if r(N)>0 {
			di as text "(`r(N)' missing values generated)"
		}
	}
	

end

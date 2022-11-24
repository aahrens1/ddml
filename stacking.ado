* add check that depvar matches dep vars of learners

program define stacking, eclass
	version 14
	
	if ~replay() {
		_stacking `0'
	}
	
	// save for display results
	tempname weights_mat
	mat `weights_mat'=e(weights)
	local base_est `e(base_est)'
	local nlearners	= e(mcount)
	
	// display results
	di
	di as text "{hline 17}{c TT}{hline 21}"
	di as text "  Method" _c
	di as text _col(18) "{c |}      Weight"
	di as text "{hline 17}{c +}{hline 21}"
	forvalues j=1/`nlearners' {
		local b : word `j' of `base_est'
		di as text "  `b'" _c
		di as text _col(18) "{c |}" _c
		di as res %15.7f el(`weights_mat',`j',1)
	}

	// parse and check for graph/table options
	// code borrowed from pystacked - needed to accommodate syntax #2
	if ~replay() {
		tokenize "`0'", parse(",")
		local beforecomma `1'
		macro shift
		local restargs `*'
		tokenize `beforecomma', parse("|")
		local mainargs `1'
		local 0 "`mainargs' `restargs'"
	}
	syntax [anything]  [if] [in] [aweight fweight] , 	///
				[										///
					GRAPH1								/// vanilla option, abbreviates to "graph"
					HISTogram							/// report histogram instead of default ROC
					graph(string asis)					/// for passing options to graph combine
					lgraph(string asis)					/// for passing options to the graphs of the learners
					TABle								/// 
					HOLDOUT1							/// vanilla option, abbreviates to "holdout"
					holdout(varname)					///
					CValid								///
					*									///
				]
	
	// graph/table block
	if `"`graph'`graph1'`lgraph'`histogram'`table'"' ~= "" {
		stacking_graph_table,							///
			`holdout1' holdout(`holdout')				///
			`cvalid'									///
			`graph1'									///
			`histogram'									///
			goptions(`graph') lgoptions(`lgraph')		///
			`table'
	}
	
	// print RMSPE table
	if "`table'" ~= "" {
		tempname m w
		mat `m' = r(m)
			
		di
		di as text "RMSPE: In-Sample, CV, Holdout"
		di as text "{hline 17}{c TT}{hline 47}"
		di as text "  Method" _c
		di as text _col(18) "{c |} Weight   In-Sample        CV         Holdout"
		di as text "{hline 17}{c +}{hline 47}"
		
		di as text "  STACKING" _c
		di as text _col(18) "{c |}" _c
		di as text "    .  " _c
		di as res  _col(30) %7.3f el(`m',1,1) _col(43) %7.3f el(`m',1,2) _col(56) %7.3f el(`m',1,3)
		
		forvalues j=1/`nlearners' {
			local b : word `j' of `base_est'
			di as text "  `b'" _c
			di as text _col(18) "{c |}" _c
			di as res _col(20) %5.3f el(`weights_mat',`j',1) _c
			di as res _col(30) %7.3f el(`m',`j'+1,1) _col(43) %7.3f el(`m',`j'+1,2) _col(56) %7.3f el(`m',`j'+1,3)
		}

		// add to estimation macros
		ereturn mat rmspe = `m'
	}
	
end

program define _stacking, eclass
	syntax varname [if] [in] ,						/// 
							[						///
							foldvar(varlist)		/// one fold var per resample
							kfolds(integer 5)		/// ignored if foldvars provided
							NORANDOM				/// first fold ID uses obs in existing order
							estring(string asis)	/// estimation string
							CValid					///
							replace					///
							NOIsily					///
							NOPREFIX				/// no prefix appended
							prefix(name)			///
							*						///
							]

	// check that ddml is installed
	cap ddml, version
	if _rc>0 {
		di as err "error - stacking requires ddml package to be installed"
		exit 199
	}
	
	// for clarity
	local depvar	`varlist'
	local cvalidflag	=("`cvalid'"~="")
	
	// default prefix; overrides prefix option
	if "`noprefix'"~="" {
		local prefix
	}
	else if "`prefix'"=="" {
		local prefix _stack_
	}

	if "`noisily'"=="" {
		local qui quietly
	}
	
	marksample touse
	// if dep var is missing, automatically not in estimation sample
	markout `touse' `depvar'

	// parse estring
	local doparse = 1
	local vnum = 0
	// copy of estring
	local testring `estring'
	// vt_list is the list of learners, each prefixed by Y1_, Y2_, ...
	// yhat_list is the same with optional prefix `prefix' and with suffix _yhat
	// cvalid_list is the same with optional prefix and with the suffix _cv
	// tvt_list is the list of temp names that will have the CV predictions
	while `doparse' {
		local ++vnum
		tokenize `"`testring'"', parse("||")
		
		// catch special case - a single | appears inside the estimation string
		if "`2'"=="|" & "`3'"~="|" {
			local firstpart `1' `2' `3'
			mac shift 2
			local testring `*'
			tokenize `"`testring'"', parse("||")
			local 1 `firstpart'
		}
		local vtilde		: word 1 of `1'
		local vtilde		Y`vnum'_`vtilde'
		local vt_list		`vt_list' `vtilde'
		local yhat_list		`yhat_list' `prefix'`vtilde'_yhat
		local cvalid_list	 `cvalid_list' `prefix'`vtilde'_cv
		// tvt is tempname for cvalid predictions; tvh is for full-sample learner predictions
		tempvar tvt tvh
		local tvt_list `tvt_list' `tvt'
		local tvh_list `tvh_list' `tvh'
		// save estimation string in a local
		local estring_`vnum' `1'
		if "`2'"~="|" & "`3'"~="|" {
			// done parsing
			local doparse = 0
		}

		mac shift 3
		local testring `*'
	}
	
	// vnum is the number of learners
	local nlearners `vnum'
	
	`qui' crossfit if `touse',		///
		estring(`estring')			///
		kfolds(`kfolds')			///
		`norandom'					///
		generate(`tvt_list')		///
		`noisily'
	
	// crossfit appends the resample number (here, _1) so remove it
	foreach vn in `tvt_list' {
		rename `vn'_1 `vn'
	}

	// get saved results from crossfit and store for later posting
	// returned but not re-saved: N_folds mse_folds 
	tempname N_folds_list mse_folds_list N_list mse_list
	local N					=r(N)
	foreach m in N_folds_list mse_folds_list N_list mse_list {
		mat ``m''	=r(`m')
		mat rownames ``m'' = `vt_list'
	}
	
	`qui' di
	`qui' di as text "Stacking NNLS (additive model):"
	`qui' _ddml_nnls `depvar' `tvt_list'
	tempname s_weights
	mat `s_weights' = e(b)
	mat `s_weights' = `s_weights''
	mat rownames `s_weights' = `vt_list'
	
	// save cross-validated predictions
	if `cvalidflag' {
		forvalues i=1/`nlearners' {
			local vn	: word `i' of `cvalid_list'
			local tvn	: word `i' of `tvt_list'
			cap gen double `vn' = `tvn'
			if _rc>0 & "`replace'"=="" {
				di as err "error - variable `vn' already exists; use replace option"
				exit 110
			}
			else if _rc>0 {
				drop `vn'
				qui gen double `vn' = `tvn'
			}
		}
	}
	
	// estimate individual base learners and obtain predicted values
	forvalues i=1/`nlearners' {
		local vn	: word `i' of `yhat_list'
		local tvh	: word `i' of `tvh_list'
		local 0 `estring_`i''
		syntax [anything] [if] [in] , [*]
		local est_main `anything'
		local est_options `options'
		// before estimating, check whether variable already exists and drop if it does
		// because of possible varname abbreviations, do this using replace
		cap replace `vn' = .
		if _rc==0 {
			// variable exists
			if `i'==1 {
				di as text "warning - overwriting existing variables with base learner predictions"
			}
			drop `vn'
		}
		// estimate using estimation sample
		`qui' `est_main' if `touse', `est_options'
		// use all observations
		qui predict double `vn'
	}
	
	ereturn clear
	ereturn post, depname(`depvar') obs(`N') esample(`touse')
	ereturn local cmd				stacking
	ereturn local predict			stacking_p
	ereturn local base_yhat			`yhat_list'
	ereturn local base_cv			`cvalid_list'
	ereturn local base_est			`vt_list'
	ereturn scalar mcount			=`nlearners'
	foreach m in N_folds_list N_list mse_folds_list mse_list {
		ereturn matrix `m'			=``m''
	}
	ereturn matrix weights			=`s_weights'
	ereturn scalar cvalid			=`cvalidflag'
	forvalues i=`nlearners'(-1)1 {
		ereturn local estring_`i'	`estring_`i''
	}
	ereturn local estring			`estring'

end


// graph and/or table
program define stacking_graph_table, rclass
	version 14
	syntax ,								///
				[							///
					HOLDOUT1				/// vanilla option, abbreviates to "holdout"
					holdout(varname)		///
					CValid					/// cross-validated
					GRAPH1					/// vanilla option, abbreviates to "graph"
					HISTogram				/// report histogram instead of default ROC
					goptions(string asis)	///
					lgoptions(string asis)	///
					table					/// 
				]

	// any graph options implies graph
	local graphflag = `"`graph1'`goptions'`lgoptions'"'~=""
	
	local cvalidflag	=("`cvalid'"~="")

	if "`holdout'`holdout1'"=="" & `cvalidflag'==0 {
		local title In-sample
		tempvar touse
		qui gen `touse' = e(sample)
	}
	else if "`holdout'`holdout1'"=="" {
		local title CV
		tempvar touse
		qui gen `touse' = e(sample)
	}
	else {
		local title Holdout
		// holdout variable provided, or default = not-in-sample?
		if "`holdout'"=="" {
			// default
			tempvar touse
			qui gen `touse' = 1-e(sample)
			// check number of holdout obs
			qui count if `touse'
			if r(N)==0 {
				di as err "error - no observations in holdout sample"
				exit 198
			}
			di
			di as text "Number of holdout observations:" as res %5.0f r(N)
		}
		else {
			// check that holdout variable doesn't overlap with e(sample)
			qui count if e(sample) & `holdout' > 0
			if r(N) > 0 {
				di as err "error - holdout and estimation samples overlap"
				exit 198
			}
			qui count if `holdout' > 0 & `holdout' < .
			if r(N) == 0 {
				di as err "error - no observations in holdout sample"
				exit 198
			}
			di
			di as text "Number of holdout observations:" as res %5.0f r(N)
			local touse `holdout'
		}
	}

	local nlearners	= e(mcount)
	local learners	`e(base_est)'
	local y			`e(depvar)'
	// weights
	tempname weights
	mat `weights'	= e(weights)

	// complete graph title
	local title `title' Predictions
		
	tempvar stacking_p stacking_r cvalid_p cvalid_r
	// stacked
	qui predict double `stacking_p'
	label var `stacking_p' "Prediction: Stacking Regressor"
	qui gen double `stacking_r' = `y' - `stacking_p'
	if e(cvalid) {
		qui predict double `cvalid_p', cvalid
	}
	else {
		qui gen `cvalid_p' = .
	}
	label var `cvalid_p' "Prediction (CV): Stacking Regressor"
	qui gen double `cvalid_r' = `y' - `cvalid_p'
	// base learners
	qui predict double `stacking_p', transform
	forvalues i=1/`nlearners' {
		local lname : word `i' of `learners'
		label var `stacking_p'`i' "Prediction: `lname'"
		tempvar stacking_r`i'
		qui gen double `stacking_r`i'' = `y' - `stacking_p'`i'
	}
	if e(cvalid) {
		qui predict double `cvalid_p', transform cvalid
	}
	else {
		forvalues i=1/`nlearners' {
			qui gen `cvalid_p'`i'=.
		}
	}
	forvalues i=1/`nlearners' {
		local lname : word `i' of `learners'
		label var `cvalid_p'`i' "Prediction (CV): `lname'"
		tempvar cvalid_r`i'
		qui gen double `cvalid_r`i'' = `y' - `cvalid_p'`i'
	}

	// graph variables
	if `cvalidflag'==0 {
		local xvar stacking_p
	}
	else {
		local xvar cvalid_p
	}
	tempname g0
	if `graphflag' {
		twoway (scatter ``xvar'' `y') (line `y' `y') if `touse'		///
			,															///
			legend(off)													///
			title("STACKING")											///
			`lgoptions'													///
			nodraw														///
			name(`g0', replace)
		local glist `g0'
		forvalues i=1/`nlearners' {
			tempname g`i'
			local lname : word `i' of `learners'
			local w : di %5.3f el(`weights',`i',1)
			twoway (scatter ``xvar''`i' `y') (line `y' `y') if `touse'	///
				,															///
				legend(off)													///
				title("Learner: `lname'")									///
				`lgoptions'													///
				subtitle("weight = `w'")									///
				nodraw														///
				name(`g`i'', replace)
			local glist `glist' `g`i''
		}
	
		graph combine `glist'										///
						,											///
						title("`title'")							///
						`goptions'
	}
	
	if "`table'"~="" {
		
		// save in matrix
		tempname m m_in m_cv m_out
		// column for in-sample RMSPE
		qui sum `stacking_r' if e(sample)
		mat `m_in' = r(sd) * sqrt( (r(N)-1)/r(N) )
		forvalues i=1/`nlearners' {
			qui sum `stacking_r`i'' if e(sample)
			mat `m_in' = `m_in' \ (r(sd) * sqrt( (r(N)-1)/r(N) ))
		}
		
		// column for CV RMSPE
		qui sum `cvalid_r' if e(sample)
		mat `m_cv' = r(sd) * sqrt( (r(N)-1)/r(N) )
		forvalues i=1/`nlearners' {
			qui sum `cvalid_r`i'' if e(sample)
			mat `m_cv' = `m_cv' \ (r(sd) * sqrt( (r(N)-1)/r(N) ))
		}
		
		// column for holdout RMSPE
		if "`holdout'`holdout1'"~="" {
			// touse is the holdout indicator
			qui sum `stacking_r' if `touse'
			mat `m_out' = r(sd) * sqrt( (r(N)-1)/r(N) )
			forvalues i=1/`nlearners' {
				qui sum `stacking_r`i'' if `touse'
				mat `m_out' = `m_out' \ (r(sd) * sqrt( (r(N)-1)/r(N) ))
			}
		}
		else {
			mat `m_out' = J(`nlearners'+1,1,.)
		}
		
		mat `m' = `m_in' , `m_cv', `m_out'
		mat colnames `m' = RMSPE_in RMSPE_cv RMSPE_out
		mat rownames `m' = STACKING `learners'
		
		return matrix m = `m'

	}
	
end


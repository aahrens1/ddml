program define stacking, eclass
	
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
					*									///
				]
	
	// graph/table block
	if `"`graph'`graph1'`lgraph'`histogram'`table'"' ~= "" {
		stacking_graph_table,							///
			`holdout1' holdout(`holdout')				///
			`graph1'									///
			`histogram'									///
			goptions(`graph') lgoptions(`lgraph')		///
			`table'
	}
	
	// print MSPE table
	if "`table'" ~= "" {
		tempname m w
		mat `m' = r(m)
			
		di
		di as text "MSPE: In-Sample and Out-of-Sample"
		di as text "{hline 17}{c TT}{hline 35}"
		di as text "  Method" _c
		di as text _col(18) "{c |} Weight   In-Sample   Out-of-Sample"
		di as text "{hline 17}{c +}{hline 35}"
		
		di as text "  STACKING" _c
		di as text _col(18) "{c |}" _c
		di as text "    .  " _c
		di as res  _col(30) %7.3f el(`m',1,1) _col(44) %7.3f el(`m',1,2)
		
		forvalues j=1/`nlearners' {
			local b : word `j' of `base_est'
			di as text "  `b'" _c
			di as text _col(18) "{c |}" _c
			di as res _col(20) %5.3f el(`weights_mat',`j',1) _c
			di as res _col(30) %7.3f el(`m',`j'+1,1) _col(44) %7.3f el(`m',`j'+1,2)
		}

		// add to estimation macros
		ereturn mat mspe = `m'
	}
	
end

program define _stacking, eclass
	syntax varname [if] [in] ,						/// 
							[						///
							foldvar(varlist)		/// one fold var per resample
							kfolds(integer 5)		/// ignored if foldvars provided
							estring(string asis)	/// estimation string
							vtilde(namelist)		/// name(s) of fitted variable(s)
							stack(name)				/// varname for stacking predictions (optional)
							eqn(name)				///
							NOIsily					///
							]

	// renaming for clarity
	local vtlist `vtilde'
	// clear the local macro
	local vtilde
	
	if "`noisily'"=="" {
		local qui quietly
	}
	
	marksample touse
	// if dep var is missing, automatically not in estimation sample
	markout `touse' `vname'
	
	if "`stack'"=="" {
		// name for in-sample stacking predictions not provided, so drop when done
		local stack __stacking
		local dropstack 1
	}
	else {
		local dropstack 0
	}

	// blank eqn - declare this way so that it's a struct and not transmorphic
	if "`eqn'" == "" {
		local eqn stacking_eqn
	}
	mata: `eqn' = init_eStruct()
	
	tempname mname
	ddml init partial, mname(`mname')
	_ddml_sample if `touse' , mname(`mname')
	// sample, id, fid not needed
	cap drop `mname'*
	
	local doparse = 1
	local vnum = 1
	local hasvtlist = "`vtlist'"~=""
	while `doparse' {
		tokenize `"`estring'"', parse("||")
		
		// catch special case - a single | appears inside the estimation string
		if "`2'"=="|" & "`3'"~="|" {
			local firstpart `1' `2' `3'
			mac shift 2
			local estring `*'
			tokenize `"`estring'"', parse("||")
			local 1 `firstpart'
		}
		
		if `hasvtlist' {
			local vtilde : word `vnum' of `vtlist'
		}
		else {
			local vtilde : word 1 of `1'
			local vtilde Y`vnum'_`vtilde'
			local vtlist `vtlist' `vtilde'
		}
		
		// put in a yeq
		`qui' ddml yeq, learner(`vtilde') mname(`mname') vname(`varlist') : `1'
	
		if "`2'"~="|" & "`3'"~="|" {
			// done parsing
			local doparse = 0
		}

		mac shift 3
		
		local estring `*'
		
		local ++vnum
	}
	
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
	mata: `eqn'.shortstack = "`stack'"
	// model struct no longer needed
	mata: mata drop `mname'
	
	`qui' crossfit if `touse',		///
		ename(`eqn') noreplace		///
		`noisily'
	
	mata: st_local("nlearners", strofreal(`eqn'.nlearners))
	
	// stacking weights
	tempname ss_weights
	mata: st_matrix("`ss_weights'",return_result_item(`eqn',"`stack'","ss_weights", "1"))

	foreach var of varlist `vtilde' {
		qui replace `touse' = 0 if `var'==.
	}
	qui count if `touse'
	local N = r(N)
	
	// housekeeping
	// drop created predicted values
	foreach vtilde in `vtlist' {
		cap drop `vtilde'_1
	}
	if `dropstack' {
		drop `stack'_1
	}
	else {
		rename `stack'_1 `stack'
		label var `stack' "Stacking out-of-sample (fold) predictions"
	}
	
	ereturn clear
	ereturn post, depname(`varlist') obs(`N') esample(`touse')
	ereturn local cmd stacking
	ereturn local stack `stack'
	ereturn local predict stacking_p
	ereturn local eqn `eqn'
	ereturn local base_est `vtlist'
	ereturn scalar mcount = `nlearners'
	mat `ss_weights' = `ss_weights''
	ereturn matrix weights = `ss_weights'

end


// graph and/or table
program define stacking_graph_table, rclass
	version 16.0
	syntax ,								///
				[							///
					HOLDOUT1				/// vanilla option, abbreviates to "holdout"
					holdout(varname)		///
					GRAPH1					/// vanilla option, abbreviates to "graph"
					HISTogram				/// report histogram instead of default ROC
					goptions(string asis)	///
					lgoptions(string asis)	///
					table					/// 
				]

	// any graph options implies graph
	local graphflag = `"`graph1'`histogram'`goptions'`lgoptions'"'~=""
	
	if "`holdout'`holdout1'"=="" {
		local title In-sample
		tempvar touse
		qui gen `touse' = e(sample)
	}
	else {
		local title Out-of-sample
		// holdout variable provided, or default = not-in-sample?
		if "`holdout'"=="" {
			// default
			tempvar touse
			qui gen `touse' = 1-e(sample)
			// check number of OOS obs
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
		
	tempvar stacking_p stacking_r
	qui predict double `stacking_p'
	label var `stacking_p' "Prediction: Stacking Regressor"
	qui gen double `stacking_r' = `y' - `stacking_p'
	qui predict double `stacking_p', transform
	forvalues i=1/`nlearners' {
		local lname : word `i' of `learners'
		label var `stacking_p'`i' "Prediction: `lname'"
		tempvar stacking_r`i'
		qui gen double `stacking_r`i'' = `y' - `stacking_p'`i'
	}

	tempname g0
	if `graphflag' {
		twoway (scatter `stacking_p' `y') (line `y' `y') if `touse'		///
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
			twoway (scatter `stacking_p'`i' `y') (line `y' `y') if `touse'	///
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
		tempname m m_in m_out
		
		// column for in-sample MSPE
		qui sum `stacking_r' if e(sample)
		mat `m_in' = r(sd) * sqrt( (r(N)-1)/r(N) )
		forvalues i=1/`nlearners' {
			qui sum `stacking_r`i'' if e(sample)
			mat `m_in' = `m_in' \ (r(sd) * sqrt( (r(N)-1)/r(N) ))
		}
		
		// column for OOS MSPE
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
		
		mat `m' = `m_in' , `m_out'
		mat colnames `m' = MSPE_in MSPE_out
		mat rownames `m' = STACKING `learners'
		
		return matrix m = `m'

	}
	
end

